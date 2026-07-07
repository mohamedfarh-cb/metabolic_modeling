"""Per-model resource registry for multi-GEM support (PART A).

The pipeline can run against several base genome-scale models. EVERY model-specific
resource for the Human models comes from the Human-GEM repo checked out at a pinned
release tag (resources/Human-GEM@<tag>/) — never from the Recon2M.2 resources.

Backward compatibility: for ``Recon2M.2`` this returns EXACTLY the file paths the
pipeline hard-coded before, so the default run is byte-for-byte unchanged.

  | Resource            | Recon2M.2 (Entrez)   | Human1/Human2 (Ensembl, from Human-GEM@tag) |
  |---------------------|----------------------|---------------------------------------------|
  | base model (SBML)   | Recon2M.2 xml (ours) | Human-GEM xml at the tag                    |
  | Entrez<->Ensembl    | n/a                  | Human-GEM model/genes.tsv (geneEntrezID)    |
  | subsystems          | SBML r.subsystem     | Human-GEM model/reactions.tsv (subSystems)  |
  | tasks / medium      | our *_bigg / RPMI    | Human-GEM's own shipped task files          |
  | omics expression    | ours (Entrez)        | ours -> mapped to Ensembl via genes.tsv     |
"""
import csv
import glob
import logging
import os

MODEL_CHOICES = ['Recon2M.2', 'Human1', 'Human2']
MEM_CHOICES   = ['MBA', 'GIMME', 'FASTCORE', 'iMAT', 'INIT', 'tINIT', 'mCADRE']
DEFAULT_MODEL = 'Recon2M.2'
DEFAULT_MEM   = 'tINIT'

# MEM methods runnable in THIS environment (cobra 0.5.4 / py3.6, no MEM library).
# Everything else fails loud until the Phase-2 backend (troppo in a modern env) exists.
RUNNABLE_MEMS = {'tINIT'}

# Pinned Human-GEM release tags + commits (verified). Human1 = v1.0.0, Human2 = v2.0.0.
HUMAN_TAGS = {
    'Human1': {'tag': 'v1.0.0', 'commit': 'e4feea4943d807ac9e3435c10430750b76928405'},
    'Human2': {'tag': 'v2.0.0', 'commit': '635f533152dc5f7290ce04d12700eaa882273c3e'},
}
HUMAN_GEM_URL = 'https://github.com/SysBioChalmers/Human-GEM.git'

# Biomass/objective reaction per tag (verified via the inspection report; validated
# to exist in the model at pre-flight time — never used blindly).
HUMAN_BIOMASS = {'v1.0.0': 'biomass_components', 'v2.0.0': 'MAR04413'}

# Fail loudly if fewer than this fraction of the model's genes are covered by the
# (translated) omics matrix — a wrong/empty gene-ID namespace would sail through
# otherwise and silently produce a meaningless reconstruction.
MIN_GENE_COVERAGE = 0.10


def _pkg_dir():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # .../tf_metabolism


def _repo_root():
    return os.path.dirname(_pkg_dir())  # .../metabolic_modeling


def data_path(fname):
    return os.path.join(_pkg_dir(), 'data', fname)


def resources_dir():
    return os.path.join(_repo_root(), 'resources')


def human_gem_dir(name):
    return os.path.join(resources_dir(), 'Human-GEM@%s' % HUMAN_TAGS[name]['tag'])


def normalize_choices(values, allowed, what):
    """Split comma/space multi-values, validate against ``allowed`` (fail loud),
    de-duplicate preserving order."""
    tokens = []
    for v in values:
        tokens += [t.strip() for t in str(v).split(',') if t.strip()]
    bad = [t for t in tokens if t not in allowed]
    if bad:
        raise ValueError("Unknown --%s value(s): %s. Choose from: %s"
                         % (what, ', '.join(bad), ', '.join(allowed)))
    seen, out = set(), []
    for t in tokens:
        if t not in seen:
            seen.add(t); out.append(t)
    return out or [allowed[0]]


def check_mem_methods(mems):
    """Fail loud on any requested MEM that this environment cannot run. Only tINIT
    is wired; the rest need a MEM library (troppo / COBRA Toolbox) that is not
    installable here (cobra 0.5.4 / Python 3.6, no MATLAB). Never fall back."""
    unsupported = [m for m in mems if m not in RUNNABLE_MEMS]
    if unsupported:
        raise NotImplementedError(
            "MEM method(s) %s require a model-extraction library (troppo / COBRA "
            "Toolbox) that is NOT installed in this environment (cobra 0.5.4 / "
            "Python 3.6; no MATLAB). Only 'tINIT' is supported right now. This is a "
            "Phase-2 task (separate modern-Python env with troppo, bridged via "
            "subprocess) — see FUTURE_REQUESTS.md. Not falling back."
            % ', '.join(unsupported))


# ---------------------------------------------------------------------------
# Per-model configuration
# ---------------------------------------------------------------------------

def get_model_config(name):
    """Return a resource-config dict for one base model. Raises on missing files."""
    if name == 'Recon2M.2':
        return _recon2m2_config()
    if name in ('Human1', 'Human2'):
        return _human_gem_config(name)
    raise ValueError("Unknown model '%s'. Choose from: %s" % (name, ', '.join(MODEL_CHOICES)))


def _recon2m2_config():
    """EXACTLY the pipeline's historical Recon2M.2 resources (unchanged path)."""
    return {
        'name': 'Recon2M.2',
        'tag': None, 'commit': None,
        'model_file': data_path('Recon2M.2_Entrez_Gene_BIGG.xml'),
        'gene_namespace': 'entrez',
        'biomass_reaction': 'biomass_reaction',
        'medium_file': data_path('RPMI1640_medium.txt'),
        'metabolic_task_file': data_path('MetabolicTasks_bigg.csv'),
        'present_metabolite_file': data_path('essential_metabolites_bigg.txt'),
        'essential_reaction_file': data_path('essential_reactions.txt'),
        'getpra_file': data_path('GeTPRA.txt'),
        'subsystem_source': 'sbml',          # cobra exposes r.subsystem
        'genes_tsv': None,                   # no crosswalk needed (omics already Entrez)
        'reactions_tsv': None,
        'reconstruction_ready': True,
        'not_ready_reason': None,
    }


def _human_gem_config(name):
    tag = HUMAN_TAGS[name]['tag']
    root = human_gem_dir(name)
    if not os.path.isdir(root):
        raise FileNotFoundError(
            "Human-GEM not checked out for %s (tag %s) at %s. Clone it first:\n"
            "  git clone --branch %s %s %s"
            % (name, tag, root, tag, HUMAN_GEM_URL, root))

    # Layouts differ by tag: v2.0.0 = model/Human-GEM.xml + model/genes.tsv +
    # model/reactions.tsv; v1.0.0 = ModelFiles/xml/humanGEM.xml + ComplementaryData/.
    if os.path.isfile(os.path.join(root, 'model', 'Human-GEM.xml')):          # v2 layout
        model_file    = os.path.join(root, 'model', 'Human-GEM.xml')
        genes_tsv     = os.path.join(root, 'model', 'genes.tsv')
        reactions_tsv = os.path.join(root, 'model', 'reactions.tsv')
        tasks_dir     = os.path.join(root, 'data', 'metabolicTasks')
    else:                                                                     # v1 layout
        model_file    = os.path.join(root, 'ModelFiles', 'xml', 'humanGEM.xml')
        genes_tsv     = _first_existing(root, [
            'model/genes.tsv', 'ComplementaryData/Ensembl/ensembl_ID_mapping.tsv'])
        reactions_tsv = _first_existing(root, ['model/reactions.tsv'])
        tasks_dir     = os.path.join(root, 'ComplementaryData', 'metabolicTasks')

    if not os.path.isfile(model_file):
        raise FileNotFoundError("Human model SBML not found for %s: %s" % (name, model_file))

    task_file = None
    if os.path.isdir(tasks_dir):
        cand = (glob.glob(os.path.join(tasks_dir, 'metabolicTasks_Essential.*'))
                or glob.glob(os.path.join(tasks_dir, '*Essential*')))
        task_file = cand[0] if cand else None

    return {
        'name': name,
        'tag': tag, 'commit': HUMAN_TAGS[name]['commit'],
        'model_file': model_file,
        'gene_namespace': 'ensembl',
        'biomass_reaction': HUMAN_BIOMASS.get(tag),   # validated in-model at pre-flight
        'medium_file': None,                 # Human-GEM's own exchange naming (Phase-1.5)
        'metabolic_task_file': task_file,
        'present_metabolite_file': None,
        'essential_reaction_file': None,
        'getpra_file': None,
        # Human-GEM keeps subsystems in SBML <groups> (unreadable by cobra 0.5.4) and
        # they are NOT in reactions.tsv either -> unavailable in this environment.
        # Only matters post-reconstruction (enrichment), which Phase-1 gates anyway.
        'subsystem_source': None,
        'genes_tsv': genes_tsv if genes_tsv and os.path.isfile(genes_tsv) else None,
        'reactions_tsv': reactions_tsv,
        # Full task-driven tINIT on Human-GEM needs adapters not yet built (below).
        'reconstruction_ready': False,
        'not_ready_reason': (
            "task-driven tINIT reconstruction on Human-GEM is not yet wired: the "
            "Human task file format differs from MetabolicTasks_bigg.csv, and the "
            "biomass reaction / medium / present-metabolite / essential-reaction "
            "sets must be taken from Human-GEM's own naming (not the Recon2M.2 "
            "*_bigg files). Phase-1 runs load+inspection+gene-mapping and STOPS "
            "here rather than producing a wrong reconstruction. See FUTURE_REQUESTS.md."),
    }


def _first_existing(root, patterns):
    for p in patterns:
        hits = glob.glob(os.path.join(root, p))
        if hits:
            return hits[0]
    return None


# ---------------------------------------------------------------------------
# Entrez <-> Ensembl crosswalk (Human models) — in-repo genes.tsv only
# ---------------------------------------------------------------------------

def entrez_to_ensembl_map(cfg):
    """Build {Entrez -> Ensembl(ENSG)} from Human-GEM's own model/genes.tsv.
    No external annotation download. Raises if the crosswalk is unavailable."""
    genes_tsv = cfg.get('genes_tsv')
    if not genes_tsv or not os.path.isfile(genes_tsv):
        raise FileNotFoundError(
            "Entrez<->Ensembl crosswalk (model/genes.tsv) not found for %s. It is "
            "required to translate the Entrez omics matrix onto the model's Ensembl "
            "gene space." % cfg['name'])
    mapping = {}
    with open(genes_tsv, newline='', encoding='utf-8') as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith('#')]  # v1 file has a comment line
        reader = csv.DictReader(lines, delimiter='\t')
        cols = {c.strip().strip('"').lower(): c for c in (reader.fieldnames or [])}
        # v2 genes.tsv: genes / geneEntrezID ;  v1 ensembl_ID_mapping.tsv: Gene_stable_ID / NCBI_gene_ID
        ens_col = cols.get('genes') or cols.get('gene_stable_id')
        ent_col = cols.get('geneentrezid') or cols.get('ncbi_gene_id')
        if not ens_col or not ent_col:
            raise ValueError("crosswalk missing Ensembl/Entrez columns (%s): %s"
                             % (list(cols), genes_tsv))
        for row in reader:
            ens = str(row[ens_col]).strip().strip('"')
            for ent in str(row[ent_col]).strip().strip('"').split(';'):
                ent = ent.strip()
                if ent and ent.lower() not in ('na', 'nan', '') and ens:
                    mapping.setdefault(ent, ens)   # first mapping wins (stable)
    return mapping


def ensembl_to_symbol_map(cfg):
    """Build {Ensembl(ENSG) -> HGNC symbol} from Human-GEM's own model/genes.tsv
    (geneSymbols column). Returns {} for models without a genes.tsv (e.g. Recon2M.2,
    which already uses a bundled Entrez->symbol map). No external download."""
    genes_tsv = cfg.get('genes_tsv')
    mapping = {}
    if not genes_tsv or not os.path.isfile(genes_tsv):
        return mapping
    with open(genes_tsv, newline='', encoding='utf-8') as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith('#')]
        reader = csv.DictReader(lines, delimiter='\t')
        cols = {c.strip().strip('"').lower(): c for c in (reader.fieldnames or [])}
        ens_col = cols.get('genes') or cols.get('gene_stable_id')
        sym_col = cols.get('genesymbols') or cols.get('gene_symbol') or cols.get('symbol')
        if not ens_col or not sym_col:
            return mapping
        for row in reader:
            ens = str(row[ens_col]).strip().strip('"')
            sym = str(row[sym_col]).strip().strip('"').split(';')[0].strip()   # first symbol
            if ens and sym and sym.lower() not in ('na', 'nan', ''):
                mapping.setdefault(ens, sym)
    return mapping


def _reactions_tsv_subsystems(reactions_tsv):
    """Return the set of subsystem names from Human-GEM model/reactions.tsv."""
    if not reactions_tsv or not os.path.isfile(reactions_tsv):
        return set()
    subs = set()
    with open(reactions_tsv, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        col = next((c for c in (reader.fieldnames or [])
                    if 'subsystem' in c.strip().strip('"').lower()), None)
        if not col:
            return set()
        for row in reader:
            v = str(row[col]).strip().strip('"')
            if v and v.lower() != 'na':
                subs.add(v.split(';')[0].strip())
    return subs


def inspect_model(cfg):
    """Load the model SBML and print a short report BEFORE any reconstruction:
    reaction/gene/metabolite counts, gene-ID type, subsystem count, task
    availability, biomass reaction. Returns a stats dict."""
    import warnings
    from cobra.io import read_sbml_model
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = read_sbml_model(cfg['model_file'])

    gene_ids = [g.id for g in model.genes]
    if cfg['subsystem_source'] == 'sbml':
        subs = {r.subsystem.split(';')[0].strip() for r in model.reactions if r.subsystem}
        sub_src = 'SBML'
    elif cfg['subsystem_source']:
        subs = _reactions_tsv_subsystems(cfg['subsystem_source']); sub_src = 'reactions.tsv'
    else:
        subs = set(); sub_src = 'unavailable (SBML groups; cobra 0.5.4)'
    biomass = next((r.id for r in model.reactions
                    if 'biomass' in (r.id + ' ' + (r.name or '')).lower()), None)
    task_ok = bool(cfg.get('metabolic_task_file')) and os.path.isfile(cfg['metabolic_task_file'])

    stats = {
        'name': cfg['name'], 'tag': cfg['tag'], 'commit': cfg['commit'],
        'n_reactions': len(model.reactions), 'n_genes': len(model.genes),
        'n_metabolites': len(model.metabolites),
        'gene_namespace': cfg['gene_namespace'],
        'sample_gene_ids': gene_ids[:3],
        'n_subsystems': len(subs), 'biomass_reaction': biomass,
        'task_file': cfg.get('metabolic_task_file'), 'task_available': task_ok,
    }
    logging.info('--- Model inspection: %s%s ---', cfg['name'],
                 (' (%s @ %s)' % (cfg['tag'], (cfg['commit'] or '')[:7])) if cfg['tag'] else '')
    logging.info('  reactions=%d  genes=%d  metabolites=%d',
                 stats['n_reactions'], stats['n_genes'], stats['n_metabolites'])
    logging.info('  gene-ID type=%s  e.g. %s', cfg['gene_namespace'], stats['sample_gene_ids'])
    logging.info('  subsystems=%d (source=%s)', stats['n_subsystems'], sub_src)
    logging.info('  biomass reaction=%s', biomass)
    logging.info('  metabolic tasks=%s', 'available (%s)' % os.path.basename(cfg['metabolic_task_file'])
                 if task_ok else 'NOT FOUND')
    return stats


def translate_omics_to_model(cfg, omics_file, out_file):
    """Translate an Entrez-indexed omics matrix onto the model's gene namespace.
    Recon2M.2 (Entrez): pass through unchanged. Human (Ensembl): map via the
    in-repo genes.tsv crosswalk. Reports coverage and FAILS LOUD if too few of the
    model's genes are covered. Returns a coverage dict; ``out_file`` is written for
    the Human case (or equals ``omics_file`` for Recon2M.2)."""
    import warnings
    import pandas as pd
    from cobra.io import read_sbml_model

    if cfg['gene_namespace'] == 'entrez':
        return {'model': cfg['name'], 'translated': False, 'omics_file': omics_file}

    omics = pd.read_csv(omics_file, index_col=0)
    omics.index = omics.index.map(lambda x: str(x).strip())
    mapping = entrez_to_ensembl_map(cfg)                       # Entrez -> ENSG
    n_omics = len(omics.index)
    mapped_mask = omics.index.to_series().isin(mapping)
    omics = omics.loc[mapped_mask.values]
    omics.index = [mapping[i] for i in omics.index]
    omics = omics.groupby(level=0).mean()                     # collapse Entrez->ENSG collisions

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = read_sbml_model(cfg['model_file'])
    model_genes = {g.id for g in model.genes}
    covered = model_genes & set(omics.index)
    coverage = float(len(covered)) / float(len(model_genes) or 1)

    cov = {'model': cfg['name'], 'translated': True, 'omics_file': out_file,
           'n_omics_genes': n_omics, 'n_mapped': int(mapped_mask.sum()),
           'n_model_genes': len(model_genes), 'n_model_covered': len(covered),
           'model_gene_coverage': coverage}
    logging.info('  gene-ID mapping (%s): omics genes=%d, Entrez->Ensembl mapped=%d; '
                 'model genes covered=%d/%d (%.1f%%)', cfg['name'], n_omics,
                 int(mapped_mask.sum()), len(covered), len(model_genes), 100 * coverage)
    if coverage < MIN_GENE_COVERAGE:
        raise ValueError(
            "Gene-ID mapping coverage for %s is only %.1f%% of model genes (< %.0f%% "
            "minimum). The Entrez omics IDs are not landing on the model's Ensembl "
            "namespace — refusing to run a meaningless reconstruction."
            % (cfg['name'], 100 * coverage, 100 * MIN_GENE_COVERAGE))
    omics.to_csv(out_file)
    return cov


# ---------------------------------------------------------------------------
# Human-GEM reconstruction adapters (Phase-1b) — resource-layer only.
#
# These make the pinned Human tag's OWN resources available to the tINIT path:
#   (1) metabolic tasks (RAVEN IN/OUT/EQU format, NOT the Recon2M.2 _bigg CSV),
#   (2) biomass/objective reaction, (3) exchange/boundary reactions, and
#   (4) subsystems from the SBML <groups> block (unreadable via cobra 0.5.4).
# They VALIDATE the resources; the actual RAVEN->internal task conversion + full
# reconstruction is a separate, validated session (see FUTURE_REQUESTS.md [P1]).
# ---------------------------------------------------------------------------

def parse_human_tasks(task_file):
    """Parse a Human-GEM RAVEN essential-task file into a list of task dicts.
    Skips '#' comment/section rows. Raises if the file is absent or has 0 tasks."""
    if not task_file or not os.path.isfile(task_file):
        raise FileNotFoundError("Human metabolic-task file not found: %s" % task_file)
    tasks = []
    with open(task_file, newline='', encoding='utf-8') as fh:   # task descriptions have non-ASCII
        reader = csv.DictReader(fh, delimiter='\t')
        lead = (reader.fieldnames or [''])[0]                  # leading col ('' ; '#' = comment)
        for row in reader:
            if str(row.get(lead, '')).strip() == '#':
                continue
            tid = str(row.get('ID', '')).strip()
            if not tid:
                continue
            rec = {'ID': tid, 'DESCRIPTION': str(row.get('DESCRIPTION', '')).strip(),
                   'IN': str(row.get('IN', '')).strip(),
                   'OUT': str(row.get('OUT', '')).strip(),
                   'EQU': str(row.get('EQU', '')).strip(),
                   'SHOULD_FAIL': str(row.get('SHOULD FAIL', '')).strip()}
            for col in ('IN LB', 'IN UB', 'OUT LB', 'OUT UB', 'EQU LB', 'EQU UB'):
                rec[col] = str(row.get(col, '') or '').strip()
            tasks.append(rec)
    if not tasks:
        raise ValueError("Human metabolic-task file parsed to 0 tasks: %s" % task_file)
    return tasks


def human_subsystem_map(cfg):
    """Build {reaction_id -> subsystem} from the SBML <groups> block (cobra 0.5.4
    does not expose it). Returns {} if the model has no groups."""
    import re
    xml = cfg['model_file']
    if not os.path.isfile(xml):
        return {}
    with open(xml, encoding='utf-8', errors='replace') as fh:
        text = fh.read()
    sub_map = {}
    for gm in re.finditer(r'<groups:group\b[^>]*groups:name="([^"]*)"(.*?)</groups:group>',
                          text, re.S):
        name = gm.group(1)
        for mm in re.finditer(r'groups:idRef="([^"]+)"', gm.group(2)):
            rid = mm.group(1)
            if rid.startswith('R_'):            # cobra 0.5.4 strips the SBML 'R_' prefix (v1 layout)
                rid = rid[2:]
            sub_map[rid] = name
    return sub_map


def attach_subsystems(model, sub_map):
    """Set r.subsystem on a loaded model from a {reaction_id -> subsystem} map, so
    downstream enrichment / gene_subsystem_flux (which read r.subsystem) work for
    Human models. No-op for reactions not in the map. Returns coverage fraction."""
    n = 0
    for r in model.reactions:
        s = sub_map.get(r.id)
        if s:
            r.subsystem = s; n += 1
    return float(n) / float(len(model.reactions) or 1)


def resolve_exchanges(model):
    """Exchange/boundary reaction IDs in the model's OWN namespace (Human = MAR####)."""
    return [r.id for r in model.reactions if r.boundary]


def preflight_reconstruction(cfg, output_dir):
    """Validate the four Human reconstruction adapters BEFORE any Gurobi run and
    write a pre-flight report. Blocking adapters (tasks, biomass, exchanges) FAIL
    LOUD if they cannot resolve; subsystems degrade gracefully (warn only).

    Returns (ready, report). ``ready`` is True when all blocking adapters resolve —
    but the caller still keeps the Gurobi reconstruction gated behind the separate
    RAVEN->internal task conversion ([P1]); this only proves the resources exist."""
    import warnings
    from cobra.io import read_sbml_model
    os.makedirs(output_dir, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = read_sbml_model(cfg['model_file'])
    rxn_ids = {r.id for r in model.reactions}
    fails = []

    # (1) metabolic tasks (blocking)
    try:
        tasks = parse_human_tasks(cfg['metabolic_task_file'])
        n_tasks = len(tasks)
    except (FileNotFoundError, ValueError) as e:
        n_tasks = 0; fails.append('tasks: %s' % e)

    # (2) biomass / objective reaction (blocking)
    biomass = cfg.get('biomass_reaction')
    biomass_ok = bool(biomass) and biomass in rxn_ids
    if not biomass_ok:
        fails.append("biomass: reaction '%s' not found in %s model" % (biomass, cfg['name']))

    # (3) exchange / boundary reactions (blocking)
    exchanges = resolve_exchanges(model)
    if not exchanges:
        fails.append('exchanges: no boundary reactions resolved')

    # (4) subsystems (non-blocking — enrichment degrades gracefully if absent)
    sub_map = human_subsystem_map(cfg)
    sub_cov = float(len(rxn_ids & set(sub_map))) / float(len(rxn_ids) or 1)

    # (5) growth feasibility (blocking): biomass must carry flux under the medium,
    #     else run_INIT (which forces biomass) is infeasible.
    bio_flux, blocked = biomass_feasibility(cfg, model)
    if bio_flux <= 1e-9:
        fails.append('biomass infeasible under medium (flux=%.4g); blocked precursors: %s'
                     % (bio_flux, ', '.join('%s(%s)' % b for b in blocked[:8]) or 'none'))

    ready = not fails
    report = {'model': cfg['name'], 'tag': cfg['tag'], 'commit': cfg['commit'],
              'n_tasks': n_tasks, 'biomass': biomass, 'biomass_ok': biomass_ok,
              'n_exchanges': len(exchanges), 'subsystem_coverage': sub_cov,
              'n_subsystems': len(set(sub_map.values())), 'blocking_failures': fails,
              'ready': ready}

    logging.info('--- Reconstruction pre-flight: %s (%s) ---', cfg['name'], cfg['tag'])
    logging.info('  metabolic tasks parsed: %d (%s)', n_tasks,
                 os.path.basename(cfg['metabolic_task_file']) if cfg.get('metabolic_task_file') else 'NONE')
    logging.info('  biomass reaction: %s  (in model: %s)', biomass, biomass_ok)
    logging.info('  exchange reactions: %d', len(exchanges))
    logging.info('  biomass growth: MAR flux=%.3g (%s)', bio_flux,
                 'FEASIBLE' if bio_flux > 1e-9 else 'INFEASIBLE')
    if sub_map:
        logging.info('  subsystems: %d, covering %.1f%% of reactions',
                     len(set(sub_map.values())), 100 * sub_cov)
    else:
        logging.warning('  subsystems unavailable for %s — enrichment will degrade gracefully',
                        cfg['name'])
    for f in fails:
        logging.error('  BLOCKING adapter failure -> %s', f)

    with open(os.path.join(output_dir, 'RECONSTRUCTION_PREFLIGHT.txt'), 'w') as fp:
        fp.write('model=%s  tag=%s  commit=%s\n' % (cfg['name'], cfg['tag'], cfg['commit']))
        fp.write('metabolic_tasks_parsed=%d\n' % n_tasks)
        fp.write('biomass=%s  in_model=%s\n' % (biomass, biomass_ok))
        fp.write('exchange_reactions=%d\n' % len(exchanges))
        fp.write('biomass_flux=%.4g  growth_feasible=%s\n' % (bio_flux, bio_flux > 1e-9))
        fp.write('subsystems=%d  reaction_coverage=%.3f\n' % (len(set(sub_map.values())), sub_cov))
        fp.write('adapters_ready=%s\n' % ready)
        if fails:
            fp.write('blocking_failures:\n' + ''.join('  - %s\n' % f for f in fails))
        fp.write('\nNOTE: adapters validated at the RESOURCE level. The task-driven '
                 'Gurobi reconstruction stays gated pending the RAVEN->internal task + '
                 'MAM####/MAR#### namespace conversion (FUTURE_REQUESTS.md [P1]).\n')
    return ready, report


# ---------------------------------------------------------------------------
# RAVEN -> internal task conversion + validation gate (Phase-1c)
#
# The internal task format (metabolic_task.py: Type/ID/Medium/Constraints) is a
# SINGLE-metabolite producibility check — "metabolite X is producible from medium
# M". Human-GEM RAVEN tasks are richer: EQU (forced reactions added to the model)
# and multiple SIMULTANEOUS OUT metabolites have NO field in the internal format.
# Policy: never silently drop/approximate a task -> such tasks make the conversion
# FAIL, and the reconstruction stays gated. (See FUTURE_REQUESTS.md [P1].)
# ---------------------------------------------------------------------------

def _model_name_resolver(model):
    """Return f(token='name[comp]') -> metabolite id, from the loaded model's
    metabolite names + compartments (case-insensitive)."""
    import re
    name_map = {}
    for m in model.metabolites:
        comp = getattr(m, 'compartment', '') or (m.id[-1] if m.id else '')
        nm = (m.name or '').strip().lower()
        if nm:
            name_map.setdefault((nm, comp), m.id)

    def resolve(token):
        mt = re.match(r'^(.*)\[(\w+)\]$', token.strip())
        if not mt:
            return None
        return name_map.get((mt.group(1).strip().lower(), mt.group(2)))
    return resolve


def _parse_equ(equ_str, resolve):
    """Parse a RAVEN EQU string 'A[c] + 2 B[c] => C[c]' -> (stoich{met_id:coeff},
    reversible). Fails loud on any unresolved metabolite."""
    import re
    s = equ_str.strip()
    if '<=>' in s:
        lhs, rhs = s.split('<=>', 1); reversible = True
    elif '=>' in s:
        lhs, rhs = s.split('=>', 1); reversible = False
    else:
        raise ValueError("cannot parse EQU (no arrow): %r" % equ_str)
    stoich = {}

    def add(side, sign):
        # split on ' + ' (whitespace-delimited) — metabolite names may contain '+' (e.g. H+)
        for term in re.split(r'\s+\+\s+', side.strip()):
            term = term.strip()
            if not term:
                continue
            mt = re.match(r'^(?:([0-9.]+)\s+)?(.+\[\w+\])$', term)
            if not mt:
                raise ValueError("cannot parse EQU term: %r" % term)
            coeff = float(mt.group(1)) if mt.group(1) else 1.0
            mid = resolve(mt.group(2))
            if mid is None:
                raise KeyError("unresolved EQU metabolite: %s" % mt.group(2))
            stoich[mid] = stoich.get(mid, 0.0) + sign * coeff
    add(lhs, -1.0)
    add(rhs, +1.0)
    return stoich, reversible


# Base inorganic medium always available during RAVEN/Recon task checks (matches
# the ion set in MetabolicTasks_bigg.csv's Medium string). The task's own IN adds
# the specific carbon/N substrate under test; these inorganics are never the point
# of a task but are needed for it to balance (e.g. O2 for beta-oxidation/OXPHOS).
HUMAN_BASE_MEDIUM = (
    # inorganic ions / gases (as in MetabolicTasks_bigg.csv's Medium)
    ['O2[e]', 'H2O[e]', 'CO2[e]', 'H+[e]', 'Pi[e]', 'Na+[e]', 'K+[e]', 'chloride[e]',
     'Ca2+[e]', 'Fe2+[e]', 'Fe3+[e]', 'sulfate[e]', 'NH3[e]', 'HCO3-[e]', 'H2S[e]']
    # + the 20 proteinogenic amino acids: standard building blocks a task rich medium
    #   supplies (e.g. the essential 'Protein synthesis from AAs' task needs them, since
    #   humans cannot make the essential AAs de novo).
    + ['%s[e]' % a for a in
       ('alanine', 'arginine', 'asparagine', 'aspartate', 'cysteine', 'glutamate',
        'glutamine', 'glycine', 'histidine', 'isoleucine', 'leucine', 'lysine',
        'methionine', 'phenylalanine', 'proline', 'serine', 'threonine',
        'tryptophan', 'tyrosine', 'valine')])


def evaluate_human_tasks(cfg, model=None):
    """Faithful RAVEN checkTasks evaluation of the Human-GEM essential tasks on a
    (full or context) model. Per task: close all existing boundaries; open uptake
    for IN metabolites (ALLMETSIN = open all); force simultaneous production of all
    OUT metabolites; inject EQU reactions with their bounds; solve feasibility.
    ESSENTIAL task PASSES iff feasible; SHOULD-FAIL PASSES iff INFEASIBLE. Fails
    loud on any unresolved metabolite. Returns a list of per-task result dicts."""
    import warnings
    from cobra import Reaction
    from cobra.io import read_sbml_model
    FORCE = 1e-3                    # small positive lower bound to force production/flux

    if model is None:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            model = read_sbml_model(cfg['model_file'])
    resolve = _model_name_resolver(model)
    tasks = parse_human_tasks(cfg['metabolic_task_file'])

    boundary = [r for r in model.reactions if r.boundary]
    orig = {r.id: (r.lower_bound, r.upper_bound) for r in boundary}

    def num(v, default):
        v = str(v).strip()
        return float(v) if v else default

    results = []
    for t in tasks:
        should_fail = t.get('SHOULD_FAIL', '').strip() in ('1', 'true', 'True')
        added = []
        try:
            # RAVEN checkTasks: block all UPTAKE but allow free secretion of
            # byproducts; the only inputs are the task's IN metabolites (below).
            for r in boundary:
                r.lower_bound = 0.0; r.upper_bound = 1000.0
            ins = [x.strip() for x in t['IN'].split(';') if x.strip()]
            outs = [x.strip() for x in t['OUT'].split(';') if x.strip()]
            if any(x.upper().startswith('ALLMETSIN') for x in ins):
                for r in boundary:               # ALLMETSIN in IN: open all uptake
                    r.lower_bound = -1000.0
                ins = [x for x in ins if not x.upper().startswith('ALLMETSIN')]
            if any(x.upper().startswith('ALLMETSIN') for x in outs):
                for r in boundary:               # ALLMETSIN in OUT: open all secretion
                    r.upper_bound = 1000.0
                outs = [x for x in outs if not x.upper().startswith('ALLMETSIN')]

            in_lb, in_ub = num(t.get('IN LB'), 0.0), num(t.get('IN UB'), 1000.0)
            out_lb, out_ub = num(t.get('OUT LB'), 0.0), num(t.get('OUT UB'), 1000.0)
            equ_lb, equ_ub = num(t.get('EQU LB'), 0.0), num(t.get('EQU UB'), 1000.0)

            for i, tok in enumerate(ins):        # uptake: '-> M'  (produce = enter)
                mid = resolve(tok)
                if mid is None:
                    raise KeyError("unresolved IN metabolite: %s" % tok)
                rx = Reaction('TASKIN_%d' % i)
                rx.lower_bound, rx.upper_bound = in_lb, in_ub
                rx.add_metabolites({model.metabolites.get_by_id(mid): 1.0}); added.append(rx)
            for i, tok in enumerate(HUMAN_BASE_MEDIUM):   # base inorganic medium (uptake allowed)
                mid = resolve(tok)
                if mid is None:
                    continue
                rx = Reaction('TASKMED_%d' % i)
                rx.lower_bound, rx.upper_bound = 0.0, 1000.0
                rx.add_metabolites({model.metabolites.get_by_id(mid): 1.0}); added.append(rx)
            for i, tok in enumerate(outs):       # demand: 'M ->' forced (simultaneous)
                mid = resolve(tok)
                if mid is None:
                    raise KeyError("unresolved OUT metabolite: %s" % tok)
                rx = Reaction('TASKOUT_%d' % i)
                rx.lower_bound = max(out_lb, FORCE); rx.upper_bound = out_ub
                rx.add_metabolites({model.metabolites.get_by_id(mid): -1.0}); added.append(rx)
            if t.get('EQU'):
                for j, equ in enumerate([e for e in t['EQU'].split(';') if e.strip()]):
                    stoich, rev = _parse_equ(equ, resolve)
                    rx = Reaction('TASKEQU_%d' % j)
                    rx.lower_bound = (-1000.0 if rev else 0.0) if equ_lb == 0.0 else equ_lb
                    rx.upper_bound = equ_ub
                    rx.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoich.items()})
                    added.append(rx)

            model.add_reactions(added)
            model.objective = added[-1] if added else None
            sol = model.optimize(solver='gurobi')
            feasible = (sol.status == 'optimal')
            passed = (not feasible) if should_fail else feasible
            results.append({'ID': t['ID'], 'desc': t['DESCRIPTION'][:40],
                            'should_fail': should_fail, 'feasible': feasible,
                            'result': 'PASSED' if passed else 'FAILED'})
        finally:
            if added:
                model.remove_reactions(added, remove_orphans=False)
            for r in boundary:                   # restore original boundary bounds
                r.lower_bound, r.upper_bound = orig[r.id]
    return results


def apply_human_medium(model, cfg):
    """Human growth medium adapter. Human-GEM ships its generic model with ALL 1660
    exchange reactions OPEN for uptake (its default medium); the only thing needed
    for growth is to OPEN the biomass reaction, which ships bounded to [0, 0] (the
    standard GEM convention of disabling the objective by default). No ad-hoc medium
    and no gap-filling. Returns the biomass reaction id. Fails loud if absent."""
    bid = cfg.get('biomass_reaction')
    if not bid or bid not in [r.id for r in model.reactions]:
        raise ValueError("biomass reaction %r not found in %s model" % (bid, cfg['name']))
    bio = model.reactions.get_by_id(bid)
    if bio.upper_bound <= 0.0:
        bio.upper_bound = 1000.0                 # enable the objective (ships closed)
    return bid


def biomass_feasibility(cfg, model=None):
    """Full-model growth oracle: apply the medium, optimize biomass, return
    (flux, blocked_precursors). flux > 0 => growth-feasible. If 0, blocked_precursors
    lists the biomass reactant metabolites that cannot be produced (for the report)."""
    import warnings
    from cobra import Reaction
    from cobra.io import read_sbml_model
    if model is None:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            model = read_sbml_model(cfg['model_file'])
    bid = apply_human_medium(model, cfg)
    bio = model.reactions.get_by_id(bid)
    model.objective = bio
    flux = model.optimize(solver='gurobi').f or 0.0
    blocked = []
    if flux <= 1e-9:
        for met in bio.reactants:
            dm = Reaction('DMprobe'); dm.lower_bound = 0.0; dm.upper_bound = 1000.0
            dm.add_metabolites({met: -1.0}); model.add_reactions([dm]); model.objective = dm
            v = model.optimize(solver='gurobi').f or 0.0
            model.remove_reactions([dm], remove_orphans=False)
            if v <= 1e-9:
                blocked.append((met.id, met.name))
    return flux, blocked


def convert_raven_tasks(cfg):
    """Classify each Human-GEM RAVEN essential task by whether it is faithfully
    representable in the internal single-metabolite format, and (for representable
    ones) resolve every metabolite token into the model's MAM####<comp> namespace.
    Returns a report dict. Raises via parse_human_tasks if no task file."""
    import warnings
    from cobra.io import read_sbml_model
    tasks = parse_human_tasks(cfg['metabolic_task_file'])

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = read_sbml_model(cfg['model_file'])
    # {(normalized_name, compartment) -> metabolite id} from the loaded model
    name_map = {}
    for m in model.metabolites:
        comp = getattr(m, 'compartment', '') or (m.id[-1] if m.id else '')
        nm = (m.name or '').strip().lower()
        if nm:
            name_map.setdefault((nm, comp), m.id)

    import re
    def _resolve(tok):
        tok = tok.strip()
        mt = re.match(r'^(.*)\[(\w+)\]$', tok)
        if not mt:
            return None
        nm, comp = mt.group(1).strip().lower(), mt.group(2)
        return name_map.get((nm, comp))

    representable, unrepresentable, unresolved = [], [], []
    for t in tasks:
        reasons = []
        if t.get('EQU'):
            reasons.append('EQU (forced reaction — no internal field)')
        outs = [o for o in t['OUT'].split(';') if o.strip()]
        if len(outs) > 1:
            reasons.append('multiple simultaneous OUT (%d — internal checks one target)' % len(outs))
        if reasons:
            unrepresentable.append((t['ID'], t['DESCRIPTION'][:40], reasons))
            continue
        # representable: resolve every IN/OUT token into MAM####<comp>
        toks = [x for x in (t['IN'].split(';') + t['OUT'].split(';')) if x.strip()]
        bad = [x.strip() for x in toks if _resolve(x) is None]
        if bad:
            unresolved += [(t['ID'], b) for b in bad]
        representable.append(t['ID'])

    passed = (not unrepresentable) and (not unresolved)
    return {'model': cfg['name'], 'n_tasks': len(tasks),
            'n_representable': len(representable), 'unrepresentable': unrepresentable,
            'unresolved_tokens': unresolved, 'passed': passed}


def validate_raven_conversion(cfg, output_dir):
    """The cheap correctness gate (STEP 1/2/3b), run on the FULL model BEFORE any
    per-patient reconstruction:
      (1/2) every RAVEN task converts + every metabolite/EQU token resolves into the
            model's MAM####<comp> namespace (``evaluate_human_tasks`` raises loud on
            any unresolved token), and
      (3b)  all ESSENTIAL tasks PASS on the complete model and all SHOULD-FAIL tasks
            FAIL — the oracle for a correct port.
    Writes RAVEN_CONVERSION_VALIDATION.txt (per-task table + overall). Returns
    (passed, report). Only overall PASS un-gates the reconstruction."""
    os.makedirs(output_dir, exist_ok=True)
    report_path = os.path.join(output_dir, 'RAVEN_CONVERSION_VALIDATION.txt')
    try:
        results = evaluate_human_tasks(cfg)          # resolves + evaluates on FULL model
    except (KeyError, ValueError) as e:              # unresolved token / unparseable EQU
        logging.error('  RAVEN conversion FAILED to resolve: %s', e)
        with open(report_path, 'w') as fp:
            fp.write('model=%s  tag=%s\noverall=FAIL (unresolved/unparseable)\n%s\n'
                     % (cfg['name'], cfg['tag'], e))
        return False, {'passed': False, 'error': str(e)}

    failed = [r for r in results if r['result'] == 'FAILED']
    passed = not failed
    n_sf = sum(1 for r in results if r['should_fail'])
    logging.info('--- RAVEN checkTasks on FULL %s model: %d tasks -> %d PASSED, %d FAILED ---',
                 cfg['name'], len(results), len(results) - len(failed), len(failed))
    for r in failed:
        logging.error('  FAILED: %s (%s) feasible=%s should_fail=%s',
                      r['ID'], r['desc'], r['feasible'], r['should_fail'])
    with open(report_path, 'w') as fp:
        fp.write('model=%s  tag=%s\n' % (cfg['name'], cfg['tag']))
        fp.write('essential_tasks=%d  should_fail=%d  passed=%d  failed=%d\n'
                 % (len(results), n_sf, len(results) - len(failed), len(failed)))
        fp.write('overall=%s\n\n' % ('PASS' if passed else 'FAIL'))
        fp.write('per-task (full-model oracle):\n')
        for r in results:
            fp.write('  [%s] %-8s %s (feasible=%s, should_fail=%s)\n'
                     % (r['result'], r['ID'], r['desc'], r['feasible'], r['should_fail']))
        if not passed:
            fp.write('\nAn ESSENTIAL task failing on the FULL model means the RAVEN->internal '
                     'port is WRONG (or an input is under-specified). Reconstruction stays '
                     'gated — do not un-gate. See FUTURE_REQUESTS.md [P1].\n')
    return passed, {'passed': passed, 'n_tasks': len(results), 'n_failed': len(failed)}
