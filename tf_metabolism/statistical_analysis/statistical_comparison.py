import warnings
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, ttest_ind, ttest_rel, ranksums, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

# Suppress scipy precision-loss warnings caused by near-identical data
warnings.filterwarnings('ignore', category=RuntimeWarning,
                        message='Precision loss occurred in moment calculation')


def _rank_biserial_corr(pop1, pop2):
    """Compute rank-biserial correlation from Mann-Whitney U statistic.

    r_rb = 1 - 2U / (n1 * n2)
    Range: [-1, +1]. |r_rb| > 0.5 corresponds to a 'large' effect size.
    Returns (U, rbc). Returns (nan, nan) on error.
    """
    n1, n2 = len(pop1), len(pop2)
    if n1 == 0 or n2 == 0:
        return np.nan, np.nan
    try:
        res = mannwhitneyu(pop1, pop2, alternative='two-sided')
        U   = res.statistic
        rbc = 1.0 - 2.0 * U / (n1 * n2)
        return U, rbc
    except Exception:
        return np.nan, np.nan


def _cohens_d(pop1, pop2):
    """Cohen's d: pooled-SD standardized mean difference (condition2 - condition1).

    d = (mean2 - mean1) / s_pooled,
    s_pooled = sqrt(((n1-1) s1^2 + (n2-1) s2^2) / (n1 + n2 - 2))  (sample SD, ddof=1)
    Range: unbounded; |d| ~ 0.2 small, 0.5 medium, 0.8 large. This is the
    parametric (t-test) analogue of the rank-biserial correlation used on the
    non-parametric branch. Returns nan if undefined (a group with n<2, or zero
    pooled SD).
    """
    a = np.asarray(pop1, dtype=float)
    b = np.asarray(pop2, dtype=float)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan
    s1, s2 = np.std(a, ddof=1), np.std(b, ddof=1)
    s_pooled = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    if s_pooled == 0 or np.isnan(s_pooled):
        return np.nan
    return (np.mean(b) - np.mean(a)) / s_pooled


def two_grouped_data_comparison(condition1_df, condition2_df,
                                related=False,
                                p_value_cutoff=0.05,
                                filter_by='FDR',
                                force_non_parametric=False,
                                rbc_cutoff=None,
                                cohens_d_cutoff=None):
    """Compare two groups of flux/omics profiles.

    Parameters
    ----------
    condition1_df, condition2_df : pd.DataFrame
        Rows = features (reactions/genes), columns = samples.
    related : bool
        Paired test if True.
    p_value_cutoff : float
        Significance threshold (default 0.05).
    filter_by : str
        'FDR'   — Benjamini-Hochberg adjusted p < p_value_cutoff  (default)
        'raw_p' — raw p < p_value_cutoff  (paper Section 2.10 for flux DF)
        'Bonferroni' — Bonferroni-corrected p
        'None'  — return all results without filtering
    force_non_parametric : bool
        Force Wilcoxon rank-sum test regardless of sample size.
    rbc_cutoff : float or None
        If set (e.g. 0.5), additionally require |rank-biserial corr| > rbc_cutoff.
        Only applied for the Wilcoxon rank-sum test (non-parametric, unrelated).
        Paper (Section 2.10) uses rbc_cutoff=0.5.
    cohens_d_cutoff : float or None
        If set (e.g. 0.8 = conventional "large"), additionally require
        |Cohen's d| > cohens_d_cutoff on the parametric (t-test) branch.

        Effect size is decoupled from the test choice and ALWAYS computed: the
        non-parametric branch reports rank-biserial correlation and the
        parametric branch reports Cohen's d (both also stored in a generic
        'Effect size' column tagged by 'Effect size type'). The two cutoffs are
        kept separate because the measures live on different scales
        (rank-biserial in [-1, 1]; Cohen's d unbounded). This keeps an
        effect-size gate active even when n>=30 switches the code to the t-test.
        Both default to None, which preserves the original behavior exactly.
    """
    while filter_by not in ['FDR', 'Bonferroni', 'None', 'raw_p']:
        print('Input argument filter_by is not correct.')
        filter_by = input("Choose ['FDR'/'raw_p'/'Bonferroni'/'None']: ")

    pure_condition1_df = condition1_df
    pure_condition2_df = condition2_df

    # --- select test method ---
    if force_non_parametric:
        auto_method = 'non_parametric'
    elif (len(list(pure_condition1_df.columns)) < 30 or
          len(list(pure_condition2_df.columns)) < 30):
        auto_method = 'non_parametric'
    else:
        auto_method = 'parametric'

    use_wilcoxon_ranksum = (auto_method == 'non_parametric' and not related)

    if auto_method == 'non_parametric':
        if related:
            method = 'Wilcoxon signed-rank test'
            target_function = wilcoxon
        else:
            method = 'Wilcoxon rank-sum test'
            target_function = ranksums        # p-value source
    else:
        if related:
            method = 'T-test on two paired samples'
            target_function = ttest_rel
        else:
            method = 'T-test on two independent samples'
            target_function = ttest_ind

    # --- build result columns ---
    # Effect size is ALWAYS computed and is decoupled from the test choice: the
    # Wilcoxon rank-sum branch reports rank-biserial correlation; the t-test
    # branch reports Cohen's d. A generic 'Effect size' (+ 'Effect size type')
    # column records whichever applies, so the effect-size gate is never lost
    # when n>=30 switches the test from Wilcoxon to t-test.
    if use_wilcoxon_ranksum:
        compute_effect, effect_type, effect_specific_col = True, 'rank_biserial', 'Rank-biserial corr'
    elif auto_method == 'parametric':
        compute_effect, effect_type, effect_specific_col = True, 'cohens_d', 'Cohen d'
    else:
        # related (paired) non-parametric: no effect size (preserves old behavior)
        compute_effect, effect_type, effect_specific_col = False, None, None

    base_cols = ['Condition1 mean', 'Condition1 std',
                 'Condition2 mean', 'Condition2 std',
                 'Method', 'P', 'log2 (condition2/condition1)', 'Status']
    if compute_effect:
        p_pos = base_cols.index('P') + 1
        base_cols[p_pos:p_pos] = [effect_specific_col, 'Effect size', 'Effect size type']

    result_df_dict = {c: [] for c in base_cols}
    result_df_index = []

    pure_index_list = list(set(pure_condition1_df.index) & set(pure_condition2_df.index))

    for each_index in pure_index_list:
        pop1 = pure_condition1_df.loc[each_index, :].dropna().values.tolist()
        pop2 = pure_condition2_df.loc[each_index, :].dropna().values.tolist()
        if len(pop1) == 0 or len(pop2) == 0:
            continue

        error = False
        effect_val = np.nan
        try:
            if use_wilcoxon_ranksum:
                # Use mannwhitneyu for both p-value and rank-biserial (same function)
                res = mannwhitneyu(pop1, pop2, alternative='two-sided')
                pvalue = res.pvalue
                U = res.statistic
                effect_val = 1.0 - 2.0 * U / (len(pop1) * len(pop2))
            else:
                _, pvalue = target_function(pop1, pop2)
                if effect_type == 'cohens_d':
                    effect_val = _cohens_d(pop1, pop2)
        except Exception:
            error = True
            pvalue = np.nan

        if np.isnan(pvalue) or error:
            continue

        mean1 = np.mean(pop1)
        mean2 = np.mean(pop2)

        if mean1 * mean2 < 0.0:
            status = 'Reversed'
        elif mean1 == mean2:
            status = 'Unchanged'
        elif abs(mean2) > abs(mean1):
            status = 'UP'
        else:
            status = 'DOWN'

        if mean1 == 0 or mean2 / mean1 <= 0:
            log2FC = np.nan
        else:
            log2FC = np.log2(mean2 / mean1)

        row = [mean1, np.std(pop1), mean2, np.std(pop2), method, pvalue, log2FC, status]

        if compute_effect:
            # specific named column + generic 'Effect size' + 'Effect size type'
            # are contiguous in base_cols (inserted right after 'P').
            pos = base_cols.index(effect_specific_col)
            row[pos:pos] = [effect_val, effect_val, effect_type]

        result_df_index.append(each_index)
        for col, val in zip(base_cols, row):
            result_df_dict[col].append(val)

    result_df_columns = list(base_cols)

    # --- p-value correction columns (conditional on filter_by) ---
    p_insert_pos = result_df_columns.index('P') + 1
    bonf_col = ''
    if filter_by == 'FDR':
        _, p_value_corrected = fdrcorrection(result_df_dict['P'])
        result_df_columns.insert(p_insert_pos, 'FDR (adjusted P)')
        result_df_dict['FDR (adjusted P)'] = list(p_value_corrected)
        fdr_col = 'FDR (alpha=%s)' % str(p_value_cutoff)
        result_df_columns.insert(p_insert_pos + 1, fdr_col)
        fdr_passed = np.array(p_value_corrected) < p_value_cutoff
        result_df_dict[fdr_col] = [
            str(i).replace('True', 'Passed').replace('False', 'Failed')
            for i in fdr_passed
        ]
    elif filter_by == 'raw_p':
        rawp_col = 'P (alpha=%s)' % str(p_value_cutoff)
        result_df_columns.insert(p_insert_pos, rawp_col)
        rawp_passed = np.array(result_df_dict['P']) < p_value_cutoff
        result_df_dict[rawp_col] = [
            str(i).replace('True', 'Passed').replace('False', 'Failed')
            for i in rawp_passed
        ]
    elif filter_by == 'Bonferroni':
        alpha = p_value_cutoff / len(result_df_index) if result_df_index else p_value_cutoff
        bonf_col = 'Bonferroni (alpha=%s)' % str(alpha)
        result_df_columns.insert(p_insert_pos, bonf_col)
        bonf_passed = np.array(result_df_dict['P']) < alpha
        result_df_dict[bonf_col] = [
            str(i).replace('True', 'Passed').replace('False', 'Failed')
            for i in bonf_passed
        ]

    # --- primary p-value filter ---
    if filter_by == 'FDR':
        filter_array = np.array(result_df_dict['FDR (adjusted P)']) < p_value_cutoff
    elif filter_by == 'raw_p':
        filter_array = np.array(result_df_dict['P']) < p_value_cutoff
    elif filter_by == 'Bonferroni':
        filter_array = np.array(result_df_dict[bonf_col]) == 'Passed'
    else:
        filter_array = np.ones(len(result_df_index), dtype=bool)

    # --- additional effect-size filter (decoupled from the test choice) ---
    # Non-parametric: |rank-biserial| > rbc_cutoff (paper Section 2.10: 0.5).
    # Parametric:     |Cohen's d|     > cohens_d_cutoff (e.g. 0.8 = "large").
    # The measures are on different scales, so the two cutoffs are never shared.
    apply_effect_filter = False
    if compute_effect:
        if effect_type == 'rank_biserial' and rbc_cutoff is not None:
            effect_cutoff = rbc_cutoff
        elif effect_type == 'cohens_d' and cohens_d_cutoff is not None:
            effect_cutoff = cohens_d_cutoff
        else:
            effect_cutoff = None
        if effect_cutoff is not None:
            eff_array = np.abs(np.array(result_df_dict['Effect size'], dtype=float))
            filter_array = filter_array & (eff_array > effect_cutoff)
            apply_effect_filter = True

    if filter_by != 'None' or apply_effect_filter:
        for col in result_df_columns:
            result_df_dict[col] = list(
                np.array(result_df_dict[col])[filter_array]
            )
        result_df_index = list(np.array(result_df_index)[filter_array])

    result_df = pd.DataFrame(
        data=result_df_dict,
        columns=result_df_columns,
        index=result_df_index
    )
    return result_df
