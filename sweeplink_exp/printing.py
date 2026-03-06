from . import config, results, metrics
import numpy as np

def print_metrics(data, threshold, score_is_pval=False):
    """
    Takes the computed dictionary from calculate_metrics and formats it for the terminal.
    """
    condition_str = f"p-value < {threshold}" if score_is_pval else f"prob >= {threshold}"
    print(f"  Runs significant ({condition_str}): {data['n_sig']}/{data['n_total']} ({int(data['perc_sig'])}%)")

    if data['mean_pred_s'] is not None:
        print(f"  Mean pred selection (of sig runs): {data['mean_pred_s']:.5f} +- {data['std_pred_s']:.5f}")

    cm_str = " ".join([("| " if s == 0.0 else "") + str(data['confusion'][s]) + (" |" if s == 0.0 else "") for s in sorted(list(set(config.SWEEPLINK_GRID)))])
    print(f"  Confusion matrix: {cm_str}")


def process_and_print_results_per_sel_coef(char_name, n_total):
    """Processes characteristic experiments (SweepLink only)."""
    char_config = config.get_char_config(char_name)
    char_values = char_config['values']
    sel_list = config.SIM_MAIN_LIST_WITH_0

    aggregated_data = {}
    for char_val in char_values:
        aggregated_data[char_val] = results.get_data_per_s_accross_repeats(
            char_name=char_name,
            char_val=char_val,
            sel_list=sel_list,
            n_total=n_total,
        )
    for sel in sel_list:
        print(f"\n\n*** Results for s = {sel} ***")

        for char_val in char_values:
            metrics_data = metrics.calculate_metrics(
                parsed_data=aggregated_data[char_val],
                sel_list=[sel],
                threshold=config.THRESHOLDS['sweeplink'],
                score_is_pval=False,
            )
            print(f"\n{char_config['x_label']}: {char_val}")
            print_metrics(
                data=metrics_data[sel],
                threshold=config.THRESHOLDS['sweeplink'],
                score_is_pval=False
            )

def process_and_print_results_per_char_value(char_name, n_total):
    """Processes characteristic experiments (SweepLink only)."""
    char_config = config.get_char_config(char_name)
    char_values = char_config['values']
    sel_list = config.SIM_MAIN_LIST_WITH_0
    
    for char_val in char_values:
        print(f"\n\n*** Results for {char_config['x_label']} = {char_val} ***")

        aggregated_data = results.get_data_per_s_accross_repeats(
            char_name=char_name,
            char_val=char_val,
            sel_list=sel_list,
            n_total=n_total,
        )

        metrics_data = metrics.calculate_metrics(
            parsed_data=aggregated_data,
            sel_list=sel_list,
            threshold=config.THRESHOLDS['sweeplink'],
            score_is_pval=False,
        )

        for sel in sel_list:
            print(f"\nSelection: {sel}")
            print_metrics(
                data=metrics_data[sel],
                threshold=config.THRESHOLDS['sweeplink'],
                score_is_pval=False
            )

def process_and_print_intermediate_results(char_name, char_val, ind, is_comparison):
    char_config = config.get_char_config(char_name)
    sel_list = config.SIM_MAIN_LIST_WITH_0

    # selection
    data = results.get_data_per_s_for_ind(char_name, char_val, sel_list, ind, is_comparison, s_trace_file=True)
    metrics_data = metrics.calculate_metrics(
        parsed_data=data,
        sel_list=sel_list,
        threshold=config.THRESHOLDS['sweeplink'],
        score_is_pval=False,
    )
    for sel in sel_list:
        print(f"\nSelection: {sel}")
        print_metrics(
            data=metrics_data[sel],
            threshold=config.THRESHOLDS['sweeplink'],
            score_is_pval=False
        )

    # Pop size
    posterior = results.get_trace_data_for_ind(char_name, char_val, ind, is_comparison)['N_pop0']
    posterior = np.array(posterior)
    posterior = posterior[~np.isnan(posterior)]
    print(f"\nPopulation size: {np.mean(posterior):.1f} +- {np.std(posterior):.1f}")
