import numpy as np
from . import config

def calculate_metrics(parsed_data, sel_list, threshold, score_is_pval=False):
    """
    Computes significance, mean predicted selection, and confusion matrices.
    Returns a dictionary of metrics keyed by true selection coefficient.
    This is perfectly structured for plotting with matplotlib/seaborn!
    """
    bounds = config.get_s_grid_bounds()
    metrics_out = {sel: [] for sel in sel_list}

    for sel in sel_list:
        assert sel in parsed_data
        #assert len(parsed_data[sel]) > 0
        #if sel not in parsed_data or len(parsed_data[sel]) == 0:
        #    continue

        runs = parsed_data[sel]
        n_total = len(runs)

        # 1. Determine Significance based on Threshold
        if score_is_pval:
            sig_runs = [r for r in runs if r['score'] < threshold]
        else:
            sig_runs = [r for r in runs if r['score'] >= threshold]

        n_sig = len(sig_runs)
        perc_sig = (n_sig / n_total) * 100 if n_total > 0 else 0

        # 2. Build Confusion Matrix & extract plotting arrays
        confusion = {s: 0 for s in set(config.SWEEPLINK_GRID)}
        pred_s_values_sig = []
        all_final_preds =[]  # Includes the ones forced to 0.0

        for r in runs:
            # Force to neutral if it didn't pass the threshold
            is_sig = (r['score'] < threshold) if score_is_pval else (r['score'] >= threshold)
            
            if not is_sig:
                final_pred_s = 0.0
            else:
                final_pred_s = r['pred_s']
                pred_s_values_sig.append(final_pred_s)
                
            all_final_preds.append(final_pred_s)

            # Bin the predicted selection into the confusion matrix
            for i, b in enumerate(bounds):
                if b[0] <= final_pred_s < b[1]:
                    confusion[config.SWEEPLINK_GRID[i]] += 1
                    break

        # 3. Store everything in a structured dictionary
        metrics_out[sel] = {
            "n_total": n_total,
            "n_sig": n_sig,
            "perc_sig": perc_sig,
            "confusion": confusion,
            "mean_pred_s": np.mean(pred_s_values_sig) if pred_s_values_sig else None,
            "std_pred_s": np.std(pred_s_values_sig) if pred_s_values_sig else None,
            
            # --- ARRAYS FOR PLOTTING ---
            "raw_scores": [r['score'] for r in runs],         # For ROC curves / p-value histograms
            "raw_sig_preds": pred_s_values_sig,               # For Boxplots of only significant runs
            "all_final_preds": all_final_preds,               # For Boxplots of ALL runs (including 0.0s)
        }
    return metrics_out

def get_error(y_pred, true_val):
    return (np.array(y_pred) - true_val) / true_val

def get_TPR_per_sel_and_FPR(parsed_data, sel_list, threshold, score_is_pval=False, return_metrics=True):
    assert 0.0 in sel_list
    metrics = calculate_metrics(
        parsed_data=parsed_data,
        sel_list=sel_list,
        threshold=threshold,
        score_is_pval=score_is_pval
    )
    # 2. get TPR for each s
    TPR = {}
    for sel in set(sel_list):
        if sel == 0.0:
            FPR = metrics[sel]["n_sig"] / metrics[sel]["n_total"] if metrics[sel]["n_total"] > 0 else 0
        else:
            TPR[sel] = metrics[sel]["n_sig"] / metrics[sel]["n_total"] if metrics[sel]["n_total"] > 0 else 0
    if return_metrics:
        return TPR, FPR, metrics
    return TPN, FPR

def get_MCC(parsed_data, sel_list, threshold, score_is_pval=False, return_metrics=False):
    assert 0.0 in sel_list

    metrics = calculate_metrics(
        parsed_data=parsed_data,
        sel_list=sel_list,
        threshold=threshold,
        score_is_pval=score_is_pval
    )

    TP, FP, TN, FN = 0, 0, 0, 0
    for sel in sel_list:
        if sel == 0.0:
            TN += metrics[sel]["n_total"] - metrics[sel]["n_sig"]
            FP += metrics[sel]["n_sig"]
        else:
            TP += metrics[sel]["n_sig"]
            FN += metrics[sel]["n_total"] - metrics[sel]["n_sig"]

    # Formula for MCC: (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
    numerator = (TP * TN) - (FP * FN)
    denominator = np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

    # Handle edge case: if denominator is 0, MCC is defined as 0
    if denominator == 0:
        return 0.0
    if return_metrics:
        return numerator / denominator, metrics
    return numerator / denominator
