import os
import numpy as np
from . import config
import pandas as pd
import tqdm

def get_true_sel_sweeplink(chrom_name):
    return float(chrom_name.split("_")[-3])

def parse_sweeplink_s_file(filepath, pos_of_sel):
    """Returns a dictionary with per-chromosome predictions."""
    parsed = {}
    if not os.path.exists(filepath): 
        return parsed
        
    with open(filepath) as f:
        try: next(f)
        except: return parsed
        
        prevpos, curchrom = -1, ""
        for line in f:
            if not line.startswith("s_"): continue
            parts = line.split()
            chrom_full = parts[0]
            if curchrom != chrom_full: prevpos = -1
            
            is_neut = float(chrom_full.split("_")[-3]) == 0
            pos = int(chrom_full.split("_")[-2])
            
            probs_len = None
            if is_neut or pos == pos_of_sel:
                probs = np.array([float(x) for x in parts[1:]])
                if probs_len is None:
                    probs_len = len(probs)
                if len(probs) != probs_len:
                    continue  # Run is not finished
                
                middle_ind = len(probs) // 2
                p_neut = probs[middle_ind - 1] + probs[middle_ind]
                p_pos = sum([p for ind, p in enumerate(probs) if ind < middle_ind - 1])
                p_neg = sum([p for ind, p in enumerate(probs) if ind > middle_ind])
                
                p_sum = p_neut + p_pos + p_neg
                if p_sum > 0:
                    p_pos, p_neg = p_pos/p_sum, p_neg/p_sum
                
                score = max(p_pos, p_neg)
                pred_s = np.sum([p*sval for p,sval in zip(probs, config.SWEEPLINK_GRID)])
                
                if chrom_full not in parsed:
                    parsed[chrom_full] = []
                parsed[chrom_full].append({
                    'score': score,
                    'pred_s': pred_s
                })
            prevpos = pos
    return parsed


def parse_sweeplink_trace_file(filepath):
    if not os.path.exists(filepath):
        return []
    df =  pd.read_csv(filepath, header=0, sep="\t", on_bad_lines='skip')
    return df.to_dict(orient="list")


def get_data_per_s_accross_repeats(char_name, char_val, sel_list,  n_total=100):
    char_config = config.get_char_config(char_name)
    print(f"Start loading results for {char_config['x_label']} = {char_val} (Reps 0 to {n_total})")
    n_full = 0
    n_empty = 0

    # Aggregate runs across all simulation folders
    aggregated_data = {sel: [] for sel in sel_list}

    for ind in tqdm.tqdm(range(n_total)):
        out_dir = os.path.join(config.get_inference_dir_for_val(char_name, char_val, tool_name="sweeplink", is_comparison=False), str(ind))
        post_file = os.path.join(out_dir, "sweepLink_s_statePosteriors.txt")
        pos_of_sel = config.get_chrom_middle_pos(char_name, char_val)
        chrom2results = parse_sweeplink_s_file(post_file, pos_of_sel=pos_of_sel)

        if len(chrom2results) == 0:
            n_empty += 1
        else:
            n_full += config.is_finished_sweeplink(out_dir)

        # Group by true selection coefficient
        for chrom_name in chrom2results:
            # the chrom name starts with additional s_ in output files
            true_s = config.get_sel_from_chrom_name(chrom_name[2:])
            assert true_s in aggregated_data
            aggregated_data[true_s].extend(chrom2results[chrom_name])
    print(f"Finished loading. Total number of loaded runs: {n_total - n_empty}, Number of finished runs: {n_full}")
    return aggregated_data


def get_trace_data_accross_repeats(char_name, char_val, n_total=100):
    all_posterior = None
    char_config = config.get_char_config(char_name)
    print(f"Start loading Ne samples for {char_config['x_label']} = {char_val} (Reps 0 to {n_total})")
    n_full = 0
    n_empty = 0
    for ind in tqdm.tqdm(range(n_total)):
        out_dir = os.path.join(config.get_inference_dir_for_val(char_name, char_val, tool_name="sweeplink", is_comparison=False), str(ind))
        post_file = os.path.join(out_dir, "sweepLink_trace.txt")
        ind_posterior = parse_sweeplink_trace_file(post_file)
        if len(ind_posterior) == 0:
            n_empty += 1
        else:
            n_full += config.is_finished_sweeplink(out_dir)
        if len(ind_posterior) == 0:
            continue
        if all_posterior is None:
            all_posterior = ind_posterior
        else:
            for key in all_posterior:
                all_posterior[key].extend(ind_posterior[key])
    # remove None in array
    for key in all_posterior:
        if len(all_posterior[key]) == 0:
            continue
        arr = np.array(all_posterior[key], dtype=float)
        all_posterior[key] = arr[~np.isnan(arr)]

    print(f"Finished loading. Total number of loaded runs: {n_total - n_empty}, Number of finished runs: {n_full}") 
    return all_posterior
