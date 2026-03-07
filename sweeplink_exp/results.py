import os
import numpy as np
from . import config, comparison
import pandas as pd
import tqdm
from scipy.stats import gamma
from scipy.stats import chi2

def get_true_sel_sweeplink(chrom_name):
    return float(chrom_name.split("_")[-3])

def parse_sweeplink_s_trace_file(filepath, pos_of_sel):
    data_dict = parse_sweeplink_trace_file(filepath)

    chrom2data = {}
    for key in data_dict:
        is_neut = config.get_sel_from_chrom_name(key[2:]) == 0
        pos = int(key.split("_")[-2])
        if is_neut or pos == pos_of_sel:
            raw_array = np.array(data_dict[key], dtype=float)
            arr = raw_array[~np.isnan(raw_array)]
            probs = np.bincount(arr.astype(int), minlength=len(config.SWEEPLINK_GRID))
            probs = probs.astype(float) / np.sum(probs)
            if key not in chrom2data:
                chrom2data[key] = []
            chrom2data[key].append(get_score_pred_s_from_probs(probs))
    return chrom2data

def get_score_pred_s_from_probs(probs):
    middle_ind = len(probs) // 2
    p_neut = probs[middle_ind - 1] + probs[middle_ind]
    p_pos = sum([p for ind, p in enumerate(probs) if ind < middle_ind - 1])
    p_neg = sum([p for ind, p in enumerate(probs) if ind > middle_ind])

    p_sum = p_neut + p_pos + p_neg
    if p_sum > 0:
        p_pos, p_neg = p_pos/p_sum, p_neg/p_sum

        score = max(p_pos, p_neg)
        pred_s = np.sum([p*sval for p,sval in zip(probs, config.SWEEPLINK_GRID)])

        return {
            'score': score,
            'pred_s': pred_s
        }
    return {"score": None, "pred_s": None}


def parse_sweeplink_s_posterior_file(filepath, pos_of_sel):
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
            # it starts with s_
            chrom_full = parts[0][2:]
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
                
                if chrom_full not in parsed:
                    parsed[chrom_full] = []
                parsed[chrom_full].append(get_score_pred_s_from_probs(probs))
            prevpos = pos
    return parsed


def parse_sweeplink_output(directory, pos_of_sel, trace_file=False):
    if trace_file:
        post_file = os.path.join(directory, "sweepLink_s_trace.txt")
        return parse_sweeplink_s_trace_file(post_file, pos_of_sel=pos_of_sel)
    post_file = os.path.join(directory, "sweepLink_s_statePosteriors.txt")
    return parse_sweeplink_s_posterior_file(post_file, pos_of_sel=pos_of_sel)


def parse_sweeplink_trace_file(filepath, header=0):
    if not os.path.exists(filepath):
        return []
    df =  pd.read_csv(filepath, header=header, sep="\t", on_bad_lines='skip')
    df = df.apply(pd.to_numeric, errors='coerce').dropna()
    return df.to_dict(orient="list")


def convert_approxwf_output(directory):
    filename = os.path.join(directory, comparison.get_tool_output_file("approxwf"))
    output_file = os.path.join(directory, "approxwf_meanVar.txt")
    if os.path.exists(output_file):
        return output_file
    print("Converting", filename)
    columns = ["Parameter", "P(x=0)", "P(x>0)", "P(x<0)", "Mean_s", "Median_s","Std_s"]
    data = {k: [] for k in columns}
    loaded_dict = parse_sweeplink_trace_file(filename)
    for key in loaded_dict:
        samples = np.array(loaded_dict[key])
        p_pos = np.sum(samples > 0) / len(samples)
        p_neg = np.sum(samples < 0) / len(samples)
        p_zero = np.sum(samples == 0) / len(samples)
        mean = np.mean(samples)
        median = np.median(samples)
        std = np.std(samples)
        vals = [key, p_zero, p_pos, p_neg, mean, median, std]
        for k, v in zip(columns, vals):
            data[k].append(v)
    df = pd.DataFrame(data)
    df.to_csv(output_file, sep='\t', index=False)
    return output_file

def parse_approxwf_output(directory, pos_of_sel):
    filename = convert_approxwf_output(directory)
    parsed = {}
    with open(filename) as f:
        for line in f:
            parameter = line.split()[0]
            if not parameter.startswith("s_"):
                continue
            #s_2_sel_0.0_1131
            chrom = "_".join(line.split()[0].split("_")[1:-1])
            is_neut = config.get_sel_from_chrom_name(chrom) == 0
            pos = int(line.split()[0].split("_")[-1])

            if is_neut or pos == pos_of_sel:
                if chrom not in parsed:
                    parsed[chrom] = []
                p_pos = float(line.split()[2])
                p_neg = float(line.split()[3])
                score = max(p_pos, p_neg)
                pred_s = float(line.split()[4])
                parsed[chrom].append({
                    'score': score,
                    'pred_s': pred_s
                })
    return parsed


def parse_diplolocus_output(directory, pos_of_sel):
    # Header
    # ID      ongrid_s2hat    ongrid_maxLogLikelihood s2hat   maxLogLikelihood        MLR     chi2_p
    filename = os.path.join(directory, comparison.get_tool_output_file("diplolocus"))
    parsed = {}
    if not os.path.exists(filename):
        return parsed
    with open(filename) as f:
        for line in f:
            if line.startswith("#") or line.startswith("ID"):
                continue
            chrom = "_".join(line.split()[0].split("_")[:-1])
            pos = int(line.split()[0].split("_")[-1])
            is_neut = config.get_sel_from_chrom_name(chrom) == 0
            if is_neut or pos == pos_of_sel:
                if chrom not in parsed:
                    parsed[chrom] = []
                pred_s = float(line.split()[3])
                score = chi2.sf(float(line.split()[5]), df=1)
                parsed[chrom].append({
                    'score': score,
                    'pred_s': pred_s
                })
    return parsed

def parse_bmws_output(directory, pos_of_sel):
    filename = os.path.join(directory, comparison.get_tool_output_file("bmws"))
    parsed = {}
    if not os.path.exists(filename):
        print(filename)
        return parsed
    # Header
    # chrom pos rs_id ref alt AF s_cap l1-norm(s), l2-norm(s) not_used not_used
    # 1_sel_0.01      50001   .       A       T       0.402   0.00456 0.00648 0.004605
    WIN_SIZE = 20
    # 1. load data for each chrom
    neut_data = {}
    with open(filename) as f:
        for line in f:
            chrom = line.split()[0]
            pos = int(line.split()[1])
            is_neut = config.get_sel_from_chrom_name(chrom) == 0
            if is_neut:
                if chrom not in neut_data:
                    neut_data[chrom] = []
                pred_s = float(line.split()[6])
                neut_data[chrom].append([pos, pred_s])
    # calculate gamma distribution from neutral
    WIN_SIZE = 20
    HALF_WIN_SIZE = WIN_SIZE // 2
    window_s_vals = []
    chrom_pos_to_val = {}
    for chrom in neut_data:
        for i in range(len(neut_data[chrom]) // HALF_WIN_SIZE - 1):
            arr = np.array(neut_data[chrom][i*HALF_WIN_SIZE:i*HALF_WIN_SIZE+WIN_SIZE])
            val = np.mean(arr[:, 1])
            window_s_vals.append(val)
            for pos, _val in arr:
                key = (chrom, pos)
                if key not in chrom_pos_to_val:
                    chrom_pos_to_val[key] = []
                chrom_pos_to_val[key].append(val)
    mu = np.mean(window_s_vals)
    var = np.var(window_s_vals)
    shape = (mu**2) / var
    scale = var / mu

    # One more time with p-value
    pos = None
    with open(filename) as f:
        for line in f:
            chrom = line.split()[0]
            if pos == int(line.split()[1]):
                continue
            pos = int(line.split()[1])
            is_neut = config.get_sel_from_chrom_name(chrom) == 0
            if is_neut or pos == pos_of_sel:
                if chrom not in parsed:
                    parsed[chrom] = []
                pred_s = float(line.split()[6])
                if pos == pos_of_sel:
                    score = gamma.sf(pred_s, a=shape, scale=scale)
                else:
                    score = None
                    key = (chrom, pos)
                    if key not in chrom_pos_to_val:
                        continue
                    for val in chrom_pos_to_val[key]:
                        pval = gamma.sf(val, a=shape, scale=scale)
                        if score is None or pval < score:
                            score = pval
                parsed[chrom].append({
                    'score': score,
                    'pred_s': pred_s
                })
    return parsed


def parse_output(char_name, char_val, ind, trace_file):
    out_dir = os.path.join(config.get_inference_dir_for_val(char_name, char_val), str(ind))
    pos_of_sel = config.get_chrom_middle_pos(char_name, char_val)

    tool_name = config.get_tool_name(char_name, char_val)
    if tool_name == "sweeplink":
        return parse_sweeplink_output(
            directory=out_dir,
            pos_of_sel=pos_of_sel,
            trace_file=trace_file
        )
    if tool_name == "approxwf":
        return parse_approxwf_output(
            directory=out_dir,
            pos_of_sel=pos_of_sel,
        )
    if tool_name == "diplolocus":
        return parse_diplolocus_output(
            directory=out_dir,
            pos_of_sel=pos_of_sel,
        )
    assert tool_name == "bmws"
    return parse_bmws_output(
        directory=out_dir,
        pos_of_sel=pos_of_sel,
    )


def get_data_per_s_for_ind(char_name, char_val, ind, s_trace_file=False):
    chrom2results = parse_output(char_name, char_val, ind, s_trace_file)
    data = {}
    for chrom_name in chrom2results:
        # the chrom name starts with additional s_ in output files
        true_s = config.get_sel_from_chrom_name(chrom_name)
        if true_s not in data:
            data[true_s] = []
        data[true_s].extend(chrom2results[chrom_name])
    return data

def get_data_per_s_accross_repeats(char_name, char_val, n_total=100):
    char_config = config.get_char_config(char_name)
    print(f"Start loading results for {char_config['x_label']} = {char_val} (Reps 0 to {n_total})")
    n_full = 0
    n_empty = 0

    # Aggregate runs across all simulation folders
    aggregated_data = {sel: [] for sel in (config.get_sel_list(char_name) + [0.0])}

    for ind in tqdm.tqdm(range(n_total)):
        sel2results = get_data_per_s_for_ind(char_name, char_val, ind)
        if len(sel2results) == 0:
            n_empty += 1
        else:
            out_dir = os.path.join(config.get_inference_dir_for_val(char_name, char_val), str(ind))
            n_full += config.is_finished_sweeplink(out_dir)

        for true_s in sel2results:
            assert true_s in aggregated_data
            aggregated_data[true_s].extend(sel2results[true_s])
    print(f"Finished loading. Total number of loaded runs: {n_total - n_empty}, Number of finished runs: {n_full}")
    return aggregated_data

def get_trace_data_for_ind(char_name, char_val, ind):
    out_dir = os.path.join(config.get_inference_dir_for_val(char_name, char_val), str(ind))
    post_file = os.path.join(out_dir, "sweepLink_trace.txt")
    ind_posterior = parse_sweeplink_trace_file(post_file)
    return ind_posterior

def get_trace_data_accross_repeats(char_name, char_val, n_total=100):
    all_posterior = None
    char_config = config.get_char_config(char_name)
    print(f"Start loading Ne samples for {char_config['x_label']} = {char_val} (Reps 0 to {n_total})")
    n_full = 0
    n_empty = 0
    for ind in tqdm.tqdm(range(n_total)):
        out_dir = os.path.join(config.get_inference_dir_for_val(char_name, char_val), str(ind))
        ind_posterior = get_trace_data_for_ind(char_name, char_val, ind)
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


# --- ApproxWF ---

def parse_approxwf_output_file(filename, pos_of_sel):
    data_dict = parse_sweeplink_trace_file(filename)

