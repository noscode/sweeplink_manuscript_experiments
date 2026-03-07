import os
from matplotlib import pyplot as plt

# --- PATHS ---
SWEEPLINK_PATH = "/home/enoskova/sweeplink/sweepLink"
APPROXWF_PATH = "/home/enoskova/approxwf/ApproxWF"
BASE_EXP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "experiments"))

# --- CLUSTER SETTINGS ---
SLURM_EXCLUDE = "binfservas[11-29],binfservas33,binfservas34,binfservas35,binfservas36,binfservas37,binfservas38,binfservas39,binfservas40"

# --- EXTERNAL STORED SIMULATION PATHS ---
# Pointing to your existing data to save disk space
EXTERNAL_SIM_PATH_1 = "/data/users/enoskova/approxwf3_experiments_final"
EXTERNAL_SIM_PATH_2 = "/data/projects/p243_ancientdna_unifr/Katya"

def get_simulation_dir(char_name, char_val, true_s):
    """
    Returns the directory containing 'logs/' and 'vcfs/' for a given setup.
    This reads directly from your previously simulated cluster paths!
    """
    # 0.03 and 0.04 were stored in the second path
    if true_s in[0.03, 0.04]:
        base_path = EXTERNAL_SIM_PATH_2
    else:
        base_path = EXTERNAL_SIM_PATH_1

    # Example output: /data/.../sample_size/sel_0.01/data_25
    if char_name == "n_grid":
        data_dir = "data"
    else:
        data_dir = f"data_{char_val}"
    return os.path.join(base_path, char_name, f"sel_{true_s}", data_dir)

# --- EXPERIMENT SETTINGS ---
SIM_MAIN_LIST = [0.01, 0.02, 0.05]
SIM_MAIN_LIST_WITH_0 = [0.0, 0.01, 0.02, 0.05]
SIM_COMPARE_LIST = [0.01, 0.02, 0.03, 0.04, 0.05]
COMPARE_SEL_LIST = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05]
COMPARE_N_SIM = 100
THRESHOLDS = {
    'sweeplink': 1.00,    # Probability threshold
    'approxwf': 0.847778096527895,     # Probability threshold
    'diplolocus': 7.56463327554629e-05,  # P-value threshold
    'bmws': 0.0031257158496882415,
}
IS_PVAL = {
    'sweeplink': False,
    'approxwf': False,
    'diplolocus': True,
    'bmws': True,
}

def get_sim_sel_list(char_name, val):
    """Only return the 5 coefficients if we are strictly running sample_size = 25."""
    if char_name == 'sample_size' and val == 25:
        return SIM_COMPARE_LIST
    return SIM_MAIN_LIST

def get_sel_str(sel, scenario):
    if scenario == "neut":
        return "0.0"
    return f"{sel}"

def get_eval_label(sel, scenario):
    return "sel_" + get_sel_str(sel, scenario)

def get_chrom_name(i_chrom, sel, scenario):
    return f"{i_chrom}_{get_eval_label(sel, scenario)}"

def get_sel_from_chrom_name(chrom_name):
    return float(chrom_name.split("_")[2])

def get_sel_scenario_combinations(char_name):
    sel_list = get_sel_list(char_name)
    sources = []
    for s in sel_list:
        sources.append({"true_s": s, "scenario": "sdn"})
        sources.append({"true_s": s, "scenario": "neut"})
    return sources


# --- CONSTANTS ---
S_POS_GRID = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06]
SWEEPLINK_GRID = [-s / (1+s) for s in reversed(S_POS_GRID)]
SWEEPLINK_GRID.append(0)  # Two zeros
SWEEPLINK_GRID.append(0)
SWEEPLINK_GRID.extend(S_POS_GRID)

def get_s_grid_bounds():
    bounds =[]
    start = -1
    for i in range(len(SWEEPLINK_GRID) - 1):
        end = (SWEEPLINK_GRID[i] + SWEEPLINK_GRID[i+1]) / 2
        bounds.append((start, end))
        start = end
    bounds.append((start, 1))
    return bounds

# --- PLOTTING ---
# For master plots (sweeplink only)
VIRIDIS_CMAP = plt.get_cmap('viridis')
SEL_COLORS = {0.0: VIRIDIS_CMAP(1.0), 0.01: VIRIDIS_CMAP(0.8), 0.02: VIRIDIS_CMAP(0.6), 0.03: VIRIDIS_CMAP(0.4), 0.04: VIRIDIS_CMAP(0.2), 0.05: VIRIDIS_CMAP(0.0)}
SEL_LABELS = {0.01: "Weak Selection", 0.02: "Moderate Selection", 0.05: "Strong Selection"}
SEL_ZORDER = {0.0: 2, 0.01: 3, 0.02: 4, 0.03: 5, 0.04: 6, 0.05: 7}

# tool comparison
TAB_CMAP = plt.get_cmap("tab10")
TOOL_LABELS = {"sweeplink": 'SweepLink', "bmws": 'bmws', "diplolocus": "diplo-locus", "approxwf": "ApproxWF"}
TOOL_COLORS =  {"sweeplink": TAB_CMAP(0), "bmws": TAB_CMAP(1), "diplolocus": TAB_CMAP(2), "approxwf": TAB_CMAP(3)}
TOOL_ZORDER = {"sweeplink": 4, "bmws": 3, "diplolocus": 2, "approxwf": 1}


# --- CHARACTERISTICS ---
characteristics_list = [
    {
        'name': 'sample_size',
        'values': [5, 10, 25, 50, 100, 250, 500],
        'x_label': 'Sample size',
        'xlog': False,
        'default_value': 25,
        'title': 'across Sample Sizes',
        'can_reuse_simulation': False,
    },
    {
        'name': 'n_loci',
        'values': [0.05, 0.1, 0.2, 0.5, 1.0],
        'val2str': {0.05: "50k", 0.1: "100k", 0.2: "200k", 0.5: "500k", 1.0: "1m"},
        'x_label': 'Sequence length of one chromosome (Mb)',
        'xlog': False,
        'default_value': 0.1,
        'title': 'across sequence length',
        'can_reuse_simulation': False,
    },
    {
        'name': 'tpoints_num',
        'values': [6, 11, 21, 41],
        'val2str': {6: "5", 11: "10", 21: "20", 41: "40"},
        'x_label': 'Number of time points',
        'xlog': False,
        'default_value': 11,
        'title': 'across number of time points',
        'can_reuse_simulation': False,
    },
    {
        'name': 'binning',
        'values': [2, 4, 8, 16, 32, 40, 80],
        'x_label': 'Binning size (generations)',
        'xlog': False,
        'default_value': 16,
        'title': 'across binning sizes',
        'can_reuse_simulation': True,
    },
    {
        'name': 'start_freq',
        'values': [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
        'x_label': 'Initial frequency',
        'xlog': False,
        'default_value': 0.3,
        'title': 'across initial frequencies',
        'can_reuse_simulation': False,
    },
    {
        'name': 'rec_rate',
        'values': [1e-9, 1e-8, 1e-7, 1e-6],
        'val2str': {1e-10: "10", 1e-9: "9", 1e-8: "8", 1e-7: "7", 1e-6: "6", 1e-5: "5"},
        'x_label': 'Recombination rate (cM / Mb)',
        'xlog': True,
        'default_value': 1e-8,
        'title': 'across recombination rates',
        'can_reuse_simulation': False,
    },
    {
        'name': 'n_grid',
        'values': [10, 20, 50, 75, 100],
        'x_label': 'Grid size',
        'xlog': False,
        'default_value': 50,
        'title': 'across grid sizes',
        'can_reuse_simulation': True,
    },
    {
        'name': 'comparison',
        'values': ["sweeplink", "approxwf", "diplolocus", "bmws"],
        'x_label': 'Tool Name',
        'xlog': False,
        'default_value': "sweeplink",
        'title': 'across tools',
        'can_reuse_simulation': True,
    },
]
for i in range(len(characteristics_list)):
    if "val2str" not in characteristics_list[i]:
        characteristics_list[i]["val2str"] = {val: str(val) for val in characteristics_list[i]["values"]}

# We will use this final structure
characteristics = {x["name"]: x for x in characteristics_list}

def get_char_config(char_name):
    if char_name not in characteristics:
        raise ValueError(f"Characteristic '{char_name}' not found.")
    return characteristics[char_name]

def get_defaults():
    return {c['name']: c['default_value'] for c in characteristics}

def get_char_name_for_comp():
    return "sample_size"

def get_char_val_for_comp():
    return characteristics["sample_size"]["default_value"]

def is_finished_sweeplink(directory):
    check_file = os.path.join(directory, "sweepLink_meanVar.txt")
    if not os.path.exists(check_file):
        return False
    with open(check_file, "r") as f:
        num_lines = len(f.readlines())
    return num_lines > 1

# --- Directories ---
def get_experiment_dir(char_name):
    return os.path.join(BASE_EXP_DIR, char_name)

def get_comparison_dir():
    return os.path.join(BASE_EXP_DIR, "tools_comparison")

def get_inference_dir(char_name):
    return os.path.join(get_experiment_dir(char_name), "inference")

def get_inference_dir_for_val(char_name, char_val):
    return os.path.join(get_inference_dir(char_name), f"val_{characteristics[char_name]['val2str'][char_val]}")

def get_simulation_dir(char_name):
    return os.path.join(get_experiment_dir(char_name), "simulation")

def get_sim_char_name_val(char_name, char_val):
    if char_name == "comparison":
        sim_char_name = "sample_size"
        return sim_char_name, get_char_config(sim_char_name)["default_value"]
    return char_name, char_val

def get_sel_list(char_name):
    if char_name == "comparison":
        return SIM_COMPARE_LIST
    return SIM_MAIN_LIST

def get_tool_name(char_name, char_val):
    if char_name == "comparison":
        return char_val
    return "sweeplink"

def get_counts_filename(char_name, char_val):
    from . import comparison
    tool_name = get_tool_name(char_name, char_val)
    return comparison.get_tool_input_file(tool_name=tool_name)

def get_generation_func_meta(char_name, char_val):
    from . import comparison
    tool_name = get_tool_name(char_name, char_val)
    return comparison.get_generation_func_meta_for_tool(tool_name=tool_name)

def get_generation_func_counts_header(char_name, char_val):
    from . import comparison
    tool_name = get_tool_name(char_name, char_val)
    return comparison.get_generation_func_counts_header(tool_name=tool_name)

def get_generation_func_counts_line(char_name, char_val):
    from . import comparison
    tool_name = get_tool_name(char_name, char_val)
    return comparison.get_generation_func_counts_line(tool_name=tool_name)

def get_generation_func_cmd(char_name, char_val):
    from . import comparison
    tool_name = get_tool_name(char_name, char_val)
    return comparison.get_generation_func_cmd(tool_name=tool_name)

def get_data_dir_name(char_name, char_val):
    sim_char_name, sim_char_val = get_sim_char_name_val(char_name, char_val)
    if characteristics[sim_char_name]["can_reuse_simulation"]:
        return "data"
    return f"data_{get_char_config(sim_char_name)['val2str'][sim_char_val]}"

def get_data_dir(char_name, char_val, sel_coef):
    #sim_dir = get_simulation_dir(char_name)
    if sel_coef in [0.03, 0.04]:
        sim_dir = f"/data/projects/p243_ancientdna_unifr/Katya/{char_name}/"
    else:
        sim_dir = f"/data/users/enoskova/approxwf3_experiments_final/{char_name}/"
    data_dir_name = get_data_dir_name(char_name, char_val)
    return os.path.join(sim_dir, f"sel_{sel_coef}", data_dir_name)

def get_path_to_log_file(char_name, char_val, sel_coef, scenario, ind):
    sim_char_name, sim_char_val = get_sim_char_name_val(char_name, char_val)
    return os.path.join(get_data_dir(sim_char_name, sim_char_val, sel_coef), "logs", scenario, f"{ind}.log")

def get_path_to_vcf_file(char_name, char_val, sel_coef, scenario, ind):
    sim_char_name, sim_char_val = get_sim_char_name_val(char_name, char_val)
    return os.path.join(get_data_dir(sim_char_name, sim_char_val, sel_coef), "vcfs", scenario, str(ind), "merged.vcf")

def get_correct_n_grid(char_name, char_val):
    if char_name != "n_grid":
        return characteristics["n_grid"]["default_value"]
    return char_val

def get_chrom_middle_pos(char_name, char_val):
    "Returns position of middle loci, where selection should have been happened"
    if char_name == "n_loci":
        len_mult = char_val
    else:
        len_mult = characteristics["n_loci"]["default_value"]
    return int(len_mult * 1_000_000 / 2 + 1)

# --- Plotting directories
def get_plotting_base_dir():
    return "figures"

def get_dir_for_master_plots():
    return os.path.join(get_plotting_base_dir(), "master_plots")

def get_filename_to_save_master_plot(char_name):
    return os.path.join(get_dir_for_master_plots(), f"{char_name}_master_plot.pdf")

def get_dir_for_mcc_plots():
    return os.path.join(get_plotting_base_dir(), "power_plots")

def get_filename_to_save_mcc_plot(char_name):
    return os.path.join(get_dir_for_mcc_plots(), f"{char_name}_power_plot.pdf")

def get_dir_for_boxplot_plots():
    return os.path.join(get_plotting_base_dir(), "boxplot_plots")

def get_filename_to_save_boxplot_plot(char_name):
    return os.path.join(get_dir_for_boxplot_plots(), f"{char_name}_boxplot_plot.pdf")
