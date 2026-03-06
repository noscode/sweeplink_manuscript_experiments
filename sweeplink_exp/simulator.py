import os
import shutil
import subprocess
from . import config

def generate_simulation_files(char_name, exp_dir, slim_template_path):
    char_values = config.get_char_config(char_name)['values']
    defaults = config.get_defaults()
    sim_base_dir = config.get_simulation_dir(char_name)
   
    if config.get_char_config(char_name)["can_reuse_simulation"]:
        # does not matter what values to use
        char_values = [char_values[0]]

    for char_val in char_values:
        sel_list = config.get_sim_sel_list(char_name, char_val)
        
        for sel in sel_list:
            # e.g., simulations/true_s_0.01/configs/
            sim_dir = os.path.join(sim_base_dir, f"true_s_{sel}", "configs")
            os.makedirs(sim_dir, exist_ok=True)
            shutil.copy(slim_template_path, os.path.join(sim_dir, "timesweeper_model.slim"))
            
            params = defaults.copy()
            params[char_name] = char_val
            
            tpoints_num, inds_per_tp = params['tpoints_num'], params['sample_size']
            sample_sizes_list =[inds_per_tp] * tpoints_num
            phys_len = int(params['n_loci'] * 1_000_000)
            sample_gens_list =[i * params['binning'] for i in range(tpoints_num)]

            data_dir = config.get_data_dir(char_name, char_val)

            yaml_content = f"""#General
work dir: ../{data_dir}
slimfile: timesweeper_model.slim 

scenarios: ["neut", "sdn"]
mut types: [2]
num_sample_points: {tpoints_num}
sample sizes: {sample_sizes_list}
inds_per_tp: {inds_per_tp}
ploidy: 2
win_size: 51
physLen: {phys_len}
selCoeff: {sel}
sampleGens: {sample_gens_list}
startFreq: {params['start_freq']}
recRate: {params['rec_rate']}

#Simulation
reps: 100
slim path: slim
"""
            with open(os.path.join(sim_dir, f"example_config_val_{char_val}.yaml"), "w") as f:
                f.write(yaml_content)
