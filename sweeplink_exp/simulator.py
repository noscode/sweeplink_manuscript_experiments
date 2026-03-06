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

def generate_simulation_slurm(char_name, n_cores=32, submit=False):
    char_values = config.get_char_config(char_name)['values']
    experiment_dir = config.get_experiment_dir(char_name)
    slurm_path = os.path.join(experiment_dir, f"run_simulations_{char_name}.sh")
    
    with open(slurm_path, "w") as f:
        f.write("#!/bin/bash -l\n")
        f.write(f"#SBATCH --error=simulation_{char_name}.err\n#SBATCH --output=simulation_{char_name}.out\n")
        f.write(f"#SBATCH --mem=80000\n#SBATCH --time=100:0:0\n")
        f.write(f"#SBATCH --job-name sim_{char_name}\n#SBATCH --cpus-per-task={n_cores}\n")
        f.write("#SBATCH --partition=pibu_el8\n\n")
        f.write(f"NTHREADS={n_cores}\nconda activate timesweeper_env\n\n")
        
        for char_val in char_values:
            sel_list = config.get_sim_sel_list(char_name, char_val)
            for sel in sel_list:
                f.write(f"cd timesweeper_sims/true_s_{sel}/configs\n")
                f.write(f"timesweeper sim_custom --threads $NTHREADS -y example_config_val_{char_val}.yaml\n")
                f.write("cd ../../..\n\n")

    if submit:
        print(f"\nSubmitting simulation job to the cluster (using {n_cores} cores)...")
        subprocess.run(["sbatch", f"--exclude={config.SLURM_EXCLUDE}", os.path.basename(slurm_path)], cwd=experiment_dir)
