import os
import subprocess
from . import config, comparison
import random

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


def generate_inference_slurm(char_name, char_values, n_start, n_sim, n_cores, tool_name, is_comparison=False, submit=False):
    inference_base = config.get_inference_dir(
        char_name=char_name,
        tool_name=tool_name,
        is_comparison=is_comparison
    )
    run_files_dir = os.path.join(inference_base, "run_files")
    os.makedirs(os.path.join(run_files_dir, "output"), exist_ok=True)

    pending_jobs =[]

    for char_val in char_values:
        out_dir = config.get_inference_dir_for_val(
            char_name=char_name,
            char_val=char_val,
            tool_name=tool_name,
            is_comparison=is_comparison
        )
        for sim in range(n_sim):
            actual_sim_idx = n_start + sim
            should_run = not comparison.is_finished(
                tool_name=tool_name,
                ind=actual_sim_idx,
                char_name=char_name,
                char_val=char_val,
                is_comparison=is_comparison
            )
            if should_run:
                sim_dir = os.path.join(out_dir, str(actual_sim_idx))
                pending_jobs.append(f"cd {sim_dir}\n./cmd > output.log")

    if not pending_jobs:
        print(f"All {len(char_values) * n_sim} SweepLink jobs are already completed!")
        return

    random.shuffle(pending_jobs)

    actual_cores = min(n_cores, len(pending_jobs))
    jobs_per_core = [[] for _ in range(actual_cores)]
    for idx, job in enumerate(pending_jobs):
        jobs_per_core[idx % actual_cores].append(job)

    generated_scripts =[]
    print(f"Distributing {len(pending_jobs)} pending jobs across {actual_cores} scripts...")
    
    for core_idx, jobs in enumerate(jobs_per_core):
        if not jobs: continue
        index = f"infer_start{n_start}_core{core_idx}"
        filename = os.path.join(run_files_dir, f'run_{index}.sh')
        generated_scripts.append(filename)
        
        with open(filename, "w") as f:
            f.write("#!/bin/bash -l\n")
            f.write(f"#SBATCH --error=output/{index}.err\n#SBATCH --output=output/{index}.out\n")
            f.write(f"#SBATCH --mem=80000\n#SBATCH --time=100:0:0\n")
            f.write(f"#SBATCH --job-name {char_name}_inf_{core_idx}\n#SBATCH --cpus-per-task=1\n")
            f.write("#SBATCH --partition=pibu_el8\n\n")
            for job in jobs: f.write(job + "\n\n")

    if submit:
        print("\nSubmitting jobs to the cluster...")
        for script in generated_scripts:
            subprocess.run(["sbatch", f"--exclude={config.SLURM_EXCLUDE}", os.path.basename(script)], cwd=run_files_dir)
