import os
import subprocess
from . import config

def generate_slurm(char_name, exp_dir, n_start, n_sim, n_cores, submit=False):
    char_values = config.get_char_config(char_name)['values']
    inference_base = os.path.join(exp_dir, "sweeplink_inference")
    run_files_dir = os.path.join(inference_base, "run_files")
    os.makedirs(os.path.join(run_files_dir, "output"), exist_ok=True)

    pending_jobs =[]

    for char_val in char_values:
        for sim in range(n_sim):
            actual_sim_idx = n_start + sim
            sim_dir = os.path.join(inference_base, f"val_{char_val}", str(actual_sim_idx))
            sL_file = os.path.join(sim_dir, "sweepLink_s_statePosteriors.txt")
            
            should_run = True
            if os.path.exists(sL_file):
                try:
                    with open(sL_file) as g:
                        if len(next(g).strip()) > 0: should_run = False
                except: pass
                    
            if should_run:
                pending_jobs.append(f"cd {sim_dir}\n./cmd > output.log")

    if not pending_jobs:
        print(f"All {len(char_values) * n_sim} SweepLink jobs are already completed!")
        return

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
