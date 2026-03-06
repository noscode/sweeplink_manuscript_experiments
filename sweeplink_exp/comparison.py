import os
import subprocess
import numpy as np
from scipy.stats import chi2
from . import config

def setup_comparison_files(sweeplink_base_dir, comp_base_dir):
    os.makedirs(comp_base_dir, exist_ok=True)
    diplo_dir, approx_dir = os.path.join(comp_base_dir, "diplolocus"), os.path.join(comp_base_dir, "approxwf")
    
    for ind in range(config.COMPARE_N_SIM):
        sl_counts = os.path.join(sweeplink_base_dir, str(ind), "sweepLink_alleleCounts.txt")
        if not os.path.exists(sl_counts): continue
            
        with open(sl_counts) as f:
            time_points =[int(el.split("_")[1]) for el in next(f).split()[2:]]

        d_out, a_out = os.path.join(diplo_dir, str(ind)), os.path.join(approx_dir, str(ind))
        os.makedirs(d_out, exist_ok=True); os.makedirs(a_out, exist_ok=True)

        with open(os.path.join(d_out, "diplolocus_alleleCounts.txt"), "w") as fd, open(os.path.join(a_out, "approxwf.loci"), "w") as fa:
            fd.write("ID" + "".join([f"\td{i+1}\tn{i+1}" for i in range(len(time_points))]) + "\n")
            fa.write("time\t" + "\t".join(map(str, time_points)) + "\n")
            
            with open(sl_counts) as g:
                next(g)
                for line in g:
                    p = line.split()
                    eval_s = p[0].split("_")[1] # 1_s_0.01 -> s_0.01
                    # DiploLocus
                    fd.write("\t".join([f"{p[0]}_{p[1]}"] + [item for el in p[2:] for item in el.split("/")]) + "\n")
                    # ApproxWF
                    if eval_s != "sel_0.0" and int(p[1]) == 50001:
                        fa.write("\t".join([f"{p[0]}_{p[1]}"] + p[2:]) + "\n")
                        
        with open(os.path.join(d_out, "cmd"), "w") as fd:
            fix_s2 = "-0.16666666666666669,-0.13043478260869565,-0.1111111111111111,-0.09090909090909091,-0.07407407407407407,-0.056603773584905655,-0.047619047619047616,-0.038461538461538464,-0.029126213592233007,-0.0196078431372549,-0.01477832512315271,-0.01234567901234568,-0.009900990099009901,-0.007936507936507936,-0.005964214711729622,-0.0049751243781094535,0,0.005,0.006,0.008,0.01,0.0125,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.125,0.15,0.2"
            fd.write(f"DiploLocus likelihood --infile diplolocus_alleleCounts.txt --sample_times {','.join(map(str, time_points))} --fix_h 0.5 --fix_s2=\"{fix_s2}\" --u01 1e-8 --u10 1e-8 --Ne 10000 --init uniform --o diplolocus_result --get_off_grid_max --get_MLR --get_chi2_pval")
        subprocess.check_call(['chmod', '+x', os.path.join(d_out, "cmd")])

        with open(os.path.join(a_out, "cmd"), "w") as fa:
            fa.write(f"{config.APPROXWF_BIN_PATH} task=estimate loci=approxwf.loci h=0.5 mutRate=1e-8 verbose=1 N=10000")
        subprocess.check_call(['chmod', '+x', os.path.join(a_out, "cmd")])

def generate_comparison_slurm(char_name, char_val, comp_base_dir, n_cores, submit=False):
    run_files_dir = os.path.join(comp_base_dir, "run_files")
    os.makedirs(os.path.join(run_files_dir, "output"), exist_ok=True)
    pending_jobs =[]

    for ind in range(config.COMPARE_N_SIM):
        for tool, out_name in[("diplolocus", "diplolocus_result_off-grid_maxLLs.txt"), ("approxwf", "approxwf_meanVar.txt")]:
            out_dir = os.path.join(comp_base_dir, tool, str(ind))
            if os.path.exists(out_dir):
                out_file = os.path.join(out_dir, out_name)
                if not os.path.exists(out_file) or os.stat(out_file).st_size == 0:
                    pending_jobs.append(f"cd {out_dir}\n./cmd > out.log")

    if not pending_jobs:
        print("All comparison jobs completed!")
        return

    actual_cores = min(n_cores, len(pending_jobs))
    jobs_per_core = [[] for _ in range(actual_cores)]
    for idx, job in enumerate(pending_jobs): jobs_per_core[idx % actual_cores].append(job)

    generated_scripts =[]
    for core_idx, jobs in enumerate(jobs_per_core):
        if not jobs: continue
        index = f"comp_core{core_idx}"
        filename = os.path.join(run_files_dir, f'run_{index}.sh')
        generated_scripts.append(filename)
        with open(filename, "w") as f:
            f.write(f"#!/bin/bash -l\n#SBATCH --error=output/{index}.err\n#SBATCH --output=output/{index}.out\n#SBATCH --mem=80000\n#SBATCH --time=100:0:0\n#SBATCH --job-name {char_name}_cmp_{char_val}_{core_idx}\n#SBATCH --cpus-per-task=1\n#SBATCH --partition=pibu_el8\n\n")
            for job in jobs: f.write(job + "\n\n")

    if submit:
        for script in generated_scripts: subprocess.run(["sbatch", f"--exclude={config.SLURM_EXCLUDE}", os.path.basename(script)], cwd=run_files_dir)

def print_diplolocus_results(comp_base_dir):
    print("\n--- DiploLocus Results ---")
    bounds, confusion, results_pval = config.get_s_grid_bounds(), {c: {s: 0 for s in config.S_GRID} for c in config.COMPARE_EVAL_COEFS}, {c:[] for c in config.COMPARE_EVAL_COEFS}
    for ind in range(config.COMPARE_N_SIM):
        fpath = os.path.join(comp_base_dir, "diplolocus", str(ind), "diplolocus_result_off-grid_maxLLs.txt")
        if not os.path.exists(fpath): continue
        with open(fpath) as f:
            try: next(f)
            except: continue
            for line in f:
                if line.startswith("#") or line.startswith("ID"): continue
                parts = line.split()
                sel_coef, pos = float(parts[0].split("_")[2]), int(parts[0].split("_")[3])
                pval = chi2.sf(float(parts[-2]), df=1)
                if sel_coef == 0.0 or pos == 50001:
                    results_pval[sel_coef].append(pval)
                    mean_s_val = float(parts[3])
                    if pval < 0.0004 and mean_s_val > 0:
                        for i, b in enumerate(bounds):
                            if b[0] <= mean_s_val < b[1]:
                                confusion[sel_coef][config.S_GRID[i]] += 1
                                break
    for sel in config.COMPARE_EVAL_COEFS:
        res = np.array(results_pval[sel])
        if len(res) == 0: continue
        print(f"\nSelection: {sel}\n  Runs with p-value < 0.0004: {np.sum(res < 0.0004)}/{len(res)} ({int(np.sum(res < 0.0004)/len(res) * 100)}%)")
        print("  Confusion matrix: | " + " ".join([str(confusion[sel][s]) + (" |" if s==0.0 else "") for s in config.S_GRID]))

def print_approxwf_results(comp_base_dir):
    print("\n--- ApproxWF Results ---")
    bounds, confusion = config.get_s_grid_bounds(), {c: {s: 0 for s in config.S_GRID} for c in config.COMPARE_EVAL_COEFS}
    results_p, results_mean = {c:[] for c in config.COMPARE_EVAL_COEFS}, {c:[] for c in config.COMPARE_EVAL_COEFS}
    for ind in range(config.COMPARE_N_SIM):
        fpath = os.path.join(comp_base_dir, "approxwf", str(ind), "approxwf_meanVar.txt")
        if not os.path.exists(fpath): continue
        with open(fpath) as f:
            try: next(f)
            except: continue
            for line in f:
                parts = line.strip().split()
                if not parts[0].startswith("s"): continue
                sel_coef, pos = float(parts[0].split("_")[2]), int(parts[0].split("_")[3])
                neut_prob, pos_prob, neg_prob, mean_s_val = float(parts[1]), float(parts[2]), float(parts[3]), float(parts[2])
                if sel_coef == 0.0 or pos == 50001:
                    results_p[sel_coef].append(max(pos_prob, neg_prob))
                    results_mean[sel_coef].append(mean_s_val)
                    if neut_prob >= 0.01: mean_s_val = 0
                    for i, b in enumerate(bounds):
                        if b[0] <= mean_s_val < b[1]:
                            confusion[sel_coef][config.S_GRID[i]] += 1
                            break
    for sel in config.COMPARE_EVAL_COEFS:
        res, res2 = np.array(results_p[sel]), np.array(results_mean[sel])
        if len(res) == 0: continue
        print(f"\nSelection: {sel}\n  Runs with >99% probability: {np.sum(res > 0.99)}/{len(res)} ({int(np.sum(res > 0.99)/len(res) * 100)}%)")
        print(f"  Mean probability: {np.mean(res):.4f} +- {np.std(res):.4f}\n  Mean pred selection: {np.mean(res2):.5f}")
        print("  Confusion matrix: | " + " ".join([str(confusion[sel][s]) + (" |" if s==0.0 else "") for s in config.S_GRID]))
