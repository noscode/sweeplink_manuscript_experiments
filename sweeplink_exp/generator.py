import os
import numpy as np
import subprocess
from . import config, comparison


def get_timepoints_and_nsamples_from_simulation(char_name, char_val, sel_scenario_el, ind):
    """
    Returns timepoints sorted from past to present (forward in time) and list of sample size
    for each timepoint
    """
    log_path = config.get_path_to_log_file(
        char_name=char_name,
        char_val=char_val,
        sel_coef=sel_scenario_el["true_s"],
        scenario=sel_scenario_el["scenario"],
        ind=ind,
    )
    with open(log_path) as f:
        line = next(f)
        points = [int(x) for x in line.split("mpGens='c(")[1].split(")'")[0].split(",")]
        samples_per_step = int(line.split('sampleSizePerStep=')[1].split()[0])
    assert points[0] == 0

    n_sample_list = [samples_per_step for _ in points]

    if char_name != "binning":
        return points, n_sample_list

    # If we have binning experiment we need to bin our samples
    # The binning size is stored in char_val and we take the same time span
    bin_size = char_val
    assert points[-1] % bin_size == 0  # otherwise we cannot do binning
    new_time_points = [i * bin_size for i in range(points[-1] // bin_size + 1)]

    # now we will create and fill number of samples per each point
    # it will be the same in the middle bins, but it will be less in first and last bins
    new_n_samples = [0 for _ in new_time_points]
    half_step = bin_size / 2.0
    for old_num, old_time in zip(n_sample_list, points):
        for i in range(len(new_time_points)):
            if new_time_points[i] - half_step < old_time < new_time_points[i] + half_step:
                new_n_samples[i] += old_num
            # If we are on the border we split samples: half goes to one bin
            # and another half to the next bin
            if old_time == new_time_points[i] - half_step:
                new_n_samples[i] += int(old_num / 2)
            if old_time == new_time_points[i] + half_step:
                new_n_samples[i] += old_num - int(old_num / 2)
    assert(np.sum(n_sample_list) == np.sum(new_n_samples))
    return new_time_points, new_n_samples

def read_counts_generator(vcf_file, n_sample_list):
    repeated_positions = []
    cur_pos = 0
    with open(vcf_file) as g:
        for line in g:
            if line.startswith("#"):
                continue
            pos = int(line.split()[1].strip())
            if cur_pos == pos:
                repeated_positions.append(pos)
            cur_pos = pos

    with open(vcf_file) as g:
        for line in g:
            if line.startswith("#"):
                continue
            chrom, cur_pos = int(line.split()[0].strip()), int(line.split()[1].strip())
            assert chrom == 1
            if cur_pos in repeated_positions:  # multiallelic
                continue
            pos = cur_pos

            vcf_data = line.strip().split()[9:]
            assert len(vcf_data) == np.sum(n_sample_list)

            ac_list = []  # allele counts per timepoint
            tc_list = []  # total counts per timepoint
            for i, n_samples in enumerate(n_sample_list):
                ac, tc = 0, 0
                start_ind = int(np.sum(n_sample_list[:i]))
                for el in vcf_data[start_ind : start_ind + n_samples]:
                    if el[0] == "1": ac += 1
                    if el[0] in ["0", "1"]: tc += 1
                    if el[2] == "1": ac += 1
                    if el[2] in ["0", "1"]: tc += 1
                ac_list.append(ac)
                tc_list.append(tc)
            yield (line, chrom, pos, ac_list, tc_list)


def get_sweeplink_command(char_name, char_val, time_points):
    numGrid = str(config.get_correct_n_grid(char_name, char_val))
    meta_file, input_file = comparison.get_files_to_check("sweeplink")
    return [
            config.SWEEPLINK_PATH, "infer",
            "--meta", meta_file,
            "--counts", input_file,
            "--nGridPoints", numGrid,
            "--numThreads", "1",
            "--numBurnin", "10",
            "--h.update", "false",
            "--numStatesAbsS", str(len(config.S_POS_GRID)),
            "--grid_abs_s", ",".join([str(s) for s in config.S_POS_GRID]),
            "--forwardPrior", "uniform",
            "--newMutPrior", "neutral",
            "--mu_a_A", "1e-8",
            "--mu_A_a", "0",
            "--N", "10000",
            "--chebychevGrid",
            "--cutNonSegregating", "1,1",
        ]

def get_approxwf_command(char_name, char_val, time_points):
    return [
            config.APPROXWF_PATH,
            "task=estimate",
            f"loci={comparison.get_tool_input_file('approxwf')}",
            "h=0.5",
            "mutRate=1e-8",
            "verbose=1",
            "N=10000",
        ]

def get_diplolocus_command(char_name, char_val, time_points):
    return [
            "DiploLocus", "likelihood",
            "--infile", comparison.get_tool_input_file("diplolocus"),
            "--sample_times", ",".join([str(_x) for _x in time_points]),
            "--fix_h", "0.5",
            '--fix_s2=' + f'{",".join([str(s) for s in sorted(list(set(config.SWEEPLINK_GRID)))])}',
            "--u01", "1e-8",
            "--u10", "0.0",
            "--Ne", "10000",
            "--init", "uniform",
            "--o", "diplolocus_result",
            "--get_off_grid_max",
            "--get_MLR",
            "--get_chi2_pval",
        ]

def get_bmws_command(char_name, char_val, time_points):
    return [
            "bmws", "analyze",
            comparison.get_tool_input_file("bmws"),
            "meta_file.meta",
            "-l", "4.5",
            "-d", "diploid",
            "-n", "10000",
            "-g", "1",
            "-o", "result.txt",
        ]


def generate_files(char_name, char_val, ind):
    vcf_combinations = config.get_sel_scenario_combinations(char_name)
    base_dir = os.path.join(config.get_inference_dir_for_val(char_name, char_val), str(ind))
    os.makedirs(base_dir, exist_ok=True)

    # 1. Get time points from first combination
    time_points, n_sample_list = get_timepoints_and_nsamples_from_simulation(
        char_name=char_name,
        char_val=char_val,
        sel_scenario_el=vcf_combinations[0],
        ind=ind,
    )
    # for sweeplink we have to create meta file
    generate_meta_file = config.get_generation_func_meta(char_name, char_val)
    if generate_meta_file is not None:
        filename, lines = generate_meta_file(time_points)
        with open(os.path.join(base_dir, filename), "w") as f:
            f.write(lines)

    # 2. Create allele counts file
    generate_counts_file_header = config.get_generation_func_counts_header(char_name, char_val)
    generate_counts_file_line = config.get_generation_func_counts_line(char_name, char_val)

    counts_filename = config.get_counts_filename(char_name, char_val)
    counts_path = os.path.join(base_dir, counts_filename)
    with open(counts_path, "w") as f:
        if generate_counts_file_header is not None:
            f.write(generate_counts_file_header(time_points).strip() + "\n")

    for i_chrom, comb in enumerate(vcf_combinations):
        vcf_path = config.get_path_to_vcf_file(
            char_name=char_name,
            char_val=char_val,
            sel_coef=comb["true_s"],
            scenario=comb["scenario"],
            ind=ind,
        )
        assert(os.path.exists(vcf_path))
        written = 0
        with open(counts_path, "a") as f:
            for line, _chr, pos, ac, tc in read_counts_generator(vcf_path, n_sample_list=n_sample_list):
                if 0 < np.sum(ac) < np.sum(tc):
                    sel, scenario = comb["true_s"], comb["scenario"]
                    line = generate_counts_file_line(line, sel, scenario, i_chrom, pos, ac, tc).strip()
                    if line!= "":
                        f.write(line + "\n")
                        written += 1
        assert written > 0

    # 3. Generate cmd file
    generate_cmd_file = config.get_generation_func_cmd(char_name, char_val)
    cmd_file = os.path.join(base_dir, "cmd")
    with open(cmd_file, "w") as f:
        cmd = generate_cmd_file(char_name, char_val, time_points)
        f.write(" ".join(cmd))
    subprocess.check_call(['chmod', '+x', cmd_file])
