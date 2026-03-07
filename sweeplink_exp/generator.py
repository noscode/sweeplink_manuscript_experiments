import os
import numpy as np
import subprocess
from . import config, comparison

def get_timepoints_from_simulation(char_name, char_val, sel_scenario_el, ind):
    """
    Returns timepoints sorted from past to present (forward in time)
    """
    log_path = config.get_path_to_log_file(
        char_name=char_name,
        char_val=char_val,
        sel_coef=sel_scenario_el["true_s"],
        scenario=sel_scenario_el["scenario"],
        ind=ind,
    )
    with open(log_path) as f:
        points = [int(x) for x in next(f).split("mpGens='c(")[1].split(")'")[0].split(",")]

    assert points[0] == 0
    return points

def read_counts_generator(vcf_file, n_timepoints, n_sample):
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
            if n_sample is None: n_sample = int(len(vcf_data) / n_timepoints)

            ac_list = []  # allele counts per timepoint
            tc_list = []  # total counts per timepoint
            for i in range(n_timepoints):
                ac, tc = 0, 0
                for el in vcf_data[i*n_sample : (i+1)*n_sample]:
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


def generate_files(inference_dir, char_name, char_val, ind, tp_func, counts_filename, get_header_func, get_ac_line_func, cmd_func, is_comparison=False):
    vcf_combinations = config.get_sel_scenario_combinations(is_comparison=is_comparison)

    # 1. Get time points from first combination
    time_points = get_timepoints_from_simulation(
        char_name=char_name,
        char_val=char_val,
        sel_scenario_el=vcf_combinations[0],
        ind=ind,
    )
    if tp_func is not None:
        tp_func(time_points)

    # 2. Create allele counts file
    counts_file = os.path.join(inference_dir, str(ind), counts_filename)
    with open(counts_file, "w") as f:
        if get_header_func is not None:
            f.write(get_header_func(time_points).strip() + "\n")

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
        with open(counts_file, "a") as f:
            for line, _chr, pos, ac, tc in read_counts_generator(vcf_path, n_timepoints=len(time_points), n_sample=None):
                if 0 < np.sum(ac) < np.sum(tc):
                    written += 1
                    sel, scenario = comb["true_s"], comb["scenario"]
                    line = get_ac_line_func(line, sel, scenario, i_chrom, pos, ac, tc).strip()
                    if line!= "":
                        f.write(line + "\n")
        assert written > 0

    # 3. Generate cmd file
    cmd_file = os.path.join(inference_dir, str(ind), "cmd")
    with open(cmd_file, "w") as f:
        cmd = cmd_func(char_name, char_val, time_points)
        f.write(" ".join(cmd))
    subprocess.check_call(['chmod', '+x', cmd_file])


def generate_sweeplink_files(char_name, char_val, ind, is_comparison=False):
    inference_dir = config.get_inference_dir_for_val(char_name, char_val, "sweeplink", is_comparison)
    os.makedirs(os.path.join(inference_dir, str(ind)), exist_ok=True)
    def tp_func(time_points):
        # write meta_file
        meta_data =[f"time_{p}_pop_0" for p in time_points]
        with open(os.path.join(inference_dir, str(ind), "sweepLink_meta.txt"), "w") as f:
            for p in time_points: f.write(f"time_{p}_pop_0\t{p}\tpop0\n")

    def get_header_func(time_points):
        meta_data =[f"time_{p}_pop_0" for p in time_points]
        return f"-\t-" + "".join(["\t" + m for m in meta_data])

    def get_ac_line(line, sel, scenario, i_chrom, pos, ac, tc):
        line = f"{config.get_chrom_name(i_chrom + 1, sel, scenario)}\t{pos}\t"
        line += "\t".join(f"{_1}/{_2}" for _1, _2 in zip(ac, tc)) + "\n"
        return line

    generate_files(
        inference_dir=inference_dir,
        char_name=char_name,
        char_val=char_val,
        ind=ind,
        tp_func=tp_func,
        counts_filename="sweepLink_alleleCounts.txt",
        get_header_func=get_header_func,
        get_ac_line_func=get_ac_line,
        cmd_func=get_sweeplink_command,
        is_comparison=is_comparison,
    )


def generate_approxwf_files(char_name, char_val, ind):
    inference_dir = config.get_inference_dir_for_val(char_name, char_val, tool_name="approxwf",is_comparison=True)
    os.makedirs(os.path.join(inference_dir, str(ind)), exist_ok=True)

    def get_header_func(time_points):
        return "time\t" + "\t".join([str(time) for time in time_points])

    def get_ac_line(line, sel, scenario, i_chrom, pos, ac, tc):
        if scenario == "neut" or pos == config.get_chrom_middle_pos(config.get_char_name_for_comp(), config.get_char_val_for_comp()):
            line = f"{config.get_chrom_name(i_chrom + 1, sel, scenario)}_{pos}\t"
            line += "\t".join(f"{_1}/{_2}" for _1, _2 in zip(ac, tc)) + "\n"
            return line
        return ""

    generate_files(
        inference_dir=inference_dir,
        char_name=char_name,
        char_val=char_val,
        ind=ind,
        tp_func=None,
        counts_filename="approxwf.loci",
        get_header_func=get_header_func,
        get_ac_line_func=get_ac_line,
        cmd_func=get_approxwf_command,
        is_comparison=True,
    )

def generate_diplolocus_files(char_name, char_val, ind):
    inference_dir = config.get_inference_dir_for_val(char_name, char_val, tool_name="diplolocus", is_comparison=True)
    os.makedirs(os.path.join(inference_dir, str(ind)), exist_ok=True)

    def get_header_func(time_points):
        line = "ID"
        for i in range(len(time_points)):
            line += f"\td{i+1}\tn{i+1}"
        return line

    def get_ac_line(line, sel, scenario, i_chrom, pos, ac, tc):
        line = f"{config.get_chrom_name(i_chrom + 1, sel, scenario)}_{pos}\t"
        line += "\t".join(f"{_1}\t{_2}" for _1, _2 in zip(ac, tc)) + "\n"
        return line

    generate_files(
        inference_dir=inference_dir,
        char_name=char_name,
        char_val=char_val,
        ind=ind,
        tp_func=None,
        counts_filename="approxwf.loci",
        get_header_func=get_header_func,
        get_ac_line_func=get_ac_line,
        cmd_func=get_diplolocus_command,
        is_comparison=True,
    )

def generate_bmws_files(char_name, char_val, ind):
    inference_dir = config.get_inference_dir_for_val(char_name, char_val, tool_name="bmws", is_comparison=True)
    os.makedirs(os.path.join(inference_dir, str(ind)), exist_ok=True)

    def get_header_func(time_points):
        vcf_path = config.get_path_to_vcf_file(
            char_name=char_name,
            char_val=char_val,
            sel_coef=0.01,
            scenario="sdn",
            ind=0,
        )
        assert(os.path.exists(vcf_path))
        header = ""
        with open(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    header += line
                break
        return header

    def get_ac_line(line, sel, scenario, i_chrom, pos, ac, tc):
        if scenario == "neut" or pos == config.get_chrom_middle_pos(config.get_char_name_for_comp(), config.get_char_val_for_comp()):
            line = line.split()
            line[0] = f"{config.get_chrom_name(i_chrom + 1, sel, scenario)}"
            return "\t".join(line)
        return ""

    generate_files(
        inference_dir=inference_dir,
        char_name=char_name,
        char_val=char_val,
        ind=ind,
        tp_func=None,
        counts_filename="data.vcf",
        get_header_func=get_header_func,
        get_ac_line_func=get_ac_line,
        cmd_func=get_bmws_command,
        is_comparison=True,
    )

def generate_files_for_comparison(tool_name, ind):
    char_name = config.get_char_name_for_comp()
    char_val = config.get_char_val_for_comp()
    if tool_name == "sweeplink":
        return generate_sweeplink_files(char_name, char_val, ind, is_comparison=True)
    if tool_name == "diplolocus":
        return generate_diplolocus_files(char_name, char_val, ind)
    if tool_name == "approxwf":
        return generate_approxwf_files(char_name, char_val, ind)
    assert tool_name == "bmws"
    return generate_bmws_files(char_name, char_val, ind)
