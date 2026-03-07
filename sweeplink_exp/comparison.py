import os
from . import config

def is_finished(char_name, char_val, ind):
    tool_name = config.get_tool_name(char_name, char_val)
    
    out_dir = config.get_inference_dir_for_val(char_name, char_val)
    run_dir = os.path.join(out_dir, str(ind))

    if tool_name == "sweeplink":
        return config.is_finished_sweeplink(run_dir)
    if tool_name == "approxwf":
        return is_finished_approxwf(run_dir)
    if tool_name == "diplolocus":
        return is_finished_diplolocus(run_dir)
    assert tool_name == "bmws"
    return is_finished_bmws(run_dir)

def is_finished_approxwf(directory):
    check_file = os.path.join(directory, get_tool_output_file("approxwf"))
    if not os.path.exists(check_file):
        return False
    with open(check_file, "r") as f:
        num_lines = len(f.readlines())
    return num_lines == 100002

def is_finished_diplolocus(directory):
    check_file = os.path.join(directory, get_tool_output_file("diplolocus"))
    if not os.path.exists(check_file):
        return False
    with open(check_file, "r") as f:
        num_lines = len(f.readlines())
    return num_lines > 20

def is_finished_bmws(directory):
    check_file = os.path.join(directory, get_tool_output_file("bmws"))
    if not os.path.exists(check_file):
        return False
    with open(check_file, "r") as f:
        num_lines = len(f.readlines())
    return num_lines > 20

def get_tool_input_file(tool_name):
    if tool_name == "sweeplink":
        return "sweepLink_alleleCounts.txt"
    if tool_name == "approxwf":
        return "approxwf.loci"
    if tool_name == "diplolocus":
        return "diplolocus_alleleCounts.txt"
    assert tool_name == "bmws"
    return "prepared_data.vcf"

def get_tool_output_file(tool_name):
    if tool_name == "sweeplink":
        return "sweepLink_s_statePosteriors.txt"
    if tool_name == "approxwf":
        return "approxwf_MCMC_output.txt"
    if tool_name == "diplolocus":
        return "diplolocus_result_off-grid_maxLLs.txt"
    assert tool_name == "bmws"
    return "result.txt"

def get_files_to_check(tool_name):
    if tool_name == "sweeplink":
        return ["sweepLink_meta.txt", get_tool_input_file(tool_name)]
    return [get_tool_input_file(tool_name)]

def get_generation_func_meta_for_tool(tool_name):
    if tool_name != "sweeplink":
        return None
    def sweeplink_func(time_points):
        # write meta_file
        meta_data = [f"time_{p}_pop_0" for p in time_points]
        lines = ""
        for p in time_points: 
            lines += f"time_{p}_pop_0\t{p}\tpop0\n"
        return get_files_to_check("sweeplink")[0], lines
    return sweeplink_func

def get_counts_filename(tool_name):
    input_files = get_files_to_check(tool_name)
    if tool_name == "sweeplink":
        return input_files[1]
    return input_files[0]

def get_generation_func_counts_header(tool_name):
    def sweeplink_func(time_points):
        meta_data =[f"time_{p}_pop_0" for p in time_points]
        return f"-\t-" + "".join(["\t" + m for m in meta_data])
    def approxwf_func(time_points):
        return "time\t" + "\t".join([str(time) for time in time_points])
    def diplolocus_func(time_points):
        line = "ID"
        for i in range(len(time_points)):
            line += f"\td{i+1}\tn{i+1}"
        return line
    def bmws_func(time_points):
        vcf_path = config.get_path_to_vcf_file(
            char_name="comparison",
            char_val="bmws",
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
    tool2func = {
        "sweeplink": sweeplink_func,
        "approxwf": approxwf_func,
        "diplolocus": diplolocus_func,
        "bmws": bmws_func,
    }
    return tool2func[tool_name]

def get_generation_func_counts_line(tool_name):
    def sweeplink_func(line, sel, scenario, i_chrom, pos, ac, tc):
        line = f"{config.get_chrom_name(i_chrom + 1, sel, scenario)}\t{pos}\t"
        line += "\t".join(f"{_1}/{_2}" for _1, _2 in zip(ac, tc)) + "\n"
        return line
    def approxwf_func(line, sel, scenario, i_chrom, pos, ac, tc):
        if scenario == "neut" or pos == config.get_chrom_middle_pos(config.get_char_name_for_comp(), config.get_char_val_for_comp()):
            line = f"{config.get_chrom_name(i_chrom + 1, sel, scenario)}_{pos}\t"
            line += "\t".join(f"{_1}/{_2}" for _1, _2 in zip(ac, tc)) + "\n"
            return line
        return ""
    def diplolocus_func(line, sel, scenario, i_chrom, pos, ac, tc):
        line = f"{config.get_chrom_name(i_chrom + 1, sel, scenario)}_{pos}\t"
        line += "\t".join(f"{_1}\t{_2}" for _1, _2 in zip(ac, tc)) + "\n"
        return line
    def bmws_func(line, sel, scenario, i_chrom, pos, ac, tc):
        if scenario == "neut" or pos == config.get_chrom_middle_pos(config.get_char_name_for_comp(), config.get_char_val_for_comp()):
            line = line.split()
            line[0] = f"{config.get_chrom_name(i_chrom + 1, sel, scenario)}"
            return "\t".join(line)
        return ""
    tool2func = {
        "sweeplink": sweeplink_func,
        "approxwf": approxwf_func,
        "diplolocus": diplolocus_func,
        "bmws": bmws_func,
    }
    return tool2func[tool_name]

def get_sweeplink_command(char_name, char_val, time_points):
    numGrid = str(config.get_correct_n_grid(char_name, char_val))
    meta_file, input_file = get_files_to_check("sweeplink")
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
            f"loci={get_tool_input_file('approxwf')}",
            "h=0.5",
            "mutRate=1e-8",
            "verbose=1",
            "N=10000",
        ]

def get_diplolocus_command(char_name, char_val, time_points):
    return [
            "DiploLocus", "likelihood",
            "--infile", get_tool_input_file("diplolocus"),
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
            get_tool_input_file("bmws"),
            "meta_file.meta",
            "-l", "4.5",
            "-d", "diploid",
            "-n", "10000",
            "-g", "1",
            "-o", "result.txt",
        ]

def get_generation_func_cmd(tool_name):
    tool2func = {
        "sweeplink": get_sweeplink_command,
        "approxwf": get_approxwf_command,
        "diplolocus": get_diplolocus_command,
        "bmws": get_bmws_command,
    }
    return tool2func[tool_name]
