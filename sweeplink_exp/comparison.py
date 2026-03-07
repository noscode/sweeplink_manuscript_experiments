import os
from . import config

def is_finished(tool_name, ind, char_name=None, char_val=None, is_comparison=False):
    if char_name is None:
        char_name = config.get_char_name_for_comp()
        char_val = config.get_char_val_for_comp()
    if not is_comparison:
        assert tool_name == "sweeplink"
    
    out_dir = config.get_inference_dir_for_val(char_name, char_val, tool_name, is_comparison=is_comparison)
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
    return "results.txt"

def get_files_to_check(tool_name):
    if tool_name == "sweeplink":
        return ["sweepLink_meta.txt", get_tool_input_file(tool_name)]
    return [get_tool_input_file(tool_name)]
