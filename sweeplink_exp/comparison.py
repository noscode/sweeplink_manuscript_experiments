import os
from . import config

def is_finished(tool_name, ind, char_name=None, char_val=None, is_comparison=False):
    if char_name is None:
        char_name = "sample_size"
        char_val = config.get_char_config()["default_value"]
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
    check_file = os.path.join(directory, "approxwf_MCMC_output.txt")
    if not os.path.exists(check_file):
        return False
    with open(check_file, "r") as f:
        num_lines = len(f.readlines())
    return num_lines == 100002

def is_finished_diplolocus(directory):
    check_file = os.path.join(directory, "diplolocus_result_off-grid_maxLLs.txt")
    if not os.path.exists(check_file):
        return False
    with open(check_file, "r") as f:
        num_lines = len(f.readlines())
    return num_lines > 20

def is_finished_diplolocus(directory):
    check_file = os.path.join(directory, "results.txt")
    if not os.path.exists(check_file):
        return False
    with open(check_file, "r") as f:
        num_lines = len(f.readlines())
    return num_lines > 20
