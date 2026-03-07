import sys
import os
import argparse
import shutil

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, comparison

def compare_files_robust(f1, f2):
    """Compares two files line by line, stripping whitespace to avoid failing on trailing spaces/newlines."""
    if not os.path.exists(f1) or not os.path.exists(f2):
        print(f1, f2)
        return False
    with open(f1) as file1, open(f2) as file2:
        # Strip newlines and spaces, and ignore completely blank lines
        l1 = [line.strip() for line in file1 if line.strip()]
        l2 = [line.strip() for line in file2 if line.strip()]
    #return True
    if len(l1) != len(l2):
        return False
    return l1 == l2


def migrate(old_val_dir, new_val_dir, info_string, inputs_to_check):
    if not os.path.exists(old_val_dir):
        print(f"[{info_string}] Skipped: Old directory ({old_val_dir}) does not exist.")
        return
        
    migrated_count = 0
    mismatch_count = 0
    missing_new_count = 0
    
    for ind in range(args.n_sim):
        old_ind_dir = os.path.join(old_val_dir, str(ind))
        new_ind_dir = os.path.join(new_val_dir, str(ind))
        
        # If the old directory doesn't have results for this index, skip it
        if not os.path.exists(old_ind_dir):
            continue
            
        # If the new directory hasn't been generated yet, we can't compare
        if not os.path.exists(new_ind_dir):
            missing_new_count += 1
            continue
            
        # 1. Verify Inputs Match
        match = True
        for inp in inputs_to_check:
            old_file = os.path.join(old_ind_dir, inp)
            new_file = os.path.join(new_ind_dir, inp)
            if not compare_files_robust(old_file, new_file):
                match = False
                break
                
        # 2. Copy Outputs if valid
        if match:
            for f in os.listdir(old_ind_dir):
                # Only copy output files (ignore the inputs we just checked, and ignore cmd)
                if f not in inputs_to_check and f != "cmd":
                    old_f_path = os.path.join(old_ind_dir, f)
                    new_f_path = os.path.join(new_ind_dir, f)
                    
                    if os.path.isfile(old_f_path):
                        shutil.copy2(old_f_path, new_f_path)
            migrated_count += 1
        else:
            mismatch_count += 1
            
    print(f"[{info_string}] Migrated: {migrated_count} | Mismatched Inputs: {mismatch_count} | Missing Target Dir: {missing_new_count}")

parser = argparse.ArgumentParser(description="Migrate old sweepLink results to the new architecture.")
parser.add_argument('--char', required=True, help="Name of characteristic")
parser.add_argument('--old_dir', required=True, help="Path to the old root directory (e.g., ../new_round)")
parser.add_argument('--n_sim', type=int, default=100, help="Number of simulations to check")
args = parser.parse_args()

if args.char == "comparison":
    print(f"Starting migration ...\n")
    for tool_name in config.THRESHOLDS:
        if tool_name == "sweeplink":
            continue
        new_dir = config.get_inference_dir_for_val(args.char, "",  tool_name=tool_name, is_comparison=True)
        tool2dir = {"sweeplink": "sweepLink", "approxwf": "approxwf_quick", "diplolocus": "diplolocus", "bmws": "bmws"}
        old_dir = os.path.join(args.old_dir, "compare_tools", tool2dir[tool_name], "results")
        print(f"Starting migration from {old_dir} -> {new_dir} ...\n")
        inputs_to_check = comparison.get_files_to_check(tool_name)

        migrate(old_dir, new_dir, info_string=f"tool {tool_name}", inputs_to_check=inputs_to_check)
    
else:
    char_values = config.get_char_config(args.char)['values']
    new_base_dir = config.get_inference_dir(args.char, tool_name="sweeplink", is_comparison=False)
    old_base_dir = os.path.join(args.old_dir, args.char, "inference")

    inputs_to_check = comparison.get_files_to_check("sweeplink")

    char_config = config.get_char_config(args.char)

    print(f"Starting migration from {old_base_dir} -> {new_base_dir} ...\n")

    for char_val in char_values:
        old_val_dir = os.path.join(old_base_dir, char_config['val2str'][char_val])
        new_val_dir = config.get_inference_dir_for_val(args.char, char_val, tool_name="sweeplink", is_comparison=False)
        migrate(old_val_dir, new_val_dir, info_string=f"{args.char}={char_val}", inputs_to_check=inputs_to_check)

print("\nMigration complete!")
