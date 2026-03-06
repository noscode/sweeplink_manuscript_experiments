import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, comparison

parser = argparse.ArgumentParser()
parser.add_argument('--char', required=True)
parser.add_argument('--char_val', required=True)
args = parser.parse_args()

char_val_parsed = config.parse_char_val(args.char, args.char_val)
comp_dir = os.path.join(config.BASE_EXP_DIR, args.char, "comparison_inference", f"val_{char_val_parsed}")
comparison.print_diplolocus_results(comp_dir)
comparison.print_approxwf_results(comp_dir)
