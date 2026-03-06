import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import printing, config

parser = argparse.ArgumentParser()
parser.add_argument('--char', required=True)
parser.add_argument('--n_total', type=int, default=100)
args = parser.parse_args()
#process_and_print_results_per_char_value(args.char, os.path.join(config.BASE_EXP_DIR, args.char), args.n_total)
printing.process_and_print_results_per_sel_coef(args.char, args.n_total)
