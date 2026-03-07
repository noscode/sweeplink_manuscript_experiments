import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import printing, config

parser = argparse.ArgumentParser()
parser.add_argument('--char', required=True)
parser.add_argument('--n_total', type=int, default=100)
args = parser.parse_args()
printing.process_and_print_results_per_sel_coef(args.char, args.n_total)
