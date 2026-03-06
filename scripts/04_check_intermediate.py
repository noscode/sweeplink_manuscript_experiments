import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import printing

parser = argparse.ArgumentParser(description="Check intermediate trace files.")
parser.add_argument('--char', required=True)
parser.add_argument('--val', required=True, help="Value of characteristic")
parser.add_argument('--sim', required=True, help="Simulation index")
args = parser.parse_args()
args.val = int(args.val)

printing.process_and_print_intermediate_results(
    char_name=args.char,
    char_val=args.val,
    ind=args.sim,
    is_comparison=False
)
