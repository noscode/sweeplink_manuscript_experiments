import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, results

parser = argparse.ArgumentParser(description="Check intermediate trace files.")
parser.add_argument('--char', required=True)
parser.add_argument('--val', required=True, help="Value of characteristic")
parser.add_argument('--sim', required=True, help="Simulation index")
args = parser.parse_args()

sim_dir = os.path.join(config.BASE_EXP_DIR, args.char, args.val, args.sim)
results.check_intermediate(sim_dir)
