import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, simulator

parser = argparse.ArgumentParser()
parser.add_argument('--char', required=True)
parser.add_argument('--n_cores', type=int, default=32)
parser.add_argument('--submit', action='store_true')
args = parser.parse_args()

slim_template = os.path.join(os.path.dirname(__file__), "..", "templates", "timesweeper_model.slim")
simulator.generate_simulation_files(args.char, slim_template)
simulator.generate_simulation_slurm(args.char, n_cores=args.n_cores, submit=args.submit)
