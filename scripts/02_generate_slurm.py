import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, slurm

parser = argparse.ArgumentParser()
parser.add_argument('--char', required=True)
parser.add_argument('--n_start', type=int, default=0)
parser.add_argument('--n_sim', type=int, default=20)
parser.add_argument('--n_cores', type=int, default=20)
parser.add_argument('--submit', action='store_true')
args = parser.parse_args()

exp_dir = os.path.join(config.BASE_EXP_DIR, args.char)
slurm.generate_slurm(args.char, exp_dir, args.n_start, args.n_sim, args.n_cores, submit=args.submit)
