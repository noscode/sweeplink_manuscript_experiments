import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, generator

parser = argparse.ArgumentParser()
parser.add_argument('--char', required=True)
parser.add_argument('--n_start', type=int, default=0)
parser.add_argument('--n_sim', type=int, default=20)
args = parser.parse_args()

exp_dir = os.path.join(config.BASE_EXP_DIR, args.char)
for char_val in config.get_char_config(args.char)['values']:
    print(f"Generating for {args.char}={char_val} (Sims {args.n_start} to {args.n_start+args.n_sim-1})")
    for i in range(args.n_start, args.n_start + args.n_sim):
        generator.generate_sweeplink_files(args.char, char_val, i, is_comparison=False)
