import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, generator

parser = argparse.ArgumentParser()
parser.add_argument('--tools', type=str, default="sweeplink,approxwf")
parser.add_argument('--n_start', type=int, default=0)
parser.add_argument('--n_sim', type=int, default=20)
args = parser.parse_args()

for tool_name in args.tools.split(","):
    print(f"Generating for {tool_name} (Sims {args.n_start} to {args.n_start+args.n_sim-1})")
    for i in range(args.n_start, args.n_start + args.n_sim):
        generator.generate_files_for_comparison(tool_name, i)
