import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, generator, slurm

parser = argparse.ArgumentParser()
parser.add_argument('--char', required=True)
parser.add_argument('--val', type=str, default=None)
parser.add_argument('--n_start', type=int, default=0)
parser.add_argument('--n_sim', type=int, default=20)
parser.add_argument('--n_cores', type=int, default=20)
parser.add_argument('--submit', action='store_true')
parser.add_argument('--skip_generation', action='store_true')
args = parser.parse_args()

all_values = config.get_char_config(args.char)["values"]
if args.val is not None:
    d = {str(val): val for val in all_values}
    args.val = d[args.val]
    val_list = [args.val]
else:
    val_list = all_values

if not args.skip_generation:
    for char_val in val_list:
        print(f"Generating for {args.char}={char_val} (Sims {args.n_start} to {args.n_start+args.n_sim-1})")
        for i in range(args.n_start, args.n_start + args.n_sim):
            generator.generate_files(args.char, char_val, i)

slurm.generate_inference_slurm(
    char_name=args.char,
    char_values=val_list,
    n_start=args.n_start,
    n_sim=args.n_sim,
    n_cores=args.n_cores,
    submit=args.submit
)
