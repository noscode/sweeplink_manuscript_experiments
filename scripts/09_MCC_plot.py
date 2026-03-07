import sys, os, argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sweeplink_exp import config, plotting

parser = argparse.ArgumentParser()
parser.add_argument('--char', required=True)
parser.add_argument('--n_total', type=int, default=100)
parser.add_argument('--out', type=str, default=None)
args = parser.parse_args()
if args.out is None:
    os.makedirs(config.get_plotting_base_dir(), exist_ok=True)
    os.makedirs(config.get_dir_for_mcc_plots(), exist_ok=True)
    args.out = config.get_filename_to_save_mcc_plot(args.char)

plotting.draw_MCC_plot(args.char, save_file=args.out, n_total=args.n_total)
