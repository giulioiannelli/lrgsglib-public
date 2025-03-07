from parsers.L2D_IsingDynamics import parse_arguments, parser
from kernels.L2D_IsingDynamics import *

def main():
    args = parse_arguments(parser)
    ic_gs = args.init_cond.startswith('ground_state') or args.init_cond.startswith('gs')
    number = int(args.init_cond.split('_')[-1]) if ic_gs else 0
    out_suffix = get_out_suffix(args)
    run_simulation(args, number, args.remove_files, out_suffix)
    #
    if args.print_chrono:
        Chronometer.print_all_chronometers()
    if args.verbose:
        print("Done!")

if __name__ == "__main__":
    main()
