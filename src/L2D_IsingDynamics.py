from parsers.L2D_IsingDynamics import parse_arguments, parser
from kernels.L2D_IsingDynamics import *

def main():
    args = parse_arguments(parser)
    ic_gs = args.init_cond.startswith('ground_state')
    number = int(args.init_cond.split('_')[-1]) if ic_gs else 0
    out_suffix = get_out_suffix(args, ic_gs, number)
    l2dDictArgs = initialize_l2d_dict_args(args)
    isingDictArgs = initialize_ising_dict_args(args, out_suffix)
    run_simulation(args, l2dDictArgs, isingDictArgs, number, args.remove_files)
    #
    if args.print_chrono:
        Chronometer.print_all_chronometers()
    if args.verbose:
        print("Done!")

if __name__ == "__main__":
    main()
