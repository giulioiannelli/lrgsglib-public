from parsers.L2D_SlaplSpect import parse_arguments, parser
from kernels.L2D_SlaplSpect import *

def main():
    args = parse_arguments(parser)
    perform_spectral_calculations(args)
    #
    if args.print_chrono:
        Chronometer.print_all_chronometers()
    if args.verbose:
        print("Done!")

if __name__ == "__main__":
    main()



# from parsers.L2D_SlaplSpect import *
# from kernels.L2D_SlaplSpect import *

# def main():
    # """Main function to parse arguments and run the process."""
    # args = parse_arguments(parser)
    # args.workDir = Lattice2D(
    #     side1=args.L, pflip=args.p, geo=args.geo, path_data=args.workDir
    # ).path_spect
    # # Determine the filename base based on the mode
    # if args.mode == "eigvec_dist":
    #     fname_base = f"dist{args.howmany}_{args.p:.3g}_{args.eigen_mode}"
    #     initial_fn = eigvec_initial_data
    #     update_fn = eigvec_update_data
    # elif args.mode == "eigval_dist":
    #     fname_base = f"dist_eigval_{args.p:.3g}_{args.cell_type}"
    #     initial_fn = eigval_initial_data
    #     update_fn = eigval_update_data
    # elif args.mode == "eigvals":
    #     fname_base = f"eigvals_{args.p:.3g}_{args.cell_type}"
    #     initial_fn = eigval_initial_data
    #     update_fn = eigval_update_data
    # # Process the eigen distribution
    # process_eigen_distribution(fname_base, initial_fn, update_fn, save_data,
    #                            args)

# if __name__ == "__main__":
#     main()
