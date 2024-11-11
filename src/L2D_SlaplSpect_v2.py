from parsers.L2D_SlaplSpect import *
from kernels.L2D_SlaplSpect import *

def main():
    """Main function to parse arguments and run the process."""
    args = parser.parse_args()
    args.workDir = Lattice2D(
        side1=args.L, pflip=args.p, geo=args.geometry, sgpathn=args.workDir
    ).spectpath

    # Determine the filename base based on the mode
    if args.mode == "eigvec_dist":
        fname_base = f"dist{args.howmany}_{args.p:.3g}_{args.eigen_mode}"
        initial_fn = eigvec_initial_data
        update_fn = eigvec_update_data
    else:
        fname_base = f"dist_eigval_{args.p:.3g}_{args.cell_type}"
        initial_fn = eigval_initial_data
        update_fn = eigval_update_data

    # Process the eigen distribution
    process_eigen_distribution(fname_base, initial_fn, update_fn, save_data,
                               args)

if __name__ == "__main__":
    main()
