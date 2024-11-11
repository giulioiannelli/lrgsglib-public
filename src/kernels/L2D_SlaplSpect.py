from lrgsglib.core import *

def process_eigen_distribution(
        fname_base, initial_data_fn, update_data_fn, save_data_fn, args):
    """
    Manages the overall process of computing and saving the eigenvalue or
    eigenvector distributions.
    """
    if args.verbose:
        print(f"Starting process_eigen_distribution for {fname_base}")

    # Check for existing files and determine the number of averages already done
    navg_done = 0
    existing_files = sorted(glob.glob(f"{args.workDir}{fname_base}_*.pkl"))
    if existing_files:
        navg_done = max(int(os.path.splitext(f.split('_')[-1])[0])
                        for f in existing_files)

    # Load existing data or compute initial data if no files exist
    if navg_done > 0:
        fname = f"{args.workDir}{fname_base}_{navg_done}.pkl"
        with open(fname, "rb") as f:
            bin_counter = pk.load(f)
        if args.mode == "eigvec_dist":
            initial_data = [val for counter in bin_counter for val in counter.elements()]
        else:
            initial_data = list(bin_counter.elements())
        if args.verbose:
            print(f"Loaded existing data from {fname}")
    else:
        fname = f"{args.workDir}{fname_base}_{navg_done}.pkl"
        initial_data = initial_data_fn(args)
        if not initial_data:
            raise ValueError("Initial data is empty. Please check the input parameters or data generation.")
        if args.mode == "eigvec_dist":
            bin_counter = [Counter() for _ in range(args.howmany)]
        else:
            bin_counter = Counter()

    # Create symmetric logarithmic bins
    bins, bin_centers = create_symmetric_log_bins(
        np.min(initial_data), np.max(initial_data), args.bins_count)
    nAvgNeed = args.navg - navg_done
    total_periods = (nAvgNeed // args.period) + bool(nAvgNeed % args.period)

    # Process each period until reaching the required number of averages
    for current_period in range(total_periods):
        if args.verbose:
            print(f"Processing period {current_period + 1}/{total_periods}")
        batch_size = min(nAvgNeed - current_period * args.period, args.period)
        bin_counter = update_data_fn(batch_size, bins, bin_centers,
                                     bin_counter, args)
        save_data_fn(args, bin_counter, fname)
        navg_done += batch_size
        new_fname = f"{args.workDir}{fname_base}_{navg_done}.pkl"
        os.rename(fname, new_fname)
        fname = new_fname

    # Rename the final file
    fname_final = f"{args.workDir}{fname_base}_{navg_done}.pkl"
    os.rename(fname, fname_final)
    if args.verbose:
        print(f"Renamed final file to {fname_final}")

def eigvec_initial_data(args):
    """Computes the initial eigenvector data."""
    if args.verbose:
        print("Computing initial eigenvector data...")
    result = np.abs(eigV_for_lattice2D_ptch(
        side=args.side, pflip=args.p, geo=args.geo, mode=args.eigen_mode,
        howmany=args.howmany))
    if args.verbose:
        print(f"Computed initial eigenvector data of size {result.shape}")
    return result

def eigval_initial_data(args):
    """Computes the initial eigenvalue data."""
    if args.verbose:
        print("Computing initial eigenvalue data...")
    result = eigv_for_lattice2D(
        side=args.side, pflip=args.p, geo=args.geo)
    if args.verbose:
        print(f"Computed initial eigenvalue data of size {result.shape}")
    return result

def eigvec_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """Updates the bin counters for eigenvector values."""
    if args.verbose:
        print(f"Updating eigenvector data for batch size {batch_size}...")
    eig_values = [[] for _ in range(args.howmany)]
    for _ in range(batch_size):
        eigV = eigV_for_lattice2D_ptch(
            side=args.side, pflip=args.p, geo=args.geo, mode=args.eigen_mode,
            howmany=args.howmany)
        for i in range(args.howmany):
            eig_values[i].append(eigV[i])
    eig_values = [np.concatenate(i) for i in eig_values]

    min_val = min(np.min(eig) for eig in eig_values)
    max_val = max(np.max(eig) for eig in eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = create_symmetric_log_bins(
            min(min_val, bins[0]), max(max_val, bins[-1]), args.bins_count)

    for i in range(args.howmany):
        bin_counter[i].update(bin_eigenvalues(eig_values[i], bins, bin_centers))
    if args.verbose:
        print("Updated eigenvector data.")
    return bin_counter

def eigval_update_data(batch_size, bins, bin_centers, bin_counter, args):
    """Updates the bin counters for eigenvalue values."""
    if args.verbose:
        print(f"Updating eigenvalue data for batch size {batch_size}...")
    eig_values = [eigv_for_lattice2D(
        side=args.side, pflip=args.p, geo=args.geo) for _ in range(batch_size)]
    eig_values = np.concatenate(eig_values)

    min_val = np.min(eig_values)
    max_val = np.max(eig_values)
    if min_val < bins[0] or max_val > bins[-1]:
        bins, bin_centers = create_symmetric_log_bins(
            min(min_val, bins[0]), max(max_val, bins[-1]), args.bins_count)

    bin_counter.update(bin_eigenvalues(eig_values, bins, bin_centers))
    if args.verbose:
        print("Updated eigenvalue data.")
    return bin_counter

def save_data(args, data, fname):
    """Saves the computed data to a file."""
    if args.verbose:
        print(f"Saving data to {fname}...")
    with open(fname, "wb") as f:
        pk.dump(data, f)
    if args.verbose:
        print(f"Data saved to {fname}")