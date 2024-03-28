from .shared import *
from .const import *
from .errwar import *
#
def sym_log_func_unsafe(x, a, b, c, d):
    """Return values from a general log function with safety checks."""
    # Ensure safe computation by adjusting values close to d
    return a * np.log(b * (np.abs(x) - d)) + c
def sym_log_func(x, a, b, c, d):
    """Return values from a general log function with safety checks."""
    # Ensure safe computation by adjusting values close to d
    safe_x = np.where(np.abs(x) > d, np.abs(x), d + 1e-10)
    return a * np.log(b * (safe_x - d)) + c

def adjust_to_even(x):
    """
    Rounds the input number to the nearest even integer.

    If the input is exactly halfway between two even numbers, it rounds up to the
    higher even number.

    Parameters
    ----------
    x : float
        The input number to be rounded.

    Returns
    -------
    int
        The nearest even integer to the input number.

    Examples
    --------
    >>> round_to_nearest_even(128 * np.sqrt(3))
    222
    >>> round_to_nearest_even(5.5)
    6
    >>> round_to_nearest_even(2.1)
    2
    """
    lower_even = int(x) - int(x) % 2
    upper_even = lower_even + 2
    return lower_even if x - lower_even < upper_even - x else upper_even

def is_in_range(number, range_start, range_end):
    """
    Checks if a given number is within a specified range, inclusive.

    Parameters:
    -----------
    number : int or float
        The number to check if it lies within the given range.
    range_start : int or float
        The starting value of the range (inclusive).
    range_end : int or float
        The ending value of the range (inclusive).

    Returns:
    --------
    bool
        Returns True if the number is within the range specified by
        range_start and range_end, inclusive; otherwise, False.

    Notes:
    ------
    This function uses Python's ability to chain comparison operators,
    making the check concise and efficient. It is versatile enough to handle
    both integer and floating-point numbers.

    Example:
    --------
    >>> is_in_range(5, 1, 10)
    True
    >>> is_in_range(15, 1, 10)
    False
    """
    return range_start <= number <= range_end


def find_matching_files(search_dir: str, pattern_str: str) -> List[str]:
    """
    Searches all files in the specified directory that match a given pattern.

    Parameters
    ----------
    search_dir : str
        The directory in which to search for files.
    pattern_str : str
        The string pattern to search for within the file names. This pattern
        will be compiled into a regular expression.

    Returns
    -------
    List[str]
        A list of file names within the search directory that match the given 
        pattern.
    """
    from os import listdir
    from re import compile
    pattern = compile(f'.*{pattern_str}.*')
    all_files = listdir(search_dir)
    matching_files = [fname for fname in all_files if pattern.match(fname)]
    return matching_files

def flatten(xs):
    """
    Recursively flattens a nested iterable into a flat iterable.

    Parameters:
    -----------
    xs : iterable
        The input nested iterable to be flattened.

    Returns:
    --------
    object
        The flattened elements from the input iterable.

    Notes:
    ------
    This function takes a nested iterable and yields each element from the
    nested structure as a flat iterable. It recursively processes nested
    iterables such as lists or tuples.

    Example:
    --------
    >>> nested_list = [1, [2, [3, 4], 5], 6]
    >>> flattened = list(flatten(nested_list))
    >>> flattened
    [1, 2, 3, 4, 5, 6]
    """
    for x in xs:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            yield from flatten(x)
        else:
            yield x

def move_to_rootf(print_tf: bool = False):
    """
    Move to the root directory of the current working directory.

    Parameters:
    -----------
    print_tf : bool, optional
        If True, print the current working directory after moving to the root.
        Default is False.

    Notes:
    ------
    - The function continuously moves up the directory hierarchy ('../') until it reaches
      a directory with the name specified by the 'PATH_ROOTF' constant.
    - If 'print_tf' is set to True, it prints the current working directory after the move.

    Example:
    --------
    To move to the root directory of the current working directory and print the path:
    move_to_rootf(print_tf=True)
    """
    while getcwd()[-len(PATH_ROOTF):] != PATH_ROOTF:
        chdir('../')
    if print_tf:
        print('cwd:', getcwd())

def boolean_overlap_fraction(boolist1, boolist2):
    """
    Calculate the fraction of overlapping True values between two boolean lists.

    Parameters:
    -----------
    boolist1 : list of bool
        The first boolean list.

    boolist2 : list of bool
        The second boolean list.

    Returns:
    --------
    float
        The fraction of overlapping True values between the two boolean lists.

    Notes:
    ------
    - The function computes the overlap fraction by first performing a bitwise XOR (^)
      operation between the two boolean lists, which results in a new boolean list where
      True represents differences and False represents matches.
    - The bitwise NOT (~) operator is then applied to invert the differences,
      turning them into True values.
    - The 'sum' function counts the number of True values in the inverted list,
      representing the number of overlapping True values.
    - Finally, this count is divided by the length of 'boolist1' to obtain the overlap fraction.

    Example:
    --------
    boolist1 = [True, False, True, False, True]
    boolist2 = [True, True, False, False, True]
    overlap_fraction = boolean_overlap_fraction(boolist1, boolist2)
    # The result is 0.4, indicating 40% overlap of True values between the two lists.

    """
    return sum(~(boolist1 ^ boolist2))/len(boolist1)

def first_index_changing_condition(condition):
    """
    Find the index of the first change in a boolean condition.

    Parameters:
    -----------
    condition : numpy.ndarray of bool
        A boolean condition represented as a NumPy array.

    Returns:
    --------
    int
        The index of the first change in the condition.

    Notes:
    ------
    - The function compares adjacent elements in the 'condition' array and returns
      the index of the first change from True to False or vice versa.
    - If there are no changes in the condition, the function returns 0.

    Example:
    --------
    condition = np.array([True, True, True, False, False, True])
    first_change_index = first_index_changing_condition(condition)
    # The result is 3, indicating the first change from True to False occurred at index 3.

    """
    return np.where(condition[:-1] != condition[1:])[0][0]

# basic math functions
def line(x, a, b):
    """
    Calculate the values of a straight line equation for given 'x' values.

    Parameters:
    -----------
    x : float or array-like
        The input values of 'x' for which the corresponding 'y' values are calculated.

    a : float
        The slope of the straight line.

    b : float
        The y-intercept of the straight line.

    Returns:
    --------
    float or numpy.ndarray
        The calculated 'y' values for the given 'x' values based on the straight line equation.

    Notes:
    ------
    - The straight line equation is 'y = a * x + b', where 'a' is the slope and 'b' is the y-intercept.
    - The function computes 'y' values for single 'x' values or arrays of 'x' values.

    Example:
    --------
    x_value = 2.0
    slope = 1.5
    y_intercept = 2.0
    y_result = line(x_value, slope, y_intercept)
    # The result is 5.0, which is the 'y' value for 'x' = 2.0 in the given line equation.

    """
    return a * x + b

#
def dv(f_x: ndarray, x: ndarray = None) -> ndarray:
    """
    Compute the derivative of an array with respect to another array.

    Parameters:
    -----------
    f_x : ndarray
        A NumPy array of N dimensions where the outermost dimension is the one
        where the derivative is computed using the `np.diff` method.

    x : ndarray, optional
        An array with the same dimensions as the outermost axis of 'f_x'.
        If not provided, it is assumed to be the range [0, len(f_x)].

    Returns:
    --------
    df_dx : ndarray
        The computed derivative of 'f_x' with respect to 'x'.

    Notes:
    ------
    - This function calculates the derivative of 'f_x' with respect to 'x' using
      the finite difference method. It computes the derivative along the outermost axis.
    - If 'x' is not provided, it is assumed to be the range [0, len(f_x)].

    Example:
    --------
    f_x = np.array([[1, 2, 3], [4, 5, 6]])
    x = np.array([0, 1, 2])
    df_dx = dv(f_x, x)
    # The result 'df_dx' contains the computed derivatives along the outermost axis.

    """
    if x is None:
        x = np.linspace(0, f_x.shape[-1], num=f_x.shape[-1])
    df_dx = np.diff(f_x, axis=-1) / np.diff(x)
    return df_dx

def lin_binning(data, binnum=20):
    """
    Perform linear binning on data and calculate the histogram.

    Parameters:
    -----------
    data : numpy.ndarray or list
        The input data for which the histogram is to be computed.

    binnum : int, optional
        The number of bins to use for the histogram. Default is 20.

    Returns:
    --------
    bin_centers : numpy.ndarray
        The center values of each bin in the histogram.

    hist : numpy.ndarray
        The histogram values for each bin.

    Notes:
    ------
    - Linear binning divides the data range into 'binnum' equal-width bins and
      calculates the histogram based on the bin counts.
    - The function returns the center values of the bins and the corresponding histogram values.

    Example:
    --------
    data = np.array([1, 2, 2, 3, 4, 4, 4, 5, 6, 7])
    bin_centers, hist = lin_binning(data, binnum=5)
    # The result provides the center values and histogram counts for 5 bins.

    """

    min_val = int(np.floor(np.min(data)))
    max_val = int(np.ceil(np.max(data)))
    bins = np.linspace(min_val, max_val, num=binnum)
    hist, _ = np.histogram(data, bins=bins)
    bin_centers = (bins[1:] + bins[:-1]) / 2.0
    return bin_centers, hist

def log_binning(data, binnum=20):
    """
    Perform logarithmic binning on data and calculate the histogram.

    Parameters:
    -----------
    data : numpy.ndarray or list
        The input data for which the histogram is to be computed.

    binnum : int, optional
        The number of bins to use for the histogram. Default is 20.

    Returns:
    --------
    bin_centers : numpy.ndarray
        The center values of each bin in the histogram.

    hist : numpy.ndarray
        The histogram values for each bin.

    bin_w : numpy.ndarray
        The bin widths for each bin in the histogram.

    Notes:
    ------
    - Logarithmic binning divides the logarithmic range of the data into 'binnum' bins and
      calculates the histogram based on the bin counts.
    - The function returns the center values of the bins, the corresponding histogram values,
      and the bin widths.

    Example:
    --------
    data = np.array([1, 10, 100, 1000, 10000])
    bin_centers, hist, bin_w = log_binning(data, binnum=5)
    # The result provides the center values, histogram counts, and bin widths for 5 log bins.

    """
    log_data = np.log10(data)
    min_val = int(np.floor(np.min(log_data)))
    max_val = int(np.ceil(np.max(log_data)))
    bins = np.logspace(min_val, max_val, num=binnum)
    hist, _ = np.histogram(data, bins=bins)
    bin_w = (bins[1:] - bins[:-1])
    bin_centers = (bins[1:] + bins[:-1]) / 2.0
    return bin_centers, hist, bin_w

def neglog_binning(data, binnum=20):
    """
    Perform negative logarithmic binning on data and calculate the histogram.

    Parameters:
    -----------
    data : numpy.ndarray or list
        The input data for which the histogram is to be computed.

    binnum : int, optional
        The number of bins to use for the histogram. Default is 20.

    Returns:
    --------
    bin_centers : numpy.ndarray
        The center values of each bin in the histogram.

    hist : numpy.ndarray
        The histogram values for each bin.

    bin_w : numpy.ndarray
        The bin widths for each bin in the histogram.

    Notes:
    ------
    - Negative logarithmic binning divides the negative logarithmic range of the absolute data
      into 'binnum' bins and calculates the histogram based on the bin counts.
    - The function returns the center values of the bins, the corresponding histogram values,
      and the bin widths.
    - Bins are reversed to have a negative x-scale for visualization.

    Example:
    --------
    data = np.array([-1, -10, -100, -1000, -10000])
    bin_centers, hist, bin_w = neglog_binning(data, binnum=5)
    # The result provides the center values, histogram counts, and bin widths for 5 neglog bins.

    """
    abs_data = np.abs(data)
    log_data = np.log10(abs_data)
    min_val = int(np.floor(np.min(log_data)))
    max_val = int(np.ceil(np.max(log_data)))
    bins = -np.logspace(min_val, max_val, num=binnum)[::-1]  # Reverse the bins for negative x-scale
    hist, _ = np.histogram(data, bins=bins)
    bin_w = (bins[1:] - bins[:-1])
    bin_centers = (bins[1:] + bins[:-1]) / 2.0
    return bin_centers, hist, bin_w


def symlog_binning(full_data, binnum=20):
    """
    Perform symmetric logarithmic binning on positive and negative data separately
    and calculate histograms.

    Parameters:
    -----------
    full_data : numpy.ndarray or list
        The input data containing both positive and negative values.

    binnum : int, optional
        The number of bins to use for the histograms. Default is 20.

    Returns:
    --------
    outp : tuple or None
        A tuple containing the result of log_binning for positive data or None
        if there are no positive values in the input.

    outm : tuple or None
        A tuple containing the result of neglog_binning for negative data or None
        if there are no negative values in the input.

    Notes:
    ------
    - The function splits the input data into positive and negative components
      and performs logarithmic binning separately on each component.
    - If there are no positive or negative values, the corresponding result is None.
    - 'outp' and 'outm' are tuples containing the center values, histogram values,
      and bin widths for positive and negative data, respectively.

    Example:
    --------
    full_data = np.array([1, 10, -0.1, -100, 10000])
    outp, outm = symlog_binning(full_data, binnum=5)
    # The result provides log_binning and neglog_binning results for positive and negative data.

    """
    datap = full_data[full_data > 0]
    datam = full_data[full_data < 0]
    if datap.size:
        outp = log_binning(datap)
    else:
        outp = None
    if datam.size:
        outm = neglog_binning(datam)
    else:
        outm = None
    return outp, outm


def round_sigfig_n(num, n: int = 1):
    """
    Round a number or array of numbers to a specified number of significant figures.

    Parameters:
    -----------
    num : float or array-like
        The number or array of numbers to be rounded to 'n' significant figures.

    n : int, optional
        The desired number of significant figures (default is 1).
        Must be within the range [1, 15].

    Returns:
    --------
    float or numpy.ndarray
        The rounded number or array of numbers with 'n' significant figures.

    Raises:
    -------
    ValueError
        If 'n' is not within the valid range [1, 15].

    Notes:
    ------
    - The function calculates the exponent required to obtain 'n' significant figures
      based on the absolute value of 'num'.
    - It then applies the rounding operation to 'num' with the calculated exponent
      to achieve the desired number of significant figures.
    - If 'num' is an array-like object, it processes each element separately.

    Example:
    --------
    num = 123.456789
    n = 3
    rounded_num = round_sigfig_n(num, n)
    # The result is 123.0, rounded to 3 significant figures.

    num_array = [123.456789, 0.00123456789]
    rounded_array = round_sigfig_n(num_array, n)
    # The result is [123.0, 0.00123], with each element rounded to 3 significant figures.

    """
    if n not in range(1, DEFAULT_MAX_DIGITS_ROUND_SIGFIG):
        raise ValueError("Significant figures number not in [1, 15].")
    expn = -np.floor(np.log10(np.abs(num))).astype('int')
    if hasattr(num, "__len__"):
        rr = np.array([np.round(nn, ee+n-1) for nn,ee in zip(num, expn)])
    else:
        rr = np.round(num, expn+n-1)
    return rr

def width_interval(a, b):
    """
    Calculate the width of an interval between two values.

    Parameters:
    -----------
    a : float or numeric
        The first value defining the interval.

    b : float or numeric
        The second value defining the interval.

    Returns:
    --------
    float or numeric
        The width of the interval, which is the absolute difference between 'a' and 'b'.

    Example:
    --------
    a = 5
    b = 8
    interval_width = width_interval(a, b)
    # The result is 3, representing the width of the interval [5, 8].

    """
    return np.abs(a - b)


def randstring(size=10, chars=string.ascii_lowercase + string.ascii_uppercase + string.digits):
    return ''.join(np.random.choice(list(chars), size))


def inf_array_regularization(arrinfs):
    """
    Regularizes an array by replacing infinite values with extremes of non-infinite values.

    Parameters:
    - arrinfs (ndarray): Input array containing infinite and finite numerical values.

    Returns:
    - ndarray: Regularized array with infinite values replaced by extremes of non-infinite values.

    Examples:
    >>> arr = np.array([1, 2, np.inf, 5, -np.inf, 7])
    >>> inf_array_regularization(arr)
    array([1., 2., 7., 5., 7., 7.])
    """
    # Filtering out infinite values from the input array (arrinfs)
    arrinfs_nnans = arrinfs[(arrinfs != np.inf) & (arrinfs != -np.inf)]

    # Regularizing infinite values in the array using nan_to_num function
    arrinfs = np.nan_to_num(arrinfs, posinf=np.max(arrinfs_nnans),
                            neginf=np.min(arrinfs_nnans))

    return arrinfs



def binder_cumulant(data):
    """
    Calculate the Binder cumulant for a set of data.

    Parameters:
    -----------
    data : array_like
        Input data representing measurements of the order parameter.

    Returns:
    --------
    float
        Binder cumulant of the data.

    Notes:
    ------
    The Binder cumulant is a statistical measure used in the analysis of 
    phase transitions in statistical physics. It is defined as 
    1 - (mean(m^4) / (3 * mean(m^2)^2)), where m is the order parameter.
    The Binder cumulant is particularly useful in finite-size scaling 
    analysis as it is dimensionless and often has a universal value at 
    the critical point for systems within the same universality class.

    Examples:
    ---------
    >>> data = np.random.normal(size=1000)
    >>> print(binder_cumulant(data))

    >>> uniform_data = np.random.uniform(-1, 1, size=1000)
    >>> print(binder_cumulant(uniform_data))
    """
    m2 = np.mean(data**2)
    m4 = np.mean(data**4)
    U4 = 1 - m4 / (3 * m2**2)
    return U4

def sum_tuples(tuple1: tuple, tuple2: tuple) -> tuple:
    """
    Sum two tuples element-wise.

    If the tuples are of different lengths, the function sums elements up to the
    length of the shorter tuple, ignoring extra elements in the longer tuple.

    Parameters:
    -----------
    tuple1 : tuple
        The first tuple to be summed.
    tuple2 : tuple
        The second tuple to be summed.

    Returns:
    --------
    tuple
        A new tuple containing the element-wise sums of `tuple1` and `tuple2`.

    Example:
    --------
    >>> sum_tuples((1, 2, 3), (4, 5, 6))
    (5, 7, 9)
    """
    return tuple(a + b for a, b in zip(tuple1, tuple2))
#
def flip_to_positive_majority(arr):
    """ Flips the elements of an array to ensure a majority of positive components.

        Given a numerical array, this function checks if the majority of its components are negative. If so, it multiplies every element by -1, effectively flipping the signs of all components to ensure a majority of positive values. This operation is intended for arrays where elements can be distinctly categorized as positive or negative (including zero as a non-negative value).

        Parameters:
        -----------
        arr : array_like
            The input array containing numerical data. This array can be a list, tuple, or any array-like object convertible to a NumPy array. The function is optimized for NumPy arrays for performance reasons.

        Returns:
        --------
        numpy.ndarray
            An array of the same shape and type as the input, with elements flipped if the original array had a majority of negative components. If the input array already had a majority of positive components, it is returned unchanged.

        Notes:
        ------
        The function utilizes NumPy for efficient computation, especially for large arrays. The determination of "majority" is based purely on the count of positive vs. negative elements, without weighting by magnitude.

        Examples:
        ---------
        >>> import numpy as np
        >>> arr = np.array([1, -2, -3, 4, -5])
        >>> flip_to_positive_majority(arr)
        array([ -1,  2,  3,  -4,  5])

        >>> arr = np.array([-1, 2, 3])
        >>> flip_to_positive_majority(arr)
        array([-1,  2,  3])
        """
    if np.sum(arr < 0) > len(arr) / 2:
        # Flip all elements by multiplying by -1
        arr = arr * -1
    return arr


def shift_with_wrap(image, shift_right, shift_down):
    """
    Shift an image with wrap-around at the edges.
    
    Parameters:
    - image: 2D numpy array to shift.
    - shift_right: Number of pixels to shift to the right.
    - shift_down: Number of pixels to shift down.
    
    Returns:
    - Shifted image with wrap-around.
    """
    # Ensure shifts are within the bounds of the image dimensions
    shift_right %= image.shape[1]
    shift_down %= image.shape[0]
    
    # Perform the shift
    shifted_image = np.roll(image, shift=shift_down, axis=0)  # Shift down
    shifted_image = np.roll(shifted_image, shift=shift_right, axis=1)  # Then shift right
    
    return shifted_image


# def unravel_1d_to_2d_nodemap(arr1d, imap, dims: Tuple = None):
#     if not dims:
#         dims = tuple(int(np.sqrt(len(arr1d))) for i in range(2))
#     arr_2d = np.empty(dims, dtype=arr1d.dtype)
#     for idx_1d, (row, col) in imap.items():
#         arr_2d[col, row] = arr1d[idx_1d]
#     return arr_2d

def unravel_1d_to_2d_nodemap(arr1d: np.ndarray, imap: Dict[int, Tuple[int, int]], dims: Tuple[int, int] = None) -> np.ndarray:
    """
    Transforms a 1D array into a 2D array based on a given index mapping and optional dimensions.

    Parameters:
    ---------------
    arr1d : np.ndarray
        The 1D input numpy array to be transformed.
    imap : Dict[int, Tuple[int, int]]
        A dictionary where keys are indices in the 1D array and values are (row, col) positions in the 2D output.
    dims : Tuple[int, int], optional
        The dimensions of the 2D output array. If not provided, it's assumed to be a square array.

    Returns:
    ---------------
    np.ndarray
        The 2D array obtained from rearranging `arr1d` according to `imap`.
    """
    if not dims:
        side = int(np.sqrt(len(arr1d)))
        dims = (side, side)
    arr_2d = np.empty(dims, dtype=arr1d.dtype)
    for idx_1d, (row, col) in imap.items():
        arr_2d[row, col] = arr1d[idx_1d]  # Adjusted indexing to [row, col] for clarity
    return arr_2d


def unravel_1d_to_2d_nodemap2(arr1d: np.ndarray, imap: Dict[int, Tuple[int, int]], dims: Tuple[int, int] = None) -> np.ndarray:
    """
    Transforms a 1D array into a 2D array based on a given index mapping and optional dimensions.

    Parameters:
    ---------------
    arr1d : np.ndarray
        The 1D input numpy array to be transformed.
    imap : Dict[int, Tuple[int, int]]
        A dictionary where keys are indices in the 1D array and values are (row, col) positions in the 2D output.
    dims : Tuple[int, int], optional
        The dimensions of the 2D output array. If not provided, it's assumed to be a square array.

    Returns:
    ---------------
    np.ndarray
        The 2D array obtained from rearranging `arr1d` according to `imap`.
    """
    if not dims:
        side = int(np.sqrt(len(arr1d)))
        dims = (side, side)
    arr_2d = np.empty(dims, dtype=arr1d.dtype)
    for idx_1d, (row, col) in imap.items():
        arr_2d[col, row] = arr1d[idx_1d]  # Adjusted indexing to [row, col] for clarity
    return arr_2d

def get_first_int_in_str(s):
    """
    Finds the first sequence of digits in the given string and returns it as an integer.
    If no digits are found, returns None.

    Parameters:
    ---------------
    s : str
        The string to search for digits.

    Returns:
    ---------------
    int or None
        The first integer found in the string, or None if no digits are found.
    """
    match = re.search(r'\d+', s)
    return int(match.group()) if match else None


from itertools import product

def cProd_Iter(dim: Union[int, Tuple]) -> iter:
    """
    Generates the Cartesian product for an n-dimensional space.

    This function creates an iterable over all possible combinations of coordinates
    in an n-dimensional grid, defined by the dimensions specified in the input tuple.
    Each element in the 'dim' tuple represents the size of the grid along that dimension.
    The Cartesian product is generated using a sequence of range objects, one for each
    dimension, thereby producing a sequence of tuples that represent points in the
    n-dimensional space.

    Parameters
    ----------
    dim : tuple of int
        A tuple where each element specifies the size of the grid along a particular
        dimension. For example, (2, 3) represents a 2D grid with dimensions 2x3.

    Returns
    -------
    iterator
        An iterator over tuples, where each tuple represents a point in the
        n-dimensional space defined by 'dim'. Each tuple contains integers,
        with the i-th integer representing the coordinate along the i-th dimension.

    Examples
    --------
    >>> list(cartesian_product((2, 3)))
    [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    >>> list(cartesian_product((1, 2, 3)))
    [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), (0, 1, 1), (0, 1, 2)]

    Notes
    -----
    - The order of the tuples in the output iterator follows the lexicographical order
      based on the input dimensions. This means smaller coordinates come before larger ones.
    - This function is particularly useful for iterating over multi-dimensional arrays
      or grids, where you need to visit each cell or point in the space.
    - The implementation relies on itertools.product, which is efficient and avoids
      explicitly constructing the grid in memory, making it suitable for large dimensions.

    See Also
    --------
    itertools.product : Cartesian product of input iterables.
    """
    return product(*[range(d) for d in dim])

def cProdSel_Iter(dim: Union[int, Tuple], selected_indices: Union[List, Tuple]) -> iter:
    """
    Generates the Cartesian product for selected dimensions in an n-dimensional space.

    This function creates an iterable over all possible combinations of coordinates
    in the specified dimensions of an n-dimensional grid. It allows for focusing on
    a subset of all available dimensions, which is specified by the selected_indices.

    Parameters
    ----------
    dim : tuple of int
        A tuple where each element specifies the size of the grid along a particular
        dimension. For example, (2, 3, 4) represents a 3D grid with dimensions 2x3x4.
    selected_indices : list or tuple of int
        Indices of the dimensions to include in the Cartesian product. For example,
        (1, 2) selects the second and third dimensions from 'dim'.

    Returns
    -------
    iterator
        An iterator over tuples, where each tuple represents a point in the
        specified dimensions of the n-dimensional space. Each tuple contains integers,
        with the i-th integer representing the coordinate along the selected i-th dimension.

    Examples
    --------
    >>> list(cartesian_product_selected((2, 3, 4), (0, 2)))
    [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3)]

    >>> list(cartesian_product_selected((2, 3, 4, 5), (1, 3)))
    [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4)]

    Notes
    -----
    - The function is useful for iterating over a selected subset of dimensions in a multi-dimensional space.
    - By allowing specification of dimensions of interest, it provides flexibility for applications that do not require
      the full Cartesian product of all dimensions.
    """
    # Generate ranges for selected dimensions
    ranges = [range(dim[i]) for i in selected_indices]
    
    # Return Cartesian product of selected dimensions
    return product(*ranges)

def cProd_Iter_adj(dim: Union[int, Tuple], range_adjustment: Union[int, List] = 0) -> iter:
    """
    Generates the Cartesian product for an n-dimensional space with adjustable ranges.

    This function creates an iterable over all possible combinations of coordinates
    in an n-dimensional grid, defined by the dimensions specified in the input tuple,
    adjusted by the range_adjustment which can be a single integer or a list of integers.

    Parameters
    ----------
    dim : tuple of int
        A tuple where each element specifies the size of the grid along a particular
        dimension. For example, (2, 3) represents a 2D grid with dimensions 2x3.
    range_adjustment : int or list of int, optional
        An integer to adjust all dimensions equally, or a list of integers to adjust
        each dimension individually. The default value is 0, which means no adjustment.

    Returns
    -------
    iterator
        An iterator over tuples, where each tuple represents a point in the
        adjusted n-dimensional space defined by 'dim' and 'range_adjustment'.
        Each tuple contains integers, with the i-th integer representing the
        coordinate along the i-th dimension.

    Examples
    --------
    >>> list(cProd_Iter((2, 3), 1))
    [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2), (2, 3)]

    >>> list(cProd_Iter((1, 2), [1, 2]))
    [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3)]

    Notes
    -----
    The `range_adjustment` allows for flexible configuration of the generated space,
    accommodating scenarios that require padding or extending the dimensions.
    """
    # Handle both int and list types for range_adjustment
    if isinstance(range_adjustment, int):
        adjusted_ranges = [range(d + range_adjustment) for d in dim]
    elif isinstance(range_adjustment, list):
        adjusted_ranges = [range(d + adj) for d, adj in zip(dim, range_adjustment)]
    else:
        raise ValueError("range_adjustment must be an int or a list of ints")

    return product(*adjusted_ranges)
