from ..shared import *
from .const import *
from .errwar import *
from .tools import *
#
def do_nothing(*args, **kwargs):
    pass
# Global dictionary to store timing data
chronometer_instances = {}

def time_function_accumulate(func):
    def wrapper(*args, **kwargs):
        # Create a new chronometer instance for the function
        chrono = Chronometer(func.__name__)
        chronometer_instances[func.__name__] = chrono
        
        # Execute the function
        result = func(*args, **kwargs)
        
        # End the chronometer
        chrono.end()
        
        return result
    return wrapper

def print_accumulated_timings():
    print("\nSummary of accumulated function timings:")
    for func_name, chrono in chronometer_instances.items():
        elapsed_time = chrono.get_elapsed_time()
        formatted_time = f"{elapsed_time:.4g}"
        print(f"Function '{func_name}' (Chronometer ID: {chrono.id}) elapsed time: {formatted_time} seconds")
def get_numerical_precision(dtype: Type[np.floating] = float) -> float:
    """
    Returns the smallest positive number that can be represented in floating-point arithmetic
    for the specified data type. This value, known as machine epsilon, provides a measure 
    of the numerical precision or the difference between 1 and the smallest floating point 
    number greater than 1 for the given data type.
    
    Args:
        dtype (Type[np.floating]): The floating-point data type to get the precision for 
                                   (e.g., np.float32, np.float64).
                                   Defaults to float (typically equivalent to np.float64).
    
    Returns:
        float: The machine epsilon for the given data type.
    """
    precision = np.finfo(dtype).eps
    print(f"The numerical precision on the resolution of 0 for {dtype} is: {precision}")
    return precision
def join_non_empty(symb: str, *args):
    """
    Joins non-empty strings with an underscore.
    
    Parameters:
        *args: Variable length argument list of strings.
    
    Returns:
        A single string with non-empty strings joined by an underscore.
    """
    non_empty_args = [s for s in args if s]
    
    if len(non_empty_args) <= 1:
        return ''.join(non_empty_args)
    else:
        return symb.join(non_empty_args)

def conditional_print(message, verbose, **kwargs):
    if verbose:
        print(message, **kwargs)
#
def unzip_dict_items(input_dict: Dict[Any, Any]) -> Tuple[List[Any], List[Any]]:
    """
    Unzip a dictionary into two lists containing its keys and values, respectively.

    Parameters:
    ---------------
    input_dict : Dict[Any, Any]
        The dictionary from which to extract keys and values.

    Returns:
    ---------------
    Tuple[List[Any], List[Any]]
        A tuple containing two lists: the first with keys and the second with values from the dictionary.

    Examples:
    ---------------
    >>> test_dict = {'a': 1, 'b': 2, 'c': 3}
    >>> keys, values = unzip_dict_items(test_dict)
    >>> print(keys)
    ['a', 'b', 'c']
    >>> print(values)
    [1, 2, 3]
    """
    keys, values = zip(*input_dict.items()) if input_dict else ([], [])
    return list(keys), list(values)
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

def sign_with_threshold(arr, threshold=1e-17):
    """
    Apply the sign function to an array with a threshold for zero.

    Parameters:
    - arr: NumPy array, the input array.
    - threshold: float, values with absolute value below this are set to 0.

    Returns:
    - A NumPy array where each entry is the sign of the input,
      with values below the threshold set to 0.
    """
    # Create a mask for values below the threshold
    mask = np.abs(arr) < threshold
    # Apply the mask to set these values to 0
    arr[mask] = 0
    # Use np.sign on the modified array
    return np.sign(arr)


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

def sort_array_by_column(arr: np.ndarray, column_index: int) -> np.ndarray:
    """
    Sorts a numpy array by a specific column.

    Parameters
    ----------
    arr : np.ndarray
        The array to be sorted.
    column_index : int
        The index of the column to sort by.

    Returns
    -------
    np.ndarray
        The sorted array.
    """
    return arr[arr[:, column_index].argsort()]


def find_matching_files(search_dir: str, pattern_str: str) -> Union[str, List[str]]:
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
    Union[str, List[str]]
        A single file name if only one file matches the pattern, or a list of 
        file names within the search directory that match the given pattern.
    """
    pattern = re.compile(f'.*{pattern_str}.*')
    all_files = os.listdir(search_dir)
    matching_files = [fname for fname in all_files if pattern.match(fname)]
    
    if len(matching_files) == 1:
        return matching_files[0]
    return matching_files

def extract_value_from_filename(file_name: str, value_pattern: str) -> float:
    """
    Extracts a numeric value from a single file name based on a specified regular expression pattern.
    
    Parameters
    ----------
    file_name : str
        The file name from which to extract the value.
    value_pattern : str
        The regular expression pattern used to extract the numeric value from the file name. The pattern
        should contain a capturing group for the numeric value.
        
    Returns
    -------
    float
        The extracted numeric value from the file name.
        
    Examples
    --------
    >>> file_name = "data_p=1.5.pkl"
    >>> value_pattern = r"p=([\\d.]+)"
    >>> extract_value_from_filename(file_name, value_pattern)
    1.5
    
    This function extracts the p-value from the provided file name.
    """
    match = re.search(value_pattern, file_name)
    if match:
        return float(match.group(1).rstrip('.'))
    else:
        raise ValueError("No match found in the file name.")

def extract_values_from_filenames(file_names: List[str], value_pattern: str, sort: bool = True, unique: bool = False) -> np.ndarray:
    """
    Extracts numeric values from a list of file names based on a specified regular expression pattern.
    
    Parameters
    ----------
    file_names : List[str]
        A list of file names from which to extract the values.
    value_pattern : str
        The regular expression pattern used to extract numeric values from the file names. The pattern
        should contain a capturing group for the numeric value.
    sort : bool, optional
        Specifies whether to sort the extracted values. Default is True.
        
    Returns
    -------
    numpy.ndarray
        An array of extracted numeric values from the file names. If `sort` is True, this array will be
        sorted in ascending order.
        
    Examples
    --------
    >>> file_names = ["data_p=1.5.pkl", "experiment_p=2.0.pkl", "results_p=0.5.pkl"]
    >>> value_pattern = r"p=([\\d.]+)"
    >>> extract_values_from_filenames(file_names, value_pattern)
    array([0.5, 1.5, 2.0])
    
    This function extracts the p-values from the provided file names and returns them sorted in ascending order.
    """
    # values = [re.search(value_pattern, filename).group(1) 
    #           for filename in file_names if re.search(value_pattern, filename)]
    values = [
        float(match.group(1).rstrip('.'))
        for file_name in file_names
        if (match := re.search(value_pattern, file_name))
    ]
    if sort:
        values.sort()
    if unique:
        values = uniques(values)
    return np.array(values)

# def extract_values_from_filenames(filenames, pattern):
#     return [re.search(pattern, filename).group(1) for filename in filenames if re.search(pattern, filename)]

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

def move_to_rootf(print_tf: bool = True):
    """
    Move to the root directory of the current working directory.

    Parameters:
    -----------
    print_tf : bool, optional
        If True, print the current working directory after moving to the root.
        Default is True.

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
    pcwd = getcwd()
    try:
        while getcwd()[-len(PATH_ROOTF)+1:] != PATH_ROOTF[:-1]:
            chdir('../')
            if getcwd() == '/':
                break
        if getcwd() == '/':
            raise FileNotFoundError(f"Root directory '{PATH_ROOTF}' not found.")
        if print_tf:
            print("Current working directory:", getcwd())
    except FileNotFoundError as e:
        chdir(pcwd)
        print(e, "\nCurrent working directory:", pcwd)

def extract_and_sort_values(path: str, search_pattern: str, value_pattern: str = None, sort: bool = True) -> np.ndarray:
    """
    Searches for files in a given directory that match a specified pattern, extracts numerical values from
    these file names based on another pattern (assumed to be the same as search pattern if not provided), 
    and optionally sorts these values. Raises an error if value_pattern does not contain a capturing group.
    
    Parameters
    ----------
    path : str
        The directory path to search for files.
    search_pattern : str
        The regular expression pattern to match file names for the search.
    value_pattern : str, optional
        The regular expression pattern used to extract numeric values from the matched file names. If None, 
        search_pattern is used. This pattern must contain at least one capturing group.
    sort : bool, optional
        Specifies whether to sort the extracted values. Default is True.
        
    Returns
    -------
    numpy.ndarray
        An array of extracted numeric values from the matched file names. If `sort` is True, this array will be
        sorted in ascending order.
        
    Raises
    ------
    ValueError
        If value_pattern does not contain a capturing group.
    """
    if value_pattern is None:
        value_pattern = search_pattern

    # Check if the value_pattern contains at least one capturing group
    if not re.search(r"\((?!\?:).+?\)", value_pattern):
        raise ValueError("value_pattern must contain at least one capturing group.")
    
    file_names = find_matching_files(path, search_pattern)
    return extract_values_from_filenames(file_names, value_pattern, sort)


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
def dv(f_x: NDArray, x: NDArray = None) -> NDArray:
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
def flip_to_positive_majority_adapted(arr: NDArray) -> NDArray:
    """ Flips the elements of an array to ensure a majority of positive components.
    
    This function assesses the balance of positive and negative values in the given numerical array.
    If the count of negative values exceeds the count of positive values, all elements of the array
    are multiplied by -1. This operation ensures that the majority of elements in the transformed array
    are positive. The function is particularly useful in data processing scenarios where the sign of
    data points affects subsequent analysis or visualization.

    Parameters:
    -----------
    arr : ndarray
        The input array containing numerical data. This array can be a list, tuple, or any array-like
        object convertible to a NumPy array. The function is optimized for NumPy arrays for performance reasons.

    Returns:
    --------
    ndarray
        An array of the same shape and type as the input, with elements flipped if the original array
        had a majority of negative components. If the input array already had a majority of positive
        components, it is returned unchanged.

    Examples:
    ---------
    >>> import numpy as np
    >>> arr = np.array([1, -2, -3, 4, -5])
    >>> flip_to_positive_majority_adapted(arr)
    array([-1,  2,  3, -4, 5])

    >>> arr = np.array([-1, 2, 3])
    >>> flip_to_positive_majority_adapted(arr)
    array([-1, 2, 3])
    """
    num_negatives = np.sum(arr < 0)
    num_positives = np.sum(arr > 0)

    if num_negatives > num_positives:
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



def read_files_to_2d_array(folder_path, keyword):
    """
    Read files from a folder that contain a specific keyword in their name and append each file's contents to a 2D array.

    Parameters
    ----------
    folder_path : str
        Path to the folder containing the files.
    keyword : str
        String that must be part of the file's name to be processed.

    Returns
    -------
    numpy.ndarray
        A 2D array containing the contents of each processed file, with the data type of the array elements as float.
    """
    # Initialize the 2D array
    data_2d_array = []
    
    # List all files in the given folder
    for file_name in os.listdir(folder_path):
        # Check if the file name contains the keyword
        if keyword in file_name.split('_'):
            # Construct full file path
            file_path = os.path.join(folder_path, file_name)
            # Open and read the file
            with open(file_path, 'r') as file:
                # Assuming each line of a file represents a separate data entry
                file_contents = [line.strip() for line in file.readlines()]
                data_2d_array.append(file_contents)
    
    return np.array(data_2d_array).astype(float)


def find_shared_p_values(pattern: str, pathdir: str, extension: str = '.pkl') -> List[Tuple[float, int]]:
    """
    Identifies and counts p-values found in file paths that appear in at least two different subdirectories,
    considering files with a specific extension.
    
    This function searches for files with the given extension within subdirectories of a specified path, 
    extracts p-values based on a provided pattern, and counts the occurrences across different subdirectories. 
    It returns a sorted list of unique p-values that are present in two or more subdirectories, along with 
    their occurrence count.

    Parameters
    ----------
    pattern : str
        The regular expression pattern used to extract p-values from file paths. 
        The pattern should contain a capturing group for the p-value.
    pathdir : str
        The root directory path under which to search for files with the specified extension. 
        This path is expected to contain subdirectories where the files are located.
    extension : str, optional
        The file extension to look for (default is '.pkl').

    Returns
    -------
    list of tuple
        A sorted list of tuples, where each tuple contains a p-value and the count of subdirectories 
        in which the p-value is found. The list is sorted by p-value in ascending order.
        
    Examples
    --------
    >>> find_shared_p_values(r"p=([\\d.]+)", "/path/to/data/")
    [(1.0, 2), (2.5, 3)]
    
    This example indicates that p-value 1.0 was found in 2 subdirectories, and p-value 2.5 was found in 3 subdirectories.
    """
    
    # Dictionary to hold sets of subdirectories for each found p value
    p_values_dirs = defaultdict(set)

    # Use glob to iterate over all files with the specified extension in subfolders
    for filepath in glob.glob(f'{pathdir}*/*{extension}'):
        match = re.search(pattern, filepath)
        if match:
            # Extract the p value
            p_value = float(match.group(1))
            # Extract subdirectory from the filepath
            # Adjust the split index based on your path structure
            subdirectory = filepath.split('/')[4]
            # Add the subdirectory to the set for this p value
            p_values_dirs[p_value].add(subdirectory)

    # Prepare a list to hold p values and the number of sharing subdirectories
    p_values_shared_count = []

    # Filter and count p values that appear in at least two different subdirectories
    for p_value, dirs in p_values_dirs.items():
        num_shared = len(dirs)
        if num_shared >= 2:
            p_values_shared_count.append((p_value, num_shared))

    # Sort the list by p value
    p_values_shared_count.sort()
    
    return p_values_shared_count


def remove_empty_dirs(path: str):
    """
    Remove all empty subdirectories within a given directory path.

    This function walks through the directory tree, starting from the bottom. It attempts to remove each directory
    if it is empty. If a directory is not empty or another error occurs (e.g., insufficient permissions), the error
    is caught and printed, but does not interrupt the process.

    Parameters:
    - path (str): The path to the directory from which empty subdirectories should be removed.

    Returns:
    None. Outputs to the console the path of each directory that is removed, or any errors encountered.
    """
    import os

    for root, dirs, files in os.walk(path, topdown=False):
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            try:
                os.rmdir(dir_path)
                print(f"Removed empty directory: {dir_path}")
            except OSError as e:
                print(e)


def create_symmetric_log_bins(min_val, max_val, num_bins, incMagn=2):
    """Creates symmetric logarithmic bins and their centers."""
    bins_positive = np.logspace(np.log10(min_val)-incMagn, np.log10(max_val)+incMagn, num_bins//2 + 1)
    bins_negative = -np.flip(bins_positive[:-1])
    bins = np.concatenate((bins_negative, [0], bins_positive))
    bin_centers = (bins[:-1] + bins[1:]) / 2
    return bins, bin_centers
def bin_eigenvalues(eig_values, bins, bin_centers):
    """Bins eigenvalues and counts occurrences in each bin, using bin centers as keys."""
    indices = np.digitize(eig_values, bins, right=True) - 1
    indices = np.clip(indices, 0, len(bin_centers) - 1)  # Ensure indices are within the valid range
    bin_keys = [bin_centers[index] for index in indices]
    return Counter(bin_keys)



class UnionFind:
    """
    A Union-Find or Disjoint Set Union (DSU) data structure with path compression and union by rank.
    
    It provides an efficient way to manage a partition of a set into disjoint subsets and is useful 
    for dealing with connectivity queries, particularly in graph algorithms.

    Attributes:
    -----------
    parent (List[int]): Stores the parent of each element. Initially, each element is its own parent.
    rank (List[int]): Represents the rank of each element, used to keep the tree flat by attaching 
                       the root of the smaller tree under the root of the larger tree.
    """

    def __init__(self, n):
        """
        Initializes the UnionFind structure with n elements.

        Parameters:
        -----------
        n (int): The number of elements in the Union-Find structure.
        
        Returns:
        --------
        None
        """
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, p):
        """
        Finds the representative of the set containing 'p' with path compression.
        Path compression flattens the structure of the tree whenever `find` is used,
        so that future operations also benefit from the flatter tree.

        Parameters:
        -----------
        p (int): The element whose set representative is to be found.
        
        Returns:
        --------
        int: The representative of the set containing 'p'.
        """
        if self.parent[p] != p:
            self.parent[p] = self.find(self.parent[p])
        return self.parent[p]

    def union(self, p, q):
        """
        Merges the sets containing 'p' and 'q'. It uses the union by rank strategy,
        which attaches the tree with less depth (smaller rank) under the root of the deeper tree
        (larger rank) to keep the tree flat.

        Parameters:
        -----------
        p (int): An element of the first set.
        q (int): An element of the second set.
        
        Returns:
        --------
        None
        """
        rootP = self.find(p)
        rootQ = self.find(q)
        if rootP != rootQ:
            if self.rank[rootP] > self.rank[rootQ]:
                self.parent[rootQ] = rootP
            elif self.rank[rootP] < self.rank[rootQ]:
                self.parent[rootP] = rootQ
            else:
                self.parent[rootQ] = rootP
                self.rank[rootP] += 1


def find_largest_cluster_circle2D(circles, radius):
    """
    Identifies the largest cluster of overlapping circles given their centers and a common radius.

    Parameters:
    -----------
    circles (np.array): A numpy array of tuples where each tuple represents the (x, y) coordinates of a circle's center.
    radius (float): The radius of each circle.

    Returns:
    --------
    list: A list of tuples representing the centers of the circles in the largest cluster.

    Description:
    ------------
    The function utilizes a KDTree for efficient spatial queries to detect overlapping circles.
    It employs a union-find data structure to group and identify clusters of overlapping circles.
    The function returns the largest cluster found.
    """
    from scipy.spatial import KDTree
    tree = KDTree(circles)
    uf = UnionFind(len(circles))
    threshold = 2 * radius

    for i in range(len(circles)):
        neighbors = tree.query_ball_point(circles[i], r=threshold)
        for j in neighbors:
            if i != j:
                uf.union(i, j)

    clusters = defaultdict(list)
    for i in range(len(circles)):
        root = uf.find(i)
        clusters[root].append(tuple(circles[i]))

    largest_cluster = max(clusters.values(), key=len)
    return largest_cluster

def compose(f: callable, g: callable, g_args: tuple = (), g_kwargs: dict = None) -> callable:
    """
    Compose two functions f and g to create a new function that represents g(f(x)), 
    optionally passing additional arguments to g.

    Parameters:
    - f (callable): The first function to apply.
    - g (callable): The second function to apply to the result of f.
    - g_args (tuple, optional): Positional arguments to pass to function g.
    - g_kwargs (dict, optional): Keyword arguments to pass to function g.

    Returns:
    - callable: A new function that takes an input x and returns g(f(x), *g_args, **g_kwargs).

    Example:
    >>> def f1(x):
    ...     return x + 1
    >>> def f2(x, factor=2):
    ...     return x * factor
    >>> composed_function = compose(f1, f2, g_args=(), g_kwargs={'factor': 3})
    >>> composed_function(3)
    12  # f2(f1(3), factor=3) = f2(4, factor=3) = 12
    """
    g_kwargs = g_kwargs or {}

    def composed_function(*args, **kwargs):
        f_result = f(*args, **kwargs)  # Execute the first function with all arguments
        return g(f_result, *g_args, **g_kwargs)  # Pass the result of f to g, with any additional args/kwargs

    return composed_function
def regbin_ndarr(eigV: NDArray) -> NDArray:
    """
    Regularizes and binarizes a NumPy array by setting all zero elements to +1 and taking the sign of each element.

    For each element in the input array `eigV`:
    - If the element is zero, it is replaced with +1.
    - Otherwise, the sign of the element is taken (resulting in -1 for negative values and +1 for positive values).

    Parameters
    ----------
    eigV : NDArray
        A NumPy array of numerical values to be regularized and binarized.

    Returns
    -------
    NDArray
        A NumPy array of the same shape as `eigV`, where:
        - All zero elements are set to +1.
        - All positive elements are set to +1.
        - All negative elements are set to -1.

    Example
    -------
    >>> import numpy as np
    >>> eigV = np.array([-2.5, 0.0, 3.1, -0.7, 0.0])
    >>> regbin_ndarr(eigV)
    array([-1.,  1.,  1., -1.,  1.])

    Notes
    -----
    - This function is useful for converting continuous data into binary form, especially in contexts where the presence or absence (or positive/negative sign) of a feature is significant.
    - The function handles multi-dimensional arrays as well.

    See Also
    --------
    numpy.sign : Returns an element-wise indication of the sign of a number.
    numpy.where : Return elements chosen from `eigV` or +1 depending on the condition.

    """
    return np.sign(np.where(eigV == 0, +1, eigV))

def project_3d_to_2d(x, y, z, theta=0., phi=0.):
    """
    Projects a 3D point (x, y, z) onto a 2D plane using specified rotation angles.

    Parameters
    ----------
    x : float
        The x-coordinate of the 3D point.
    y : float
        The y-coordinate of the 3D point.
    z : float
        The z-coordinate of the 3D point.
    theta : float, optional
        Rotation angle around the y-axis, in radians. If not provided, uses self.theta.
    phi : float, optional
        Rotation angle around the x-axis, in radians. If not provided, uses self.phi.

    Returns
    -------
    tuple of float
        The (x, y) coordinates of the point projected onto the 2D plane.
    """
    # Rotation matrix around the y-axis (theta)
    R_theta = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

    # Rotation matrix around the x-axis (phi)
    R_phi = np.array([
        [1, 0, 0],
        [0, np.cos(phi), -np.sin(phi)],
        [0, np.sin(phi), np.cos(phi)]
    ])

    # Initial position vector
    position = np.array([x, y, z])

    # Apply rotations (order matters)
    position_rotated = R_phi @ R_theta @ position

    # Project onto 2D plane (ignore z after rotation)
    x2, y2 = position_rotated[0], position_rotated[1]

    return x2, y2

def interpolate_grid_data(
    x: np.ndarray, 
    y: np.ndarray, 
    z: np.ndarray, 
    num_points: int = 1000
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Performs interpolation of z-data over a uniformly spaced grid defined by x and y coordinates.

    Parameters
    ----------
    x : np.ndarray
        Meshgrid of x-coordinate values.
    y : np.ndarray
        Meshgrid of y-coordinate values.
    z : np.ndarray
        Matrix of computed z-values corresponding to each (x, y) pair.
    num_points : int, optional
        Number of points along each axis for the new grid. Defaults to 1000.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        grid_x : np.ndarray
            Meshgrid of interpolated x-coordinate values.
        grid_y : np.ndarray
            Meshgrid of interpolated y-coordinate values.
        z_new : np.ndarray
            Interpolated z-values on the new grid.
    """
    points = np.column_stack((x.ravel(), y.ravel()))
    grid_x, grid_y = np.meshgrid(
        np.linspace(x.min(), x.max(), num_points),
        np.linspace(y.min(), y.max(), num_points)
    )
    z_new = griddata(
        points, z.ravel(), (grid_x, grid_y), method='cubic', fill_value=np.nan
    )
    return grid_x, grid_y, z_new
from typing import List, Any

def uniques(lst: List[Any]) -> List[Any]:
    """
    Returns a list of unique elements from the input list. This function leverages Python's 
    built-in `set` data structure to eliminate duplicate entries efficiently. Note that the 
    original order of elements is not preserved.

    Parameters
    ----------
    lst : List[Any]
        The input list from which to extract unique elements. The list can contain elements of 
        any data type that is hashable.

    Returns
    -------
    List[Any]
        A new list containing only the unique elements from the input list, with duplicates removed.

    Examples
    --------
    ```python
    >>> uniques([1, 2, 2, 3, 4, 4, 5])
    [1, 2, 3, 4, 5]

    >>> uniques(['apple', 'banana', 'apple', 'cherry'])
    ['apple', 'banana', 'cherry']
    ```

    Notes
    -----
    - **Order Preservation**: This function does **not** preserve the original order of elements. 
      If maintaining order is essential, consider using alternative methods such as `dict.fromkeys` 
      or `collections.OrderedDict`.
    - **Hashable Elements**: All elements in the input list must be hashable. Unhashable elements 
      (e.g., lists, dictionaries) will raise a `TypeError`.

    See Also
    --------
    dict.fromkeys : Create a dictionary with keys from the input list, preserving order (Python 3.7+).
    collections.OrderedDict : Ordered dictionary for maintaining element order (pre-Python 3.7).
    """
    return list(set(lst))
