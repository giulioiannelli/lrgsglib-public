from .shared import *
from .core import Lattice2D, SignedGraph
def create_Lattice2D_WeigvV(cell: str, **kwargs: Any) -> Lattice2D:
    """
    Create and process a Lattice2D object.

    This function initializes a Lattice2D object using the provided arguments,
    flips selected edges, and computes eigenvalues using the additional keyword arguments.

    Parameters
    ----------
    cell : str
        The cell index for selecting edges in the network dictionary.
    **kwargs : Any
        Additional keyword arguments for Lattice2D and compute_k_eigvV methods.
        The function dynamically separates these arguments based on their relevance
        to the Lattice2D constructor or the compute_k_eigvV method.

    Returns
    -------
    Lattice2D
        The processed Lattice2D object with the selected edges flipped and eigenvalues computed.

    Examples
    --------
    Create and process a Lattice2D object:

    >>> lattice = create_Lattice2D_eigvV(
    ...     cell='0',
    ...     side=10,
    ...     geo='square',
    ...     pflip=0.1,
    ...     sgpath='/path/to/subgraph/data',
    ...     pbc=False,
    ...     init_nw_dict=True,
    ...     with_positions=True,
    ...     some_optional_arg=42
    ... )
    """
    from inspect import signature

    # Get the parameters of the Lattice2D and SignedGraph constructors
    l2D_params = signature(Lattice2D).parameters
    sGraph_params = signature(SignedGraph).parameters

    # Combine parameter keys from both Lattice2D and SignedGraph
    l2D_keys = set(l2D_params.keys()).union(sGraph_params.keys())

    # Separate kwargs for constructors and compute_k_eigvV
    kwargs_l = {k: v for k, v in kwargs.items() if k in l2D_keys}
    kwargs_c = {k: v for k, v in kwargs.items() if k not in l2D_keys}

    # Create the Lattice2D object with additional kwargs
    l = Lattice2D(**kwargs_l)

    # Process the Lattice2D object with additional compute_k_eigvV kwargs
    l.flip_sel_edges(l.nwDict[cell]['G'])
    l.compute_k_eigvV(**kwargs_c)

    return l