import numpy as np
from . import ArrayLike


def rpplot_prepare(sa: ArrayLike, counts: ArrayLike):
    """
    Prepare data for range pair plot by sorting the input arrays.

    Parameters
    ----------
    sa : ArrayLike
        The stress amplitude array.
    counts : ArrayLike
        The counts array corresponding to the stress amplitudes.

    Returns
    -------
    dict
        A dictionary with keys "sa" and "counts", containing the sorted and formatted data,
        and "sa_max" with the maximum value that "sa" contains.
    """
    sa = np.asarray(sa).flatten()
    counts = np.asarray(counts).flatten()

    # Filter out elements with non-zero counts
    mask = counts > 0
    sa, counts = sa[mask], counts[mask]
    sa_max = np.nan
    if np.any(mask):
        # Sort descending by stress amplitude
        i = np.argsort(-sa)
        sa, counts = sa[i], counts[i]

        # First stress amplitude is the maximum
        sa_max = sa[0]

    return {"sa": sa, "counts": counts, "sa_max": sa_max}
