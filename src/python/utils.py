import numpy as np
from . import ArrayLike


def rpplot_prepare(sa: ArrayLike, counts: ArrayLike):
    """
    Prepare data for rainflow plot by sorting and formatting the input arrays.

    Parameters
    ----------
    sa : ArrayLike
        The stress amplitude array.
    counts : ArrayLike
        The counts array corresponding to the stress amplitudes.

    Returns
    -------
    dict
        A dictionary with keys "sa" and "counts", containing the sorted and formatted data.
    """
    # Convert inputs to flattened arrays
    sa = np.asarray(sa).flatten()
    counts = np.asarray(counts).flatten()

    # Stack and sort data in descending order by stress amplitude
    x = np.vstack((sa, counts)).T
    x = x[np.argsort(-x[:, 0])]

    # Add leading and trailing points
    x = np.vstack(((x[0, 0], 0), x, (0, 0.01)))

    return {"sa": x[:, 0], "counts": x[:, 1]}
