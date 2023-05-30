import numpy as np


def rpplot_prepare(sa, counts):
    x = np.vstack((np.array(sa).flatten(), np.array(counts).flatten())).T
    i = np.argsort(x[:, 0])
    x = x[i[::-1], :]
    x = np.vstack(((x[0, 0], 0), x, (0, 0.01)))

    return {"sa": x[:, 0], "counts": x[:, 1]}
