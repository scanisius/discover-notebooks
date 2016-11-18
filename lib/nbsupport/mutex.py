import numpy
import scipy.stats


def mutex(events):
    col_marginals = events.sum(0)
    return max(
        scipy.stats.fisher_exact(numpy.histogram2d(
            events[i],
            (col_marginals - events[i] > 0).astype(int), [0, 0.5, 1.0])[0], "less")[1]
        for i in xrange(events.shape[0]))
