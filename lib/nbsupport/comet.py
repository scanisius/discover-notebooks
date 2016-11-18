import numpy
import comet

from cComet_with_timeout import exact_test, binom_test, precompute_factorials


def comet_test(events, timeout=60):
    table = numpy.histogramdd(
        events.T, [(0, 0.5, 1)] * events.shape[0])[0].ravel().astype(int)
    precompute_factorials(events.shape[1])
    num_tables, p = exact_test(events.shape[0], events.shape[1],
                               table.tolist(), 1.1, timeout)
    if p < 0:
        # If the CoMEt test timed out, return the binomial approximation
        p = binom_test(events.shape[0], events.shape[1], table.tolist(), 1.1)
    return p
