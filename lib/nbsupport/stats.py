import numpy


def cummin(x):
    result = numpy.empty(len(x))
    result[0] = x[0]
    for i in xrange(1, len(x)):
        result[i] = min(result[i - 1], x[i])

    return result


def fdr(p, pi0=1.0):
    if not 0 <= pi0 <= 1:
        raise ValueError("Invalid value for pi0: %s. Legal values are between 0 and 1" % pi0)

    nna = ~numpy.isnan(p)
    q = numpy.repeat(numpy.nan, len(p))
    p = p[nna]

    i = numpy.arange(len(p), 0, -1)
    o = numpy.argsort(p)[::-1]
    ro = numpy.argsort(o)

    q[nna] = numpy.minimum(1, cummin(float(pi0) * len(p) / i * p[o])[ro])
    return q

