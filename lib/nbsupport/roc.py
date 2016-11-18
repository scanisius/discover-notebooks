import numpy

from scipy.interpolate import interp1d


def pvalue_roc_curve(labels, pvalues):
    sorted_pvalues = numpy.sort(pvalues)
    sorted_labels = labels[pvalues.argsort()]

    unique_thresholds = numpy.unique(numpy.r_[0, numpy.logspace(-70, -1, 1000),
                                              numpy.linspace(0.1, 1, 1000), 1.1])

    indices = sorted_pvalues.searchsorted(unique_thresholds)
    tp = numpy.append(0, sorted_labels.cumsum())[indices]
    fp = indices - tp

    if fp[0] != 0:
        tp = numpy.r_[0, tp]
        fp = numpy.r_[0, fp]
        unique_thresholds = numpy.r_[unique_thresholds[0] - 1, unique_thresholds]

    p = labels.sum()
    n = len(pvalues) - p

    tpr = 1.0 * tp / p
    fpr = 1.0 * fp / n

    return fpr, tpr, unique_thresholds


def auc(fpr, tpr):
    return numpy.trapz(tpr, fpr)


def average_roc_curves(coords):
    fprs = []
    tprs = []

    for fpr, tpr in coords:
        # get largest tpr for each fpr
        i = fpr.searchsorted(numpy.unique(fpr), "right") - 1
        fprs.append(fpr[i])
        tprs.append(tpr[i])

    funcs = map(interp1d, fprs, tprs)
    x = numpy.unique(numpy.concatenate(fprs))
    y = numpy.mean([f(x) for f in funcs], 0)
    
    if y[0] > 0:
        y = numpy.r_[0, y]
        x = numpy.r_[0, x]

    return x, y


def average_calibration_curves(coords):
    fprs, tprs, ps = map(list, zip(*coords))
    
    for i in xrange(len(fprs)):
        fprs[i] = numpy.r_[fprs[i], 1]
        ps[i] = numpy.r_[ps[i], 1]
        
        if ps[i][0] < 0:
            fprs[i] = fprs[i][1:]
            ps[i] = ps[i][1:]
        
        if ps[i][0] > 0:
            fprs[i] = numpy.r_[0, fprs[i]]
            ps[i] = numpy.r_[0, ps[i]]
        
    
    funcs = map(interp1d, ps, fprs)
    x = numpy.unique(numpy.concatenate(ps))
    x = x[x <= 1]
    y = numpy.mean([f(x) for f in funcs], 0)

    return x, y


def average_sensitivity_curves(coords):
    fprs, tprs, ps = map(list, zip(*coords))
    
    for i in xrange(len(tprs)):
        tprs[i] = numpy.r_[tprs[i], 1]
        ps[i] = numpy.r_[ps[i], 1]
        
        if ps[i][0] < 0:
            tprs[i] = tprs[i][1:]
            ps[i] = ps[i][1:]
        
        if ps[i][0] > 0:
            tprs[i] = numpy.r_[0, tprs[i]]
            ps[i] = numpy.r_[0, ps[i]]
        
    
    funcs = map(interp1d, ps, tprs)
    x = numpy.unique(numpy.concatenate(ps))
    x = x[x <= 1]
    y = numpy.mean([f(x) for f in funcs], 0)

    return x, y
