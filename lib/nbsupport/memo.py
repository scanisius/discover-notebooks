import numpy
import pandas

import switching


def memo_test(events, selected_genes, groups, permutations=10000):
    groups_memo = [pandas.match(group, selected_genes) for group in groups]
    events_selected = events[selected_genes]
    sampler = switching.EventMatrixSampler(events_selected.astype(int), "gobbi")
    
    coverages = numpy.array([events_selected[i].any(0).sum() for i in groups_memo])
    higher_coverage = numpy.zeros_like(coverages)
    
    for i in xrange(permutations):
        null_sample = sampler.sample()
        higher_coverage += (numpy.array([null_sample[i].any(0).sum() for i in groups_memo]) >= coverages).astype(int)
    
    return (higher_coverage + 1.0) / (permutations + 1.0)

