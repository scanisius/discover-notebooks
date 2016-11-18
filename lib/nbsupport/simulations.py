import networkx
import numpy
import scipy.stats

from itertools import combinations
from statsmodels.regression import quantile_regression


####
#
# Simulate independent alteration data
#

def generate_independent_alterations(geneMarginals, sampleMarginals):
    graph = networkx.bipartite.configuration_model(
        geneMarginals, sampleMarginals, create_using=networkx.Graph())
    return networkx.algorithms.bipartite.biadjacency_matrix(
        graph, [n for n in graph if graph.node[n]["bipartite"] == 0]).toarray()


####
#
# Simulate mutual exclusivities
#

def generate_mutex(events, genes, qreg):
    # 1. Pick a random gene
    i = numpy.random.choice(genes)
    
    # 2. Sample the coverage of the second gene
    p = numpy.random.choice(events[genes].sum(1))
   
    # 3. Generate the second gene vector
    x = numpy.zeros_like(events[i])
    subset = events[i] == 0
    weights = events.sum(0)[subset]
    k = numpy.random.choice(subset.nonzero()[0], p, False, 1.0 * weights / weights.sum())
    x[k] = 1
        
    cov = 1.0 * (events[i].sum() + x.sum()) / events.shape[1]
    
    # 4. Sample the impurity parameter
    impurity = numpy.random.random() * qreg.predict(cov)
    
    x[events[i] == 1] = scipy.stats.bernoulli(impurity / (1.0 * events[i].sum() / (events[i].sum() + x.sum()))).rvs(numpy.sum(events[i] == 1))
    
    return i, x


def add_mutual_exclusivities(events, numMutex=500):
    events = numpy.asarray(events)
    
    genesMinFreq50 = numpy.where(events.sum(1) >= 50)[0]
    from discover.stats import mutexParams

    mutexStats = numpy.rec.fromrecords(
        [mutexParams(events[[i, j]])
         for i, j
         in combinations(genesMinFreq50, 2)],
        names=["coverage", "impurity", "balance"])

    mutexQReg = quantile_regression.QuantReg(
        mutexStats["impurity"], mutexStats["coverage"]).fit(0.01)

    mutexGeneIndexes = []
    mutexGeneEvents = []

    for i in xrange(numMutex):
        j, e = generate_mutex(events, genesMinFreq50, mutexQReg)
        mutexGeneIndexes.append((j, events.shape[0] + i))
        mutexGeneEvents.append(e)

    eventsWithMutex = numpy.row_stack([events, mutexGeneEvents])

    return eventsWithMutex, numpy.array(mutexGeneIndexes)


def generate_mutex_group(filtered_gene_counts, sample_marginals, coverage=0.5, impurity=0.05):
    num_samples = len(sample_marginals)
    num_covered_samples = int(coverage * num_samples)

    gene_coverages = []
    while sum(gene_coverages) < num_covered_samples:
        gene_coverages.append(numpy.random.choice(filtered_gene_counts))

    p = 1.0 * sample_marginals / sample_marginals.sum()
    altered_samples = numpy.random.choice(num_samples, num_covered_samples, False, p)

    num_genes = len(gene_coverages)
    events = numpy.zeros((num_genes, num_samples))

    gene_weights = numpy.array(gene_coverages, dtype=float)
    for sample in altered_samples:
        gene = numpy.random.choice(num_genes, p=gene_weights/gene_weights.sum())
        events[gene, sample] = 1

    current_impurity = 1.0 * (events.sum(0) > 1).sum() / events.any(0).sum()
    while current_impurity < impurity:
        i = numpy.random.randint(events.shape[0])
        j = numpy.random.choice(altered_samples, p=p[altered_samples] / p[altered_samples].sum())
        events[i, j] = 1
        current_impurity = 1.0 * (events.sum(0) > 1).sum() / events.any(0).sum()

    return events


def add_mutex_groups(events, numSets=100, impurityParams=[0.02, 0.05, 0.08], minGenes=3, maxGenes=6):
    sampleMarginals = events.sum(0)
    geneMarginals = events.sum(1)
    filteredGeneMarginals = geneMarginals[geneMarginals > 50]
    
    # generate positive mutex groups
    impurities = numpy.random.choice(impurityParams, numSets)

    posGroups = []
    for impurity in impurities:
        numGenes = 0
        while not (minGenes <= numGenes <= maxGenes):
            x = generate_mutex_group(
                    filteredGeneMarginals, sampleMarginals,
                    scipy.stats.truncnorm((0.2 - 0.4) / 0.2, (0.8 - 0.4) / 0.2, 0.4, 0.2).rvs(),
                    impurity)
            numGenes = x.shape[0]
        posGroups.append(x)
 
    indices = events.shape[0] + numpy.array([[0, len(x)] for x in posGroups]).ravel().cumsum().reshape((len(posGroups), 2))
    posGenes = [numpy.arange(start, end) for start, end in indices]
    
    events = numpy.row_stack([events] + posGroups)

    # generate negative mutex groups
    dist = scipy.stats.lognorm(*scipy.stats.lognorm.fit(events.sum(1)[numpy.unique(numpy.concatenate(posGenes))]))
    gene_weights = dist.pdf(events.sum(1))
    gene_weights[events.sum(1) < 50] = 0
    gene_weights[numpy.concatenate(posGenes)] = 0
    gene_weights /= gene_weights.sum()
    
    negGenes = [numpy.random.choice(gene_weights.shape[0], len(g), replace=False, p=gene_weights)
                for g in posGenes]
    
    return events, posGenes, negGenes


####
#
# Simulate co-occurrences
#

def generate_cooc(events, genes, qreg):
    i = numpy.random.choice(genes)
    
    k = scipy.stats.geom(0.2).rvs()
    coverage = events[i].sum() + k
    
    minOverlap = qreg.predict(coverage)[0]
    overlap = minOverlap + scipy.stats.beta(5, 1).rvs() * (events[i].sum() - minOverlap)
    
    x = numpy.zeros_like(events[i])
    x[numpy.random.choice(events[i].nonzero()[0], overlap, False)] = 1
    x[numpy.random.choice((events[i] == 0).nonzero()[0], coverage - events[i].sum(), False)] = 1
    
    return i, x


def add_cooccurrences(events, numCooc=500):
    events = numpy.asarray(events)

    genesMinFreq50 = numpy.where(events.sum(1) >= 50)[0]

    coocStats = numpy.rec.fromrecords(
        [[events[[i, j]].all(0).sum(), events[[i, j]].any(0).sum()]
         for i, j in combinations(genesMinFreq50, 2)],
        names=["overlap", "coverage"])

    coocQReg = quantile_regression.QuantReg(
        coocStats["overlap"], coocStats["coverage"]).fit(0.99)

    coocGeneIndexes = []
    coocGeneEvents = []

    for i in xrange(numCooc):
        j, e = generate_cooc(events, genesMinFreq50, coocQReg)
        coocGeneIndexes.append((j, events.shape[0] + i))
        coocGeneEvents.append(e)

    eventsWithCooc = numpy.row_stack([events, coocGeneEvents])

    return eventsWithCooc, numpy.array(coocGeneIndexes)
