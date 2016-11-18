import matplotlib
import matplotlib.pyplot as plt
import numpy
import pandas

import nbsupport.roc

from contextlib import contextmanager
from matplotlib.transforms import Bbox


DEFAULT_MPL_PARAMS = {
    'axes.linewidth': 1.5,
    'font.family': 'arial',
    'font.size': 14,
    'xtick.direction': 'out',
    'ytick.direction': 'out'}


@contextmanager
def use_custom_style(params={}):
    plotParams = DEFAULT_MPL_PARAMS.copy()
    plotParams.update(params)

    with plt.rc_context(plotParams):
        yield
        for ax in plt.gcf().axes:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')
            ax.spines['bottom'].set_position(('outward', 10))
            ax.spines['left'].set_position(('outward', 15))


def event_plot(events, rowHeight=100):
    genes = numpy.unique(numpy.char.partition(events.rownames.astype(str), "_")[:, 0])
    x = pandas.DataFrame(numpy.zeros((len(genes), events.shape[1]), dtype=int), index=numpy.array(genes), columns=events.colnames)

    fig = plt.gcf()
    bbox = Bbox([[0.8 * fig.dpi,
                  0.5 * fig.dpi], 
                 [(fig.get_figwidth() - 0.5) * fig.dpi,
                  0.3 * len(x) * fig.dpi]])
    plt.axes(bbox.transformed(plt.gcf().transFigure.inverted()).bounds)
    
    for feature in events.rownames:
        gene, eventType = feature.split("_", 1)
        i = genes.searchsorted(gene)
        if eventType == "mut":
            shift = 0
        elif eventType == "gain":
            shift = 1
        elif eventType == "loss":
            shift = 2

        x.values[i] |= events[[feature]].events[0].astype(int) << shift


    geneOrder = (x.values > 0).sum(1).argsort()
    sampleOrder = numpy.lexsort(x.values[geneOrder])[::-1]

    colours = matplotlib.colors.ListedColormap([
            "#dddddd",
            "#377eb8",
            "#2ca25f",
            "#756bb1",
            "#de2d26",
            "#ff7f00"])
    
    plt.imshow(
        x.values[geneOrder][:, sampleOrder], aspect="auto", interpolation="none",
        cmap=colours , origin="lower", vmin=0, vmax=5)

    plt.grid(ls="-", c="white", axis="x")
    plt.setp(plt.gca().spines.values(), color="#888888", alpha=0)

    for tic in plt.gca().xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    for tic in plt.gca().yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    try:
        plt.yticks(numpy.arange(len(x)), x.index[geneOrder])
    except AttributeError:
        pass

    plt.gca().yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(
        numpy.linspace(-0.5, len(x) - 0.5, x.shape[0] + 1)))
    plt.grid(ls="-", c="white", axis="y", lw=4, alpha=1, which="minor")
    plt.grid(ls="None", axis="y", which="major")
    for tic in plt.gca().yaxis.get_minor_ticks():
        tic.tick1On = tic.tick2On = False

    plt.gcf().patch.set_facecolor("white")

    plt.xlim(0, x.ix[genes].any(0).sum())


def plot_average_roc_curves(pvalues, labels, method_params):
    with use_custom_style({"font.size": 12,
                           "legend.fontsize": 12,
                           "lines.linewidth": 2}):
        figure = plt.figure(figsize=(4.0, 5.0))

        ax = plt.gca()

        for method, colour, label in method_params:
            coords = [
                nbsupport.roc.pvalue_roc_curve(labels, numpy.array(pvalues[i][method]))[:2]
                for i in xrange(len(pvalues))]

            fpr, tpr = nbsupport.roc.average_roc_curves(coords)
            plt.plot(fpr, tpr, label=label, c=colour)
            print "%s:" % method, nbsupport.roc.auc(fpr, tpr)

        ax.legend(frameon=False, ncol=1, loc="lower right")

        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")


def plot_average_calibration_curves(pvalues, labels, method_params):
    with use_custom_style({"font.size": 12,
                           "legend.fontsize": 12,
                           "lines.linewidth": 2}):
        figure = plt.figure(figsize=(4.0, 5.0))

        ax = plt.gca()

        for method, colour, label in method_params:
            fprs_tprs_thresholds = [
                nbsupport.roc.pvalue_roc_curve(labels, numpy.array(pvalues[i][method]))
                for i in xrange(len(pvalues))]

            threshold, fpr = nbsupport.roc.average_calibration_curves(fprs_tprs_thresholds)
            plt.plot(threshold, fpr, label=label, c=colour)

        ax.legend(frameon=False, ncol=1, loc="lower right")

        ax.set_xlabel("Significance level $\\alpha$")
        ax.set_ylabel("False positive rate")


def plot_average_sensitivity_curves(pvalues, labels, method_params):
    with use_custom_style({"font.size": 12,
                           "legend.fontsize": 12,
                           "lines.linewidth": 2}):
        figure = plt.figure(figsize=(4.0, 5.0))

        ax = plt.gca()

        for method, colour, label in method_params:
            fprs_tprs_thresholds = [
                nbsupport.roc.pvalue_roc_curve(labels, numpy.array(pvalues[i][method]))
                for i in xrange(len(pvalues))]

            threshold, tpr = nbsupport.roc.average_sensitivity_curves(fprs_tprs_thresholds)
            plt.plot(threshold, tpr, label=label, c=colour)

        ax.legend(frameon=False, ncol=1, loc="lower right")

        ax.set_xlabel("Significance level $\\alpha$")
        ax.set_ylabel("True positive rate")
        plt.gca().set_xscale("log")
