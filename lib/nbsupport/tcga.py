import numpy
import pandas


PANCAN12_STUDIES = ['COAD', 'LUSC', 'READ', 'GBM', 'LAML', 'HNSC',
                    'BLCA', 'UCEC', 'LUAD', 'OV','BRCA', 'KIRC']


def read_copynumber_data(stream):
    header = stream.readline().rstrip().split("\t")
    assert header[:3] == ["Gene Symbol", "Locus ID", "Cytoband"]
    assert header[3].startswith("TCGA-")
    stream.seek(0)
    return pandas.read_table(stream, index_col=0, usecols=header[:1] + header[3:])


def read_mutation_data(stream):
    maf = pandas.read_table(stream, index_col=False, usecols=["Tumor_Sample_Barcode", "Hugo_Symbol"])
    return numpy.sign(maf.pivot_table(
        index=["Hugo_Symbol"],
        columns=["Tumor_Sample_Barcode"],
        aggfunc=len,
        fill_value=0))


def read_gistic_output(filename):
    with open(filename) as stream:
        header = stream.readline().rstrip("\n").split("\t")
        assert header[0] == "cytoband"

        segments = { name: [] for name in header[1:] if name != "" }

        for i in xrange(3):
            stream.readline()

        line = stream.readline().rstrip("\n").split("\t")
        assert line[0] == "genes in wide peak"
        while len(line) == len(header):
            for i, gene in enumerate(line[1:], 1):
                if gene != "":
                    segments[header[i]].append(gene)

            line = stream.readline().rstrip("\n").split("\t")

    return segments
