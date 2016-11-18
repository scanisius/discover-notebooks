import hashlib
import numpy
import random


NUMPY_RANDOM_SEED = 1234
PYTHON_RANDOM_SEED = 1234


def set_random_seed():
    random.seed(PYTHON_RANDOM_SEED)
    numpy.random.seed(NUMPY_RANDOM_SEED)


def align_columns(*dataframes):
    def alignment_helper(args):
        orders = map(numpy.argsort, args)
        indexes = numpy.zeros(len(args), dtype=int)

        while all(i < len(x) for i, x in zip(indexes, args)):
            if all(args[i - 1][orders[i - 1][indexes[i - 1]]] == args[i][orders[i][indexes[i]]]
                   for i in xrange(1, len(args))):
                yield [orders[i][indexes[i]] for i in xrange(len(args))]
                indexes += 1
            else:
                values = numpy.array([args[i][orders[i][indexes[i]]]
                                      for i in xrange(len(args))])
                indexes[values < max(values)] += 1

    column_labels = [df.columns for df in dataframes]
    return tuple(dataframes[i].iloc[:, selection]
                 for i, selection
                 in enumerate(zip(*alignment_helper(column_labels))))


def check_digest(filename, md5sum, buffer_size = 1024**2):
    md5 = hashlib.md5()

    with open(filename) as stream:
        buf = stream.read(buffer_size)
        while buf != "":
            md5.update(buf)
            buf = stream.read(buffer_size)
       
    if md5.hexdigest() != md5sum:
        raise RuntimeError("%s differs from the file used for the original analysis" % filename)

