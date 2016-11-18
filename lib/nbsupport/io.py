import discover
import pandas


def save_discover_matrix(path_or_buf, key, matrix):
    matrix._events.to_hdf(path_or_buf, key + "/events", complevel=9, complib="bzip2")
    matrix._bg.to_hdf(path_or_buf, key + "/bg", complevel=9, complib="bzip2")


def load_discover_matrix(path_or_buf, key):
    return discover.DiscoverMatrix(
        pandas.read_hdf(path_or_buf, key + "/events"),
        pandas.read_hdf(path_or_buf, key + "/bg"))


def save_pairwise_result(path_or_buf, key, result):
    def f(store):
        store.put(key + "/pvalues", result.pvalues)
        store.put(key + "/qvalues", result.qvalues)
        attrs = store.get_node(key)._v_attrs
        attrs["pi0"] = result.pi0
        attrs["alternative"] = result.alternative

    if isinstance(path_or_buf, pandas.HDFStore):
        f(path_or_buf)
    else:
        with pandas.HDFStore(path_or_buf, complevel=9, complib="bzip2") as store:
            f(store)


def load_pairwise_result(path_or_buf, key):
    def f(store):
        attrs = store.get_node(key)._v_attrs
        return discover.pairwise.PairwiseDiscoverResult(
            store.get(key + "/pvalues"),
            store.get(key + "/qvalues"),
            attrs["pi0"],
            attrs["alternative"])

    if isinstance(path_or_buf, pandas.HDFStore):
        return f(path_or_buf)
    else:
        with pandas.HDFStore(path_or_buf, "r", complevel=9, complib="bzip2") as store:
            return f(store)
