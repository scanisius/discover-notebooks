import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

from rpy2 import robjects
from rpy2.robjects.packages import importr

timex_package = importr("TiMEx")
__version__ = timex_package.___NAMESPACE___["spec"].rx2("version")[0]

def timex_test(events, groups):
    events_r = robjects.Matrix(events.T)
    return [
        timex_package.testCliqueAsGroup(robjects.IntVector(group), events_r).rx2("pvalueLRT")[0]
        for group in groups]
