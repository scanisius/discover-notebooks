from __future__ import absolute_import

import io
import os
import urllib
import zipfile

import scipy.stats

import nbsupport.util

import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

from rpy2 import robjects


MEGSA_URL = "http://dceg.cancer.gov/tools/analysis/megsa/MEGSA_beta.zip"


def install(download_dir):
    local_filename = "%s/MEGSA_beta.zip" % download_dir
    if not os.path.exists(local_filename):
        filename, response = urllib.urlretrieve(MEGSA_URL, local_filename)
    else:
        filename = local_filename
    nbsupport.util.check_digest(filename, "655e1ec48a67530672303d8ac7fbc925")

    with zipfile.ZipFile(local_filename) as archive:
        with io.TextIOWrapper(archive.open("version beta/MEGSA.R")) as stream:
            robjects.reval(stream.read())

    global megsa
    def megsa(events):
        s = robjects.r.funEstimate(events.T).rx2("S")[0]
        return 0.5 * scipy.stats.chisqprob(s, 1) + 0.5 * int(s == 0)
