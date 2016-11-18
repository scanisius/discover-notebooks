from collections import namedtuple

MuexResult = namedtuple("MuexResult", "gamma delta alpha beta pvalue statistic")


try:
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    
    from rpy2.robjects.packages import importr

    muex_package = importr("muex")

    __version__ = muex_package.___NAMESPACE___["spec"].rx2("version")[0]

    def muex(*args):
        result = muex_package.muex(*args)
        return MuexResult(**dict(result.items()))

except:
    import warnings
    warnings.warn("Unable to import the muex R package. "
                  "Using a Python-based reimplementation instead. "
                  "This should give the same results, but note "
                  "that the results for the DISCOVER paper were "
                  "generated using the R package.")
    
    import numpy
    import scipy.stats


    __version__ = "1.0 (Python-based reimplementation)"


    def muex(Y, err=False, theta_given=None, i_no=20, iter=5000, verbose=False):
        if err:
            raise RuntimeError("muex with non-zero error rates is not supported by this Python version")
        
        theta = muexEstNoErr(Y)
        if theta_given is not None:
            if verbose:
                print "Note: The values of alpha and beta are assumed to be 0, and their values provided in theta_given will be ignored."
            if "gamma" in theta_given:
                theta["gamma"] = theta_given["gamma"]
            if "delta" in theta_given:
                theta["delta"] = theta_given["delta"]

        theta["alpha"] = 0
        theta["beta"] = 0
        test_result = muexTestNoErr(Y, theta)

        return MuexResult(
            gamma=theta["gamma"], delta=theta["delta"], alpha=theta["alpha"],
            beta=theta["beta"], pvalue=test_result["pvalue"], statistic=test_result["statistic"])


    def muexEstNoErr(Y):
        ks = Y.sum(1)
        m, n = Y.shape
        gamma = (ks > 0).mean()
        delta = max(0, 1.0 * Y.sum() / ((n - 1) * m * gamma) - 1.0 / (n - 1))
        return dict(gamma=gamma, delta=delta)


    def muexTestNoErr(Y, theta):
        n = Y.shape[1]
        pi_est = F1_i_parest(Y)
        lik_ind = F1_i_llik(Y, pi_est)
        lo1 = numpy.apply_along_axis(muexLlikNoErr, 1, Y, theta)
        lik_g1 = lo1.sum()
        lo2 = F1_i_llik_o(Y, pi_est)
        vuong = F1_vuong(lo1, lo2, coef_no1=2, coef_no2=n)
        return dict(pvalue=vuong["pvalue"], statistic=vuong["statistic"])


    def F1_i_parest(Y):
        return Y.mean(0)


    def F1_i_llik(Y, pi):
        kg = Y.sum(0)
        m = Y.shape[0]
        return numpy.sum(numpy.log(pi**kg) + numpy.log((1 - pi)**(m - kg)))


    def muexLlikNoErr(y, theta):
        k = y.sum()
        n = len(y)
        if theta["gamma"] == 0 or theta["delta"] == 0:
            if k > 0:
                p = theta["gamma"] * 1.0 * k / n * theta["delta"]**(k - 1) * (1 - theta["delta"])**(n - k)
            else:
                p = 1 - theta["gamma"]
            l = numpy.log(p)
        else:
            if k > 0:
                l = numpy.log(theta["gamma"]) + numpy.log(1.0 * k / n) + (k - 1) * numpy.log(theta["delta"]) + (n - k) * numpy.log(1 - theta["delta"])
            else:
                l = numpy.log(1 - theta["gamma"])
        return l


    def F1_i_llik_o(Y, pi):
        return numpy.apply_along_axis(F1_i_llik_y, 1, Y, pi)


    def F1_i_llik_y(y, pi):
        lpi = numpy.log(pi**y)
        lpi1 = numpy.log((1 - pi)**(1 - y))
        return numpy.sum(lpi + lpi1)


    def F1_vuong(lo1, lo2, coef_no1, coef_no2, verbose=False):
        lr_o = lo1 - lo2
        lr_full = lo1.sum() - lo2.sum()
        s = numpy.sqrt(numpy.mean(lr_o**2) - (numpy.mean(lr_o))**2)
        m = len(lr_o)
        correct = (numpy.log(m)/2) * (coef_no1 - coef_no2)
        v = (lr_full - correct)/(numpy.sqrt(m) * s)
        if verbose:
            print "Vuong Non-Nested Hypothesis Test-Statistic:", v
            if v > 0:
                print "model1 > model2, with p-value", 1 - scipy.stats.norm.cdf(v)
            else:
                print "model2 > model1, with p-value", scipy.stats.norm.cdf(v)

        return dict(pvalue=1 - scipy.stats.norm.cdf(v), statistic=v)
