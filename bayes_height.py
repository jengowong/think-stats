"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: bayes_height.py
"""

import matplotlib.pyplot as pyplot
import math
import numpy
import cPickle
import random
import correlation
import _03_thinkstats
import _04_Pmf
import _05_myplot
import _13_Cdf
import _17_rankit
import _19_brfss


def _make_uniform_prior(t, num_points, label, spread=3.0):
    """
    Makes a prior distribution for mu and sigma based on a sample.

    Args:
        t:          sample
        num_points: number of values in each dimension
        label:      string label for the new Pmf
        spread:     number of standard errors to include

    Returns:
        Pmf that maps from (mu, sigma) to prob.
    """
    # estimate mean and stddev of t
    n = len(t)
    xbar, S2 = _03_thinkstats._mean_var(t)
    sighat = math.sqrt(S2)

    print(xbar, sighat, sighat / xbar)

    # compute standard error for mu and the range of ms
    stderr_xbar = sighat / math.sqrt(n)
    mspread = spread * stderr_xbar
    ms = numpy.linspace(xbar - mspread, xbar + mspread, num_points)

    # compute standard error for sigma and the range of ss
    stderr_sighat = sighat / math.sqrt(2 * (n - 1))
    sspread = spread * stderr_sighat
    ss = numpy.linspace(sighat - sspread, sighat + sspread, num_points)

    # populate the PMF
    pmf = _04_Pmf.Pmf(name=label)
    for m in ms:
        for s in ss:
            pmf._set((m, s), 1)
    return ms, ss, pmf


def _log_update(suite, evidence):
    """
    Updates a suite of hypotheses based on new evidence.

    Modifies the suite directly; if you want to keep the original, make a copy.

    Args:
        suite:    Pmf object
        evidence: whatever kind of object Likelihood expects
    """
    for hypo in suite._values():
        likelihood = _log_likelihood(evidence, hypo)
        suite._incr(hypo, likelihood)
    print(suite._total())


def _log_likelihood(evidence, hypo):
    """
    Computes the log likelihood of the evidence under the hypothesis.

    Args:
        evidence: a list of values
        hypo:     tuple of hypothetical mu and sigma

    Returns:
        log likelihood of the sample given mu and sigma
    """
    t = evidence
    mu, sigma = hypo

    total = _summation(t, mu)
    return -len(t) * math.log(sigma) - total / 2 / sigma ** 2


def _summation(t, mu, cache={}):
    """
    Computes the sum of (x-mu)**2 for x in t.

    Caches previous results.

    Args:
        t:     tuple of values
        mu:    hypothetical mean
        cache: cache of previous results
    """
    try:
        return cache[t, mu]
    except KeyError:
        ds = [(x - mu) ** 2 for x in t]
        total = sum(ds)
        cache[t, mu] = total
        return total


def _estimate_parameters(t, label, num_points=31):
    """
    Computes the posterior distibution of mu and sigma.

    Args:
        t:          sequence of values
        label:      string label for the suite of hypotheses
        num_points: number of values in each dimension

    Returns
        xs:    sequence of hypothetical values for mu
        ys:    sequence of hypothetical values for sigma
        suite: Pmf that maps from (mu, sigma) to prob
    """
    xs, ys, suite = _make_uniform_prior(t, num_points, label)

    suite._log()
    _log_update(suite, tuple(t))
    suite._exp()
    suite._normalize()

    return xs, ys, suite


def _compute_marginals(suite):
    """
    Computes the marginal distributions for mu and sigma.

    Args:
        suite: Pmf that maps (x, y) to z

    Returns:
        Pmf objects for mu and sigma
    """
    pmf_m = _04_Pmf.Pmf()
    pmf_s = _04_Pmf.Pmf()
    for (m, s), p in suite._items():
        pmf_m._incr(m, p)
        pmf_s._incr(s, p)
    return pmf_m, pmf_s


def _compute_coef_variation(suite):
    """
    Computes the distribution of CV.

    Args:
        suite: Pmf that maps (x, y) to z

    Returns:
        Pmf object for CV.
    """
    pmf = _04_Pmf.Pmf()
    for (m, s), p in suite._items():
        pmf._incr(s / m, p)
    return pmf


def _prob_bigger(pmf1, pmf2):
    """Returns the probability that a value from one pmf exceeds another."""
    total = 0.0
    for v1, p1 in pmf1._items():
        for v2, p2 in pmf2._items():
            if v1 > v2:
                total += p1 * p2
    return total


def _plot_posterior(xs, ys, suite, pcolor=False, contour=True):
    """
    Makes a contour plot.

    Args:
        xs:    sequence of values
        ys:    sequence of values
        suite: Pmf that maps (x, y) to z
    """
    X, Y = numpy.meshgrid(xs, ys)
    func = lambda x, y: suite._prob((x, y))
    prob = numpy.vectorize(func)
    Z = prob(X, Y)

    pyplot.clf()
    if pcolor:
        pyplot.pcolor(X, Y, Z)
    if contour:
        pyplot.contour(X, Y, Z)

    _05_myplot._save(root='bayes_height_posterior_%s' % suite.name,
                     title='Posterior joint distribution',
                     xlabel='Mean height (cm)',
                     ylabel='Stddev (cm)')


def _plot_coef_variation(suites):
    """
    Plot the posterior distributions for CV.

    Args:
        suites: map from label to Pmf of CVs.
    """
    pyplot.clf()

    pmfs = {}
    for label, suite in suites.iteritems():
        pmf = _compute_coef_variation(suite)
        cdf = _13_Cdf._make_cdf_from_pmf(pmf, label)
        _05_myplot._cdf(cdf)

        pmfs[label] = pmf

    _05_myplot._save(root='bayes_height_cv',
                     title='Coefficient of variation',
                     xlabel='cv',
                     ylabel='CDF')

    print('female bigger', _prob_bigger(pmfs['female'], pmfs['male']))
    print('male bigger', _prob_bigger(pmfs['male'], pmfs['female']))


def _plot_cdfs(samples):
    """Make CDFs showing the distribution of outliers."""
    cdfs = []
    for label, sample in samples.iteritems():
        outliers = [x for x in sample if x < 150]

        cdf = _13_Cdf._make_cdf_from_list(outliers, label)
        cdfs.append(cdf)

    _05_myplot._clf()
    _05_myplot._cdfs(cdfs)
    _05_myplot._save(root='bayes_height_cdfs',
                     title='CDF of height',
                     xlabel='Reported height (cm)',
                     ylabel='CDF')


def _normal_prob_plot(samples):
    """Makes a normal probability plot for each sample in samples."""
    pyplot.clf()

    markers = dict(male='b', female='g')

    for label, sample in samples.iteritems():
        _normal_plot(sample, label, markers[label], jitter=0.0)

    _05_myplot._save(show=True,
                     # root='bayes_height_normal',
                     title='Normal probability plot',
                     xlabel='Standard normal',
                     ylabel='Reported height (cm)')


def _normal_plot(ys, label, color='b', jitter=0.0, **line_options):
    """
    Makes a normal probability plot.
    
    Args:
        ys:           sequence of values
        label:        string label for the plotted line
        color:        color string passed along to pyplot.plot
        jitter:       float magnitude of jitter added to the ys
        line_options: dictionary of options for pyplot.plot        
    """
    n = len(ys)
    xs = [random.gauss(0.0, 1.0) for i in range(n)]
    xs.sort()
    ys = [y + random.uniform(-jitter, +jitter) for y in ys]
    ys.sort()

    inter, slope = correlation._least_squares(xs, ys)
    fit = correlation._fit_line(xs, inter, slope)
    pyplot.plot(*fit, color=color, linewidth=0.5, alpha=0.5)

    pyplot.plot(sorted(xs), sorted(ys),
                color=color,
                marker='.',
                label=label,
                markersize=3,
                alpha=0.1,
                **line_options)


def _plot_marginals(suite):
    """Plot the marginal distributions for a 2-D joint distribution."""
    pmf_m, pmf_s = _compute_marginals(suite)

    pyplot.clf()
    pyplot.figure(1, figsize=(7, 4))

    pyplot.subplot(1, 2, 1)
    cdf_m = _13_Cdf._make_cdf_from_pmf(pmf_m, 'mu')
    _05_myplot._cdf(cdf_m)
    pyplot.xlabel('Mean height (cm)')
    pyplot.ylabel('CDF')

    pyplot.subplot(1, 2, 2)
    cdf_s = _13_Cdf._make_cdf_from_pmf(pmf_s, 'sigma')
    _05_myplot._cdf(cdf_s)
    pyplot.xlabel('Std Dev height (cm)')
    pyplot.ylabel('CDF')

    _05_myplot._save(root='bayes_height_marginals_%s' % suite.name)


def _plot_ages(resp):
    """Plot the distribution of ages."""
    ages = [r.age for r in resp.records]
    cdf = _13_Cdf._make_cdf_from_list(ages)
    _05_myplot._clf()
    _05_myplot._cdf(cdf)
    _05_myplot._show()


def _dump_heights(data_dir='.', n=10000):
    """Read the BRFSS dataset, extract the heights and pickle them."""
    resp = _19_brfss.Respondents()
    resp._read_records(data_dir, n)

    # PlotAges(resp)

    d = {1: [], 2: []}
    [d[r.sex].append(r.htm3) for r in resp.records if r.htm3 != 'NA']

    fp = open('bayes_height_data.pkl', 'wb')
    cPickle.dump(d, fp)
    fp.close()


def _load_heights():
    """Read the pickled height data."""
    fp = open('bayes_height_data.pkl', 'r')
    d = cPickle.load(fp)
    fp.close()
    return d


def _winsorize(xs, p=0.01):
    """Compresses outliers."""
    cdf = _13_Cdf._make_cdf_from_list(xs)
    low, high = cdf._value(p), cdf._value(1 - p)
    print(low, high)

    outliers = [x for x in xs if x < low or x > high]
    outliers.sort()
    print(outliers)

    wxs = [min(max(low, x), high) for x in xs]
    return wxs


def main():
    if False:
        random.seed(16)
        t = [random.gauss(3, 5) for i in range(100000)]
        _estimate_parameters(t)
        return

    # DumpHeights(n=1000000)
    d = _load_heights()

    labels = {1: 'male', 2: 'female'}

    samples = {}
    suites = {}
    for key, t in d.iteritems():
        label = labels[key]
        print(label, len(t))

        t = _winsorize(t, 0.0001)
        samples[label] = t

        xs, ys, suite = _estimate_parameters(t, label)
        suites[label] = suite

        _plot_posterior(xs, ys, suite)
        _plot_marginals(suite)

    # PlotCdfs(samples)
    # NormalProbPlot(samples)
    _plot_coef_variation(suites)


if __name__ == '__main__':
    main()
