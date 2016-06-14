"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: khan.py
"""

import matplotlib.pyplot as pyplot
import scipy
import scipy.stats
import random
import numpy
import _03_thinkstats
import _04_Pmf
import _05_myplot
import _13_Cdf


def _chi_squared(expected, observed):
    """
    Compute the Chi-squared statistic for two tables.
    
    Args:
        expected: Hist of expected values
        observed: Hist of observed values
      
    Returns:
        float chi-squared statistic
    """
    total = 0.0
    for x, exp in expected._items():
        obs = observed._freq(x)
        total += (obs - exp) ** 2 / exp
    return total


def _simulate(pa, q, n):
    """
    Run a simulation.

    Args:
        pa: probability of showing version A
        q:  probability of success for both A and B
        n:  number of trials
    
    Returns: tuple of
        diff_list:   simulated differences
        chi2_list:   simulated chi-square values
        pvalue_list: simulated p-values
    """
    hist = _04_Pmf.Hist()
    diff_list = []
    chi2_list = []
    pvalue_list = []

    for i in range(1, n + 1):
        version = _flip(pa, 'A', 'B')
        outcome = _flip(q, 'Y', 'N')
        hist._incr((version, outcome))

        expected = _expected(pa, q, i)
        try:
            rate_a, rate_b = _compute_rates(hist)
            diff = rate_b - rate_a
            chi2 = _chi_squared(expected, hist)
            pvalue = _pvalue(chi2, df=3)

            diff_list.append((i, diff))
            chi2_list.append((i, chi2))
            pvalue_list.append((i, pvalue))
        except ZeroDivisionError:
            pass

    return diff_list, chi2_list, pvalue_list


def _flip(p, y='Y', n='N'):
    """Returns y with probability p; otherwise n."""
    return y if random.random() <= p else n


def _compute_rates(hist):
    """Returns a sequence of success rates, one for each version."""
    rates = []
    for version in ['A', 'B']:
        hits = hist._freq((version, 'Y'))
        misses = hist._freq((version, 'N'))
        rate = float(hits) / (hits + misses)
        rates.append(rate)

    return rates


def _expected(pa, q, n):
    """
    Makes a Pmf with the expected number of trials in each of four bins.

    Args:
        pa: probability of offering version A
        q:  probability of success for both A and B
        n:  number of trials

    Returns:
        hist that maps (version, outcome) to expected number,
        where version is string A or B and outcome is string Y or N.
    """
    versions = _04_Pmf._make_pmf_from_dict(dict(A=pa, B=1 - pa))
    outcomes = _04_Pmf._make_pmf_from_dict(dict(Y=q, N=1 - q))

    hist = _04_Pmf.Hist()
    for version, pp in versions._items():
        for outcome, qq in outcomes._items():
            hist._incr((version, outcome), pp * qq * n)
    return hist


def _pvalue(chi2, df):
    """
    Returns the p-value of getting chi2 from a chi-squared distribution.

    Args:
        chi2: observed chi-squared statistic
        df:   degrees of freedom
    """
    return 1 - scipy.stats.chi2.cdf(chi2, df)


def _crosses(ps, thresh):
    """Tests whether a sequence of p-values ever drops below thresh."""
    if thresh is None:
        return False

    for p in ps:
        if p <= thresh:
            return True
    return False


def _check_cdf():
    """Compare chi2 values from simulation with chi2 distributions."""
    for df in [1, 2, 3]:
        xs, ys = _chi2_cdf(df=df, high=15)
        pyplot.plot(xs, ys, label=df)

    t = [_simulate_chi2() for i in range(1000)]
    cdf = _13_Cdf._make_cdf_from_list(t)

    _05_myplot._cdf(cdf)
    _05_myplot._save(root='khan3',
                     xlabel='chi2 value',
                     ylabel="CDF",
                     formats=['png'])


def _check_cdf2():
    """Compare chi2 values from the simulation with a chi-squared dist."""
    df = 3
    t = [_simulate_chi2() for i in range(1000)]
    t2 = [scipy.stats.chi2.cdf(x, df) for x in t]
    cdf = _13_Cdf._make_cdf_from_list(t2)

    _05_myplot._cdf(cdf)
    _05_myplot._show()


def _chi2_cdf(df=2, high=5, n=100):
    """
    Evaluates the chi-squared CDF.

    Args:
        df:   degrees of freedom
        high: high end of the range of x
        n:    number of values between 0 and high
    """
    xs = numpy.linspace(0, high, n)
    ys = scipy.stats.chi2.cdf(xs, df)
    return xs, ys


def _simulate_chi2(pa=0.5, q=0.5, n=100):
    """Run a simulation and return the chi2 statistic."""
    expected = _expected(pa, q, n)

    simulated = _04_Pmf.Hist()

    for i in range(1, n + 1):
        version = _flip(pa, 'A', 'B')
        outcome = _flip(q, 'Y', 'N')
        simulated._incr((version, outcome))

    chi2 = _chi_squared(expected, simulated)
    return chi2


def _make_spaghetti(iters=1000, lines=100, n=300, thresh=0.05, index=2):
    """
    Makes a spaghetti plot of random-walk lines.

    Args:
        iters:  number of simulations to run
        lines:  number of lines to plot
        n:      number of trials to simulate
        thresh: threshold p-value
    """
    pyplot.clf()
    if thresh is not None:
        pyplot.plot([1, n], [thresh, thresh], color='red', alpha=1, linewidth=2)

    count = 0.0
    for i in range(iters):
        lists = _simulate(0.5, 0.5, n)
        pairs = lists[index]
        xs, ys = zip(*pairs)
        if _crosses(ys, thresh):
            count += 1

        if i < lines:
            pyplot.plot(xs, ys, alpha=0.2)

    print(iters, count / iters)

    labels = ['Difference in success rate', 'chi-squared stat', 'p-value']

    _05_myplot._save(root='khan%d' % index,
                     xlabel='Number of trials',
                     ylabel=labels[index],
                     title='A-B test random walk',
                     formats=['png'])


def main():
    _check_cdf()
    return

    random.seed(17)
    _make_spaghetti(10, 10, 200, index=0, thresh=None)
    _make_spaghetti(10, 10, 1000, index=1, thresh=None)
    _make_spaghetti(20, 20, 1000, index=2, thresh=0.05)


if __name__ == "__main__":
    main()
