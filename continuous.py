"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: continuous.py
"""

import matplotlib.pyplot as pyplot
import cumulative
import math
import random
import _03_thinkstats
import _05_myplot
import _13_Cdf
import _16_erf
import _17_rankit


def _expo_cdf(x, lam):
    """Evaluates CDF of the exponential distribution with parameter lam."""
    return 1 - math.exp(-lam * x)


def _pareto_cdf(x, alpha, xmin):
    """Evaluates CDF of the Pareto distribution with parameters alpha, xmin."""
    if x < xmin:
        return 0
    return 1 - pow(x / xmin, -alpha)


def _pareto_median(xmin, alpha):
    """Computes the median of a Pareto distribution."""
    return xmin * pow(2, 1 / alpha)


def _make_expo_cdf():
    """Generates a plot of the exponential CDF."""
    n = 40
    max = 2.5
    xs = [max * i / n for i in range(n)]

    lam = 2.0
    ps = [_expo_cdf(x, lam) for x in xs]

    percentile = -math.log(0.05) / lam
    print('Fraction <= ', percentile, _expo_cdf(lam, percentile))

    pyplot.clf()
    pyplot.plot(xs, ps, linewidth=2)
    _05_myplot._save('expo_cdf',
                     title='Exponential CDF',
                     xlabel='x',
                     ylabel='CDF',
                     legend=False)


def _make_pareto_cdf():
    """Generates a plot of the Pareto CDF."""
    n = 50
    max = 10.0
    xs = [max * i / n for i in range(n)]

    xmin = 0.5
    alpha = 1.0
    ps = [_pareto_cdf(x, alpha, xmin) for x in xs]
    print('Fraction <= 10', _pareto_cdf(xmin, alpha, 10))

    pyplot.clf()
    pyplot.plot(xs, ps, linewidth=2)
    _05_myplot._save('pareto_cdf',
                     title='Pareto CDF',
                     xlabel='x',
                     ylabel='CDF',
                     legend=False)


def _make_pareto_cdf2():
    """Generates a plot of the CDF of height in Pareto World."""
    n = 50
    max = 1000.0
    xs = [max * i / n for i in range(n)]

    xmin = 100
    alpha = 1.7
    ps = [_pareto_cdf(x, alpha, xmin) for x in xs]
    print('Median', _pareto_median(xmin, alpha))

    pyplot.clf()
    pyplot.plot(xs, ps, linewidth=2)
    _05_myplot._save('pareto_height',
                     title='Pareto CDF',
                     xlabel='height (cm)',
                     ylabel='CDF',
                     legend=False)


def _render_normal_cdf(mu, sigma, max, n=50):
    """Generates sequences of xs and ps for a normal CDF."""
    xs = [max * i / n for i in range(n)]
    ps = [_16_erf._normal_cdf(x, mu, sigma) for x in xs]
    return xs, ps


def _make_normal_cdf():
    """Generates a plot of the normal CDF."""
    xs, ps = _render_normal_cdf(2.0, 0.5, 4.0)

    pyplot.clf()
    pyplot.plot(xs, ps, linewidth=2)
    _05_myplot._save('normal_cdf',
                     title='Normal CDF',
                     xlabel='x',
                     ylabel='CDF',
                     legend=False)


def _make_normal_model(weights):
    """Plot the CDF of birthweights with a normal model."""

    # estimate parameters: trimming outliers yields a better fit
    mu, var = _03_thinkstats._trimmed_mean_var(weights, p=0.01)
    print('Mean, Var', mu, var)

    # plot the model
    sigma = math.sqrt(var)
    print('Sigma', sigma)
    xs, ps = _render_normal_cdf(mu, sigma, 200)

    pyplot.clf()
    pyplot.plot(xs, ps, label='model', linewidth=4, color='0.8')

    # plot the data
    cdf = _13_Cdf._make_cdf_from_list(weights)
    xs, ps = cdf._render()
    pyplot.plot(xs, ps, label='data', linewidth=2, color='blue')

    _05_myplot._save('nsfg_birthwgt_model',
                     title='Birth weights',
                     xlabel='birth weight (oz)',
                     ylabel='CDF')


def _make_normal_plot(weights):
    """Generates a normal probability plot of birth weights."""
    _17_rankit._make_normal_plot(weights,
                                 root='nsfg_birthwgt_normal',
                                 ylabel='Birth weights (oz)', )


def main():
    random.seed(17)

    # make the continuous CDFs
    _make_expo_cdf()
    _make_pareto_cdf()
    _make_pareto_cdf2()
    _make_normal_cdf()

    # test the distribution of birth weights for normality
    pool, _, _ = cumulative._make_tables()

    t = pool.weights
    _make_normal_model(t)
    _make_normal_plot(t)


if __name__ == "__main__":
    main()
