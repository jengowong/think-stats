"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: pareto_world.py
"""

import matplotlib.pyplot as pyplot
import random
import _05_myplot
import _13_Cdf


def _pareto_cdf(x, alpha, xmin):
    """Evaluates CDF of the Pareto distribution with parameters alpha, xmin."""
    if x < xmin:
        return 0
    return 1 - pow(x / xmin, -alpha)


def _pareto_median(xmin, alpha):
    """Computes the median of a Pareto distribution."""
    return xmin * pow(2, 1 / alpha)


def _make_pareto_cdf():
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
    _05_myplot._save('pareto_world1',
                     title='Pareto CDF',
                     xlabel='height (cm)',
                     ylabel='CDF',
                     legend=False)


def _make_figure(xmin=100, alpha=1.7, mu=150, sigma=25):
    """
    Makes a figure showing the CDF of height in ParetoWorld.

    Compared to a normal distribution.

    Args:
        xmin:  parameter of the Pareto distribution
        alpha: parameter of the Pareto distribution
        mu:    parameter of the Normal distribution
        sigma: parameter of the Normal distribution
    """
    t1 = [xmin * random.paretovariate(alpha) for i in range(10000)]
    cdf1 = _13_Cdf._make_cdf_from_list(t1, name='pareto')

    t2 = [random.normalvariate(mu, sigma) for i in range(10000)]
    cdf2 = _13_Cdf._make_cdf_from_list(t2, name='normal')

    _05_myplot._clf()
    _05_myplot._cdfs([cdf1, cdf2])
    _05_myplot._save(root='pareto_world2',
                     title='Pareto World',
                     xlabel='height (cm)',
                     ylabel='CDF')


def _tallest_pareto(iters=2, n=10000, xmin=100, alpha=1.7):
    """
    Find the tallest person in Pareto World.

    Args:
        iters: how many samples to generate
        n:     how many in each sample
        xmin:  parameter of the Pareto distribution
        alpha: parameter of the Pareto distribution
    """
    tallest = 0
    for i in range(iters):
        t = [xmin * random.paretovariate(alpha) for i in range(n)]
        tallest = max(max(t), tallest)
    return tallest


def main():
    _make_figure()
    _make_pareto_cdf()
    print(_tallest_pareto(iters=2))


if __name__ == "__main__":
    main()
