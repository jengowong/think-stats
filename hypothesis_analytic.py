"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: hypothesis_analytic.py
"""

import matplotlib.pyplot as pyplot
import cumulative
import hypothesis
import math
import _03_thinkstats
import _16_erf


def _test(actual1, actual2, model, iters=1000):
    """
    Estimates p-values based on differences in the mean.
    
    Args:
        actual1:
        actual2: sequences of observed values for groups 1 and 2
        model:   sequences of values from the hypothetical distribution
    """
    n = len(actual1)
    m = len(actual2)

    # compute delta
    mu1, mu2, delta = hypothesis._difference_in_mean(actual1, actual2)
    delta = abs(delta)

    print('n:', n)
    print('m:', m)
    print('mu1', mu1)
    print('mu2', mu2)
    print('delta', delta)

    # compute the expected distribution of differences in sample mean
    mu_pooled, var_pooled = _03_thinkstats._mean_var(model)
    print('(Mean, Var) of pooled data', mu_pooled, var_pooled)

    f = 1.0 / n + 1.0 / m
    mu, var = (0, f * var_pooled)
    print('Expected Mean, Var of deltas', mu, var)

    # compute the p-value of delta in the observed distribution
    sigma = math.sqrt(var)
    left = _16_erf._normal_cdf(-delta, mu, sigma)
    right = 1 - _16_erf._normal_cdf(delta, mu, sigma)
    pvalue = left + right
    print('Tails:', left, right)
    print('p-value:', pvalue)

    # compare the mean and variance of resamples differences
    deltas = [hypothesis._resample(model, model, n, m) for i in range(iters)]
    mean_var = _03_thinkstats._mean_var(deltas)
    print('(Mean, Var) of resampled deltas', mean_var)

    return pvalue


def main():
    # get the data
    pool, firsts, others = cumulative._make_tables()

    # run the test
    _test(firsts.lengths,
          others.lengths,
          pool.lengths,
          iters=1000)


if __name__ == "__main__":
    main()
