"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

Name: hypothesis.py
"""

import matplotlib.pyplot as pyplot
import cumulative
import random
import math
import _03_thinkstats
import _05_myplot
import _13_Cdf


def _run_test(root,
              pool,
              actual1,
              actual2,
              iters=1000,
              trim=False,
              partition=False):
    """
    Computes the distributions of delta under H0 and HA.
    
    Args:
        root:      string filename root for the plots
        pool:      sequence of values from the pooled distribution
        actual1:   sequence of values in group 1
        actual2:   sequence of values in group 2
        iters:     how many resamples
        trim:      whether to trim the sequences
        partition: whether to cross-validate by partitioning the data
    """
    if trim:
        pool.sort()
        actual1.sort()
        actual2.sort()
        pool = _03_thinkstats._trim(pool)
        actual1 = _03_thinkstats._trim(actual1)
        actual2 = _03_thinkstats._trim(actual2)

    if partition:
        n = len(actual1)
        m = len(actual2)
        actual1, model1 = _partition(actual1, n / 2)
        actual2, model2 = _partition(actual2, m / 2)
        pool = model1 + model2
    else:
        model1 = actual1
        model2 = actual2

    # P(E|H0)
    peh0 = _test(root + '_deltas_cdf',
                 actual1,
                 actual2,
                 pool,
                 pool,
                 iters,
                 plot=True)

    # P(E|Ha)
    peha = _test(root + '_deltas_ha_cdf',
                 actual1,
                 actual2,
                 model1,
                 model2,
                 iters)

    prior = 0.5
    pe = prior * peha + (1 - prior) * peh0
    posterior = prior * peha / pe
    print('Posterior', posterior)


def _test(root, actual1, actual2, model1, model2, iters=1000, plot=False):
    """
    Estimates p-values based on differences in the mean.
    
    Args:
        root:    string filename base for plots
        actual1:
        actual2: sequences of observed values for groups 1 and 2
        model1: 
        model2:  sequences of values from the hypothetical distributions
        iters:   how many resamples
        plot:    whether to plot the distribution of differences in the mean
    """
    n = len(actual1)
    m = len(actual2)

    mu1, mu2, delta = _difference_in_mean(actual1, actual2)
    delta = abs(delta)

    cdf, pvalue = _p_value(model1, model2, n, m, delta, iters)
    print('n:', n)
    print('m:', m)
    print('mu1', mu1)
    print('mu2', mu2)
    print('delta', delta)
    print('p-value', pvalue)

    if plot:
        _plot_cdf(root, cdf, delta)

    return pvalue


def _difference_in_mean(actual1, actual2):
    """
    Computes the difference in mean between two groups.

    Args:
        actual1: sequence of float
        actual2: sequence of float

    Returns:
        tuple of (mu1, mu2, mu1-mu2)
    """
    mu1 = _03_thinkstats._mean(actual1)
    mu2 = _03_thinkstats._mean(actual2)
    delta = mu1 - mu2
    return mu1, mu2, delta


def _p_value(model1, model2, n, m, delta, iters=1000):
    """
    Computes the distribution of deltas with the model distributions.

    And the p-value of the observed delta.

    Args:
        model1: 
        model2: sequences of values from the hypothetical distributions
        n:      sample size from model1
        m:      sample size from model2
        delta:  the observed difference in the means
        iters:  how many samples to generate
    """
    deltas = [_resample(model1, model2, n, m) for i in range(iters)]
    mean_var = _03_thinkstats._mean_var(deltas)
    print('(Mean, Var) of resampled deltas', mean_var)

    cdf = _13_Cdf._make_cdf_from_list(deltas)

    # compute the two tail probabilities
    left = cdf._prob(-delta)
    right = 1.0 - cdf._prob(delta)

    pvalue = left + right
    print('Tails (left, right, total):', left, right, left + right)

    return cdf, pvalue


def _plot_cdf(root, cdf, delta):
    """
    Draws a Cdf with vertical lines at the observed delta.

    Args:
        root:  string used to generate filenames
        cdf:   Cdf object
        delta: float observed difference in means
    """

    def _vert_line(x):
        """Draws a vertical line at x."""
        xs = [x, x]
        ys = [0, 1]
        pyplot.plot(xs, ys, linewidth=2, color='0.7')

    _vert_line(-delta)
    _vert_line(delta)

    xs, ys = cdf._render()
    pyplot.subplots_adjust(bottom=0.11)
    pyplot.plot(xs, ys, linewidth=2, color='blue')

    _05_myplot._save(root,
                     title='Resampled differences',
                     xlabel='difference in means (weeks)',
                     ylabel='CDF(x)',
                     legend=False)


def _resample(t1, t2, n, m):
    """
    Draws samples and computes the difference in mean.
    
    Args:
        t1: sequence of values
        t2: sequence of values
        n:  size of the sample to draw from t1
        m:  size of the sample to draw from t2
    """
    sample1 = _sample_with_replacement(t1, n)
    sample2 = _sample_with_replacement(t2, m)
    mu1, mu2, delta = _difference_in_mean(sample1, sample2)
    return delta


def _partition(t, n):
    """
    Splits a sequence into two random partitions.
    
    Side effect: shuffles t
    
    Args:
        t: sequence of values
        n: size of the first partition

    Returns:
        two lists of values
    """
    random.shuffle(t)
    return t[:n], t[n:]


def _sample_with_replacement(t, n):
    """
    Generates a sample with replacement.
    
    Args:
        t: sequence of values
        n: size of the sample
        
    Returns:
        list of values
    """
    return [random.choice(t) for i in range(n)]


def _sample_without_replacement(t, n):
    """
    Generates a sample without replacement.
    
    Args:
        t: sequence of values
        n: size of the sample
        
    Returns:
        list of values
    """
    return random.sample(t, n)


def main():
    random.seed(1)

    # get the data
    pool, firsts, others = cumulative._make_tables()
    mean_var = _03_thinkstats._mean_var(pool.lengths)
    print('(Mean, Var) of pooled data', mean_var)

    # run the test
    _run_test('length',
              pool.lengths,
              firsts.lengths,
              others.lengths,
              iters=1000,
              trim=False,
              partition=False)


if __name__ == "__main__":
    main()
