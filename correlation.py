"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: correlation.py
"""

import math
import random
import _03_thinkstats


def _cov(xs, ys, mux=None, muy=None):
    """
    Computes Cov(X, Y).

    Args:
        xs:  sequence of values
        ys:  sequence of values
        mux: optional float mean of xs
        muy: optional float mean of ys

    Returns:
        Cov(X, Y)
    """
    if mux is None:
        mux = _03_thinkstats._mean(xs)
    if muy is None:
        muy = _03_thinkstats._mean(ys)

    total = 0.0
    for x, y in zip(xs, ys):
        total += (x - mux) * (y - muy)

    return total / len(xs)


def _corr(xs, ys):
    """
    Computes Corr(X, Y).

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        Corr(X, Y)
    """
    xbar, varx = _03_thinkstats._mean_var(xs)
    ybar, vary = _03_thinkstats._mean_var(ys)

    corr = _cov(xs, ys, xbar, ybar) / math.sqrt(varx * vary)

    return corr


def _serial_corr(xs):
    """Computes the serial correlation of a sequence."""
    return _corr(xs[:-1], xs[1:])


def _spearman_corr(xs, ys):
    """
    Computes Spearman's rank correlation.

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        float Spearman's correlation
    """
    xranks = _map_to_ranks(xs)
    yranks = _map_to_ranks(ys)
    return _corr(xranks, yranks)


def _least_squares(xs, ys):
    """
    Computes a linear least squares fit for ys as a function of xs.

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        tuple of (intercept, slope)
    """
    xbar, varx = _03_thinkstats._mean_var(xs)
    ybar, vary = _03_thinkstats._mean_var(ys)

    slope = _cov(xs, ys, xbar, ybar) / varx
    inter = ybar - slope * xbar

    return inter, slope


def _fit_line(xs, inter, slope):
    """
    Returns the fitted line for the range of xs.

    Args:
        xs:    x values used for the fit
        slope: estimated slope
        inter: estimated intercept
    """
    fxs = min(xs), max(xs)
    fys = [x * slope + inter for x in fxs]
    return fxs, fys


def _residuals(xs, ys, inter, slope):
    """
    Computes residuals for a linear fit with parameters inter and slope.

    Args:
        xs:    independent variable
        ys:    dependent variable
        inter: float intercept
        slope: float slope

    Returns:
        list of residuals
    """
    res = [y - inter - slope * x for x, y in zip(xs, ys)]
    return res


def _coef_determination(ys, res):
    """
    Computes the coefficient of determination (R^2) for given residuals.

    Args:
        ys:  dependent variable
        res: residuals
        
    Returns:
        float coefficient of determination
    """
    ybar, vary = _03_thinkstats._mean_var(ys)
    resbar, varres = _03_thinkstats._mean_var(res)
    return 1 - varres / vary


def _map_to_ranks(t):
    """
    Returns a list of ranks corresponding to the elements in t.

    Args:
        t: sequence of numbers
    
    Returns:
        list of integer ranks, starting at 1
    """
    # pair up each value with its index
    pairs = enumerate(t)

    # sort by value
    sorted_pairs = sorted(pairs, key=lambda pair: pair[1])

    # pair up each pair with its rank
    ranked = enumerate(sorted_pairs)

    # sort by index
    resorted = sorted(ranked, key=lambda trip: trip[1][0])

    # extract the ranks
    ranks = [trip[0] + 1 for trip in resorted]
    return ranks


def _correlated_generator(rho):
    """
    Generates standard normal variates with correlation.

    Args:
        rho: target coefficient of correlation

    Returns:
        iterable
    """
    x = random.gauss(0, 1)
    yield x

    sigma = math.sqrt(1 - rho ** 2)
    while True:
        x = random.gauss(x * rho, sigma)
        yield x


def _correlated_normal_generator(mu, sigma, rho):
    """
    Generates normal variates with correlation.

    Args:
        mu:    mean of variate
        sigma: standard deviation of variate
        rho:   target coefficient of correlation

    Returns:
        iterable
    """
    for x in _correlated_generator(rho):
        yield x * sigma + mu


def main():
    pass


if __name__ == '__main__':
    main()
