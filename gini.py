"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: gini.py
"""

import math
import sys
import _04_Pmf
import _13_Cdf
import _23_irs


def _pmf_mean(pmf):
    total = 0.0
    for val, p in pmf._items():
        total += p * val
    return total


def _pmf_moment(pmf, mean=None, exponent=2):
    if mean is None:
        mean = _pmf_mean(pmf)

    total = 0.0
    for val, p in pmf._items():
        total += p * (val - mean) ** exponent
    return total


def _relative_mean_difference(pmf, mean=None):
    if mean is None:
        mean = _pmf_mean(pmf)

    diff = _04_Pmf.Pmf()
    for v1, p1 in pmf._items():
        for v2, p2 in pmf._items():
            diff._incr(abs(v1 - v2), p1 * p2)

    print(_pmf_mean(diff), mean)

    return _pmf_mean(diff) / mean


def _summarize_data(pmf, cdf):
    mean = _pmf_mean(pmf)
    print('mean:', mean)

    median = cdf._percentile(50)
    print('median:', median)

    fraction_below_mean = cdf._prob(mean)
    print('fraction below mean:', fraction_below_mean)

    m2 = _pmf_moment(pmf, mean, 2)
    m3 = _pmf_moment(pmf, mean, 3)

    sigma = math.sqrt(m2)
    print('sigma:', sigma)

    g1 = m3 / m2 ** (3 / 2)
    print('skewness:', g1)

    gp = 3 * (mean - median) / sigma
    print('Pearsons skewness:', gp)

    gini = _relative_mean_difference(pmf) / 2
    print('gini', gini)


def main(script, *args):
    data = _23_irs._read_income_file()
    hist, pmf, cdf = _23_irs._make_income_dist(data)
    _summarize_data(pmf, cdf)


if __name__ == "__main__":
    main(*sys.argv)
