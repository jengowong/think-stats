"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib.pyplot as pyplot
import math
import _01_survey
import _02_first
import _03_thinkstats
import _04_Pmf
import _05_myplot
import _06_descriptive
import _13_Cdf


def Process(table, name):
    """
    Runs various analyses on this table.

    Creates instance variables:
        weights:    sequence of int total weights in ounces
        weight_pmf: Pmf object
        weight_cdf: Cdf object
        oz_pmf:     Pmf of just the ounce field
    """
    _06_descriptive._process(table, name)

    table.weights = [p.totalwgt_oz for p in table.records if p.totalwgt_oz != 'NA']
    table.weight_pmf = _04_Pmf._make_pmf_from_list(table.weights, table.name)
    table.weight_cdf = _13_Cdf._make_cdf_from_list(table.weights, table.name)


def MakeTables(data_dir='.'):
    """Reads survey data and returns a tuple of Tables"""
    table, firsts, others = _02_first._make_tables(data_dir)
    pool = _06_descriptive._pool_records(firsts, others)

    Process(pool, 'live births')
    Process(firsts, 'first babies')
    Process(others, 'others')

    return pool, firsts, others


def Resample(cdf, n=10000):
    sample = cdf._sample(n)
    new_cdf = _13_Cdf._make_cdf_from_list(sample, 'resampled')
    _05_myplot._clf()
    _05_myplot._cdfs([cdf, new_cdf])
    _05_myplot._save(root='resample_cdf',
                     title='CDF',
                     xlabel='weight in oz',
                     ylabel='CDF(x)')


def MakeExample():
    """Make a simple example CDF."""
    t = [2, 1, 3, 2, 5]
    cdf = _13_Cdf._make_cdf_from_list(t)
    _05_myplot._clf()
    _05_myplot._cdf(cdf)
    _05_myplot._save(root='example_cdf',
                     title='CDF',
                     xlabel='x',
                     ylabel='CDF(x)',
                     axis=[0, 6, 0, 1],
                     legend=False)


def MakeFigures(pool, firsts, others):
    """Creates several figures for the book."""

    # plot PMFs of birth weights for first babies and others
    _05_myplot._clf()
    _05_myplot._hist(firsts.weight_pmf, linewidth=0, color='blue')
    _05_myplot._hist(others.weight_pmf, linewidth=0, color='orange')
    _05_myplot._save(root='nsfg_birthwgt_pmf',
                     title='Birth weight PMF',
                     xlabel='weight (ounces)',
                     ylabel='probability')

    # plot CDFs of birth weights for first babies and others
    _05_myplot._clf()
    _05_myplot._cdf(firsts.weight_cdf, linewidth=2, color='blue')
    _05_myplot._cdf(others.weight_cdf, linewidth=2, color='orange')
    _05_myplot._save(root='nsfg_birthwgt_cdf',
                     title='Birth weight CDF',
                     xlabel='weight (ounces)',
                     ylabel='probability',
                     axis=[0, 200, 0, 1])


def main(name, data_dir=''):
    MakeExample()

    pool, firsts, others = MakeTables(data_dir)
    MakeFigures(pool, firsts, others)


if __name__ == '__main__':
    import sys

    main(*sys.argv)
