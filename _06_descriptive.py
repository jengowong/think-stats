"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: _06_descriptive.py
"""

import math
import matplotlib.pyplot as pyplot
import _01_survey
import _02_first
import _03_thinkstats
import _04_Pmf
import _05_myplot


def _process(table, name):
    """Runs various analyses on this table."""
    _02_first._process(table)
    table.name = name

    table.var = _03_thinkstats._var(table.lengths, table.mu)
    table.trim = _03_thinkstats._trimmed_mean(table.lengths)

    table.hist = _04_Pmf._make_hist_from_list(table.lengths, name=name)
    table.pmf = _04_Pmf._make_pmf_from_hist(table.hist)


def _pool_records(*tables):
    """
    Construct a table with records from all tables.
    
    Args:
        constructor: init method used to make the new table
    
        tables: any number of tables

    Returns:
        new table object
    """
    pool = _01_survey.Pregnancies()
    for table in tables:
        pool._extend_records(table.records)
    return pool


def _make_tables(data_dir='.'):
    """Reads survey data and returns a tuple of Tables"""
    table, firsts, others = _02_first._make_tables(data_dir)
    pool = _pool_records(firsts, others)

    _process(pool, 'live births')
    _process(firsts, 'first babies')
    _process(others, 'others')

    return pool, firsts, others


def _summarize(pool, firsts, others):
    """Print various summary statistics."""

    print
    print('Variance')
    print('First babies', firsts.var)
    print('Others', others.var)

    diff_mu = firsts.mu - others.mu

    print('Difference in mean', diff_mu)

    sigma = math.sqrt(pool.var)

    print('Pooled mean', pool.mu)
    print('Pooled variance', pool.var)
    print('Pooled sigma', sigma)

    print(firsts.mu, others.mu)
    print(firsts.trim, others.trim)

    # live_lengths = pool.hist._get_dict().items()
    live_lengths = list(pool.hist._get_dict().items())
    # live_lengths.sort()
    sorted(live_lengths)
    print('Shortest lengths:')
    for weeks, count in live_lengths[:10]:
        print(weeks, count)

    print('Longest lengths:')
    for weeks, count in live_lengths[-10:]:
        print(weeks, count)


def _make_figures(firsts, others):
    """Plot Hists and Pmfs for the pregnancy length."""

    # bar options is a list of option dictionaries to be passed to myplot.bar
    bar_options = [
        dict(color='0.9'),
        dict(color='blue')
    ]

    # make the histogram
    axis = [23, 46, 0, 2700]
    _hists([firsts.hist, others.hist])
    _05_myplot._save(root='nsfg_hist',
                     title='Histogram',
                     xlabel='weeks',
                     ylabel='frequency',
                     axis=axis)

    # make the PMF
    axis = [23, 46, 0, 0.6]
    _hists([firsts.pmf, others.pmf])
    _05_myplot._save(root='nsfg_pmf',
                     title='PMF',
                     xlabel='weeks',
                     ylabel='probability',
                     axis=axis)


def _hists(hists):
    """
    Plot two histograms on the same axes.

    hists: list of Hist
    """
    width = 0.4
    shifts = [-width, 0.0]

    option_list = [
        dict(color='0.9'),
        dict(color='blue')
    ]

    pyplot.clf()
    for i, hist in enumerate(hists):
        xs, fs = hist._render()
        xs = _shift(xs, shifts[i])
        pyplot.bar(xs, fs, label=hist.name, width=width, **option_list[i])


def _shift(xs, shift):
    """
    Adds a constant to a sequence of values.

    Args:
        xs:    sequence of values
        shift: value to add

    Returns:
        sequence of numbers
    """
    return [x + shift for x in xs]


def _make_diff_figure(firsts, others):
    """Plot the difference between the PMFs."""

    weeks = range(35, 46)
    diffs = []
    for week in weeks:
        p1 = firsts.pmf._prob(week)
        p2 = others.pmf._prob(week)
        diff = 100 * (p1 - p2)
        diffs.append(diff)

    pyplot.clf()
    pyplot.bar(weeks, diffs, align='center')
    _05_myplot._save(root='nsfg_diffs',
                     title='Difference in PMFs',
                     xlabel='weeks',
                     ylabel='100 (PMF$_{first}$ - PMF$_{other}$)',
                     legend=False)


def main(name, data_dir=''):
    pool, firsts, others = _make_tables(data_dir)
    _summarize(pool, firsts, others)
    _make_figures(firsts, others)
    _make_diff_figure(firsts, others)


if __name__ == '__main__':
    import sys

    main(*sys.argv)
