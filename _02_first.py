"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: _02_first.py
"""

import _01_survey


# copying Mean from _03_thinkstats.py so we don't have to deal with
# importing anything in Chapter 1

def _mean(t):
    """
    Computes the mean of a sequence of numbers.

    Args:
        t: sequence of numbers

    Returns:
        float
    """
    return float(sum(t)) / len(t)


def _partition_records(table):
    """
    Divides records into two lists: first babies and others.

    Only live births are included

    Args:
        table: pregnancy Table
    """
    firsts = _01_survey.Pregnancies()
    others = _01_survey.Pregnancies()

    for p in table.records:
        # skip non-live births
        if p.outcome != 1:
            continue

        if p.birthord == 1:
            firsts._add_record(p)
        else:
            others._add_record(p)

    return firsts, others


def _process(table):
    """
    Runs analysis on the given table.
    
    Args:
        table: table object
    """
    table.lengths = [p.prglength for p in table.records]
    table.n = len(table.lengths)
    table.mu = _mean(table.lengths)


def _make_tables(data_dir='.'):
    """Reads survey data and returns tables for first babies and others."""
    table = _01_survey.Pregnancies()
    table._read_records(data_dir)

    firsts, others = _partition_records(table)

    return table, firsts, others


def _process_tables(*tables):
    """
    Processes a list of tables
    
    Args:
        tables: gathered argument tuple of Tuples
    """
    for table in tables:
        _process(table)


def _summarize(data_dir):
    """
    Prints summary statistics for first babies and others.
    
    Returns:
        tuple of Tables
    """
    table, firsts, others = _make_tables(data_dir)
    _process_tables(firsts, others)

    print('Number of first babies', firsts.n)
    print('Number of others', others.n)

    mu1, mu2 = firsts.mu, others.mu

    print('Mean gestation in weeks:')
    print('First babies', mu1)
    print('Others', mu2)

    print('Difference in days', (mu1 - mu2) * 7.0)


def main(name, data_dir='.'):
    _summarize(data_dir)


if __name__ == '__main__':
    import sys

    main(*sys.argv)
