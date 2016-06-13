"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: _08_conditional.py
"""

import matplotlib.pyplot as pyplot
import _04_Pmf
import _05_myplot
import _06_descriptive
import _07_risk


def _condition_pmf(pmf, filter_func, name='conditional'):
    """
    Computes a conditional PMF based on a filter function.
    
    Args:
        pmf:         Pmf object
        filter_func: callable that takes a value from the Pmf and returns a boolean
        name:        string name for the new pmf
        
    Returns:
        new Pmf object
    """
    cond_pmf = pmf._copy(name)

    vals = [val for val in pmf._values() if filter_func(val)]
    for val in vals:
        cond_pmf._remove(val)

    cond_pmf._normalize()
    return cond_pmf


def _condition_on_weeks(pmf, week=39, name='conditional'):
    """
    Computes a PMF conditioned on the given number of weeks.
    
    Args:
        pmf:  Pmf object
        week: the current duration of the pregnancy
        name: string name for the new pmf
        
    Returns:
        new Pmf object
    """

    def filter_func(x):
        return x < week

    cond = _condition_pmf(pmf, filter_func, name)
    return cond


def _make_figure(firsts, others):
    """Makes a figure showing..."""

    weeks = range(35, 46)

    # probs is a map from table name to list of conditional probabilities
    probs = {}
    for table in [firsts, others]:
        name = table.pmf.name
        probs[name] = []
        for week in weeks:
            cond = _condition_on_weeks(table.pmf, week)
            prob = cond._prob(week)
            print(week, prob, table.pmf.name)
            probs[name].append(prob)

    # make a plot with one line for each table
    pyplot.clf()
    # for name, ps in probs.iteritems():
    for name, ps in probs.items():
        pyplot.plot(weeks, ps, label=name)
        print(name, ps)

    _05_myplot._save(root='conditional',
                     xlabel='weeks',
                     ylabel=r'Prob{x $=$ weeks | x $\geq$ weeks}',
                     title='Conditional Probability')


def _relative_risk(first, others, week=38):
    """
    Computes relative risk of the conditional prob of having
    a baby for each week, first babies compared to others.

    Args:
        first:  Pregnancies table
        others: Pregnancies table
        week:
    """
    print(type(first))
    first_cond = _condition_on_weeks(first.pmf, week, 'first babies')
    other_cond = _condition_on_weeks(others.pmf, week, 'others')

    _07_risk._compute_relative_risk(first_cond, other_cond)


def main():
    pool, firsts, others = _06_descriptive._make_tables()
    _relative_risk(firsts, others)
    _make_figure(firsts, others)


if __name__ == "__main__":
    main()
