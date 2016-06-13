"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import _04_Pmf
import _06_descriptive


def _prob_range(pmf, low, high):
    """
    Computes the total probability between low and high, inclusive.
    
    Args:
        pmf:  Pmf object
        low:  low value
        high: high ValueError
        
    Returns:
        float probability
    """
    total = 0.0
    for week in range(low, high + 1):
        total += pmf._prob(week)
    return total


def _prob_early(pmf):
    """
    Computes the probability of a birth in Week 37 or earlier.
    
    Args:
        pmf: Pmf object
        
    Returns:
        float probability
    """
    return _prob_range(pmf, 0, 37)


def _prob_on_time(pmf):
    """
    Computes the probability of a birth in Weeks 38, 39 and 40.
    
    Args:
        pmf: Pmf object
        
    Returns:
        float probability
    """
    return _prob_range(pmf, 38, 40)


def _prob_late(pmf):
    """
    Computes the probability of a birth in Week 41 or later.
    
    Args:
        pmf: Pmf object
        
    Returns:
        float probability
    """
    return _prob_range(pmf, 41, 50)


def _compute_relative_risk(first_pmf, other_pmf):
    """
    Computes relative risks for two PMFs.

    Args:
        first_pmf: Pmf object
        other_pmf: Pmf object
    """

    print('Risks:')
    funcs = [_prob_early, _prob_on_time, _prob_late]
    risks = {}
    for func in funcs:
        for pmf in [first_pmf, other_pmf]:
            prob = func(pmf)
            risks[func.__name__, pmf.name] = prob
            print(func.__name__, pmf.name, prob)

    print
    print('Risk ratios (first babies / others):')
    for func in funcs:
        try:
            ratio = (risks[func.__name__, 'first babies'] / risks[func.__name__, 'others'])
            print(func.__name__, ratio)
        except ZeroDivisionError:
            pass


def main():
    pool, firsts, others = _06_descriptive._make_tables()
    _compute_relative_risk(firsts.pmf, others.pmf)


if __name__ == "__main__":
    main()
