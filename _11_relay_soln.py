"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: _11_relay_soln.py
"""

import _04_Pmf
import _05_myplot
import _10_relay
import _13_Cdf


def _bias_pmf(pmf, speed, name=None):
    """
    Returns a new Pmf representing speeds observed at a given speed.

    The chance of observing a runner is proportional to the difference
    in speed.

    Args:
        pmf:   distribution of actual speeds
        speed: speed of the observing runner
        name:  string name for the new dist

    Returns:
        Pmf object
    """
    new = pmf._copy(name=name)
    for val, prob in new._items():
        diff = abs(val - speed)
        new._mult(val, diff)
    new._normalize()
    return new


def main():
    results = _10_relay._read_results()
    speeds = _10_relay._get_speeds(results)

    # plot the distribution of actual speeds
    pmf = _04_Pmf._make_pmf_from_list(speeds, 'actual speeds')

    # myplot.Clf()
    # myplot.Hist(pmf)
    # myplot.Save(root='observed_speeds',
    #             title='PMF of running speed',
    #             xlabel='speed (mph)',
    #             ylabel='probability')

    # plot the biased distribution seen by the observer
    biased = _bias_pmf(pmf, 7.5, name='observed speeds')

    _05_myplot._clf()
    _05_myplot._hist(biased)
    _05_myplot._save(root='observed_speeds',
                     title='PMF of running speed',
                     xlabel='speed (mph)',
                     ylabel='probability')

    cdf = _13_Cdf._make_cdf_from_pmf(biased)

    _05_myplot._clf()
    _05_myplot._cdf(cdf)
    _05_myplot._save(root='observed_speeds_cdf',
                     title='CDF of running speed',
                     xlabel='speed (mph)',
                     ylabel='cumulative probability')


if __name__ == '__main__':
    main()
