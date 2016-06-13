"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import _05_myplot
import _10_relay
import _13_Cdf


def main():
    results = _10_relay._read_results()
    speeds = _10_relay._get_speeds(results)

    # plot the distribution of actual speeds
    cdf = _13_Cdf.MakeCdfFromList(speeds, 'speeds')

    _05_myplot._cdf(cdf)
    _05_myplot._save(root='relay_cdf',
                     title='CDF of running speed',
                     xlabel='speed (mph)',
                     ylabel='probability')


if __name__ == '__main__':
    main()
