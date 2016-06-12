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
    results = _10_relay.ReadResults()
    speeds = _10_relay.GetSpeeds(results)

    # plot the distribution of actual speeds
    cdf = _13_Cdf.MakeCdfFromList(speeds, 'speeds')

    _05_myplot.Cdf(cdf)
    _05_myplot.Save(root='relay_cdf',
                    title='CDF of running speed',
                    xlabel='speed (mph)',
                    ylabel='probability')


if __name__ == '__main__':
    main()