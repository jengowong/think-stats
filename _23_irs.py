"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import sys
import _04_Pmf
import _05_myplot
import _13_Cdf

"""
Results: on a log-log scale the tail of the CCDF is a straight line,
which suggests that the Pareto distribution is a good model for this data,
at least for people with taxable income above the median.
"""


def ReadIncomeFile(filename='08in11si.csv'):
    """
    Reads a data file from the IRS and returns the first two columns.

    Skips the header and returns only the first table (non-cumulative).

    Args:
        filename: string data file

    Returns:
        list of string pairs
    """
    reader = csv.reader(open(filename))
    for line in reader:
        if line[0] == 'All returns':
            break

    t = []
    for line in reader:
        if line[0].startswith('Accumulated'):
            break
        t.append(line[0:2])

    return t


def MakeIncomeDist(data):
    """
    Converts the strings from the IRS file to a Hist, Pmf and Cdf.

    Args:
        data: list of (dollar range, number of returns) string pairs

    Returns:
        tuple of (Hist, Pmf, Cdf) representing the number of returns in each bin
    """

    def clean(s):
        """Converts dollar amounts to integers."""
        try:
            return int(s.lstrip('$'))
        except ValueError:
            if s in ['No', 'income']:
                return 0
            if s == 'more':
                return -1
            return None

    def midpoint(low, high):
        """Finds the midpoint of a range."""
        if high == -1:
            return low * 3 / 2
        else:
            return (low + high) / 2

    hist = _04_Pmf.Hist()

    for column, number in data:
        # convert the number of returns
        number = number.replace(',', '')
        number = int(number)

        # convert the income range
        column = column.replace(',', '')
        t = column.split()
        low, high = t[0], t[-1]
        low, high = clean(low), clean(high)

        # add to the histogram
        x = midpoint(low, high)
        hist._incr(x, number)
        print(x, number)

    pmf = _04_Pmf._make_pmf_from_hist(hist)
    cdf = _13_Cdf.MakeCdfFromDict(pmf._get_dict())
    return hist, pmf, cdf


def main(script, *args):
    data = ReadIncomeFile()
    hist, pmf, cdf = MakeIncomeDist(data)

    # plot the CDF on a log-x scale
    _05_myplot._clf()
    _05_myplot._cdf(cdf)
    _05_myplot._save(root='income_logx',
                     xscale='log',
                     xlabel='income',
                     ylabel='CDF')

    # plot the complementary CDF on a log-log scale
    _05_myplot._clf()
    _05_myplot._cdf(cdf, complement=True)
    _05_myplot._save(root='income_loglog',
                     complement=True,
                     xscale='log',
                     yscale='log',
                     xlabel='income',
                     ylabel='complementary CDF')


if __name__ == "__main__":
    main(*sys.argv)
