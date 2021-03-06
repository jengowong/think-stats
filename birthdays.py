"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: birthdays.py
"""

import csv
import datetime
import sys
import _05_myplot
import _13_Cdf


def _read_birthdays(filename='birthdays.csv'):
    """
    Reads a CSV file of birthdays and returns a list of date objects.

    The first column of the file must contain dates in MM-DD format.

    Args:
        filename: string filename

    Returns:
        list of datetime.date objects
    """

    fp = open(filename)
    reader = csv.reader(fp)
    bdays = []

    for t in reader:
        bday = t[0]
        month, day = [int(x) for x in bday.split('-')]
        date = datetime.date(2010, month, day)
        bdays.append(date)

    return bdays


def _diff(t):
    """
    Computes the differences between the adjacent elements of a sequence.

    Args:
        t: sequence of anything subtractable

    Returns:
        list of whatever is in t
    """

    diffs = []
    for i in range(len(t) - 1):
        diff = t[i + 1] - t[i]
        diffs.append(diff)
    return diffs


def main(script):
    # read 'em and sort 'em
    birthdays = _read_birthdays()
    birthdays.sort()

    # compute the intervals in days
    deltas = _diff(birthdays)
    days = [inter.days for inter in deltas]

    # make and plot the CCDF on a log scale.
    cdf = _13_Cdf._make_cdf_from_list(days, name='intervals')
    scale = _05_myplot._cdf(cdf, transform='exponential')
    _05_myplot._save(root='intervals',
                     xlabel='days',
                     ylabel='ccdf',
                     **scale)


if __name__ == '__main__':
    main(*sys.argv)
