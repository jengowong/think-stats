"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: _10_relay.py
"""

import urllib
import _04_Pmf
import _05_myplot

results = 'http://www.coolrunning.com/results/10/ma/Apr25_27thAn_set1.shtml'

"""
Sample line.

Place Div/Tot  Div   Guntime Nettime  Pace  Name                   Ag S Race# City/state              
===== ======== ===== ======= =======  ===== ====================== == = ===== ======================= 
    1   1/362  M2039   30:43   30:42   4:57 Brian Harvey           22 M  1422 Allston MA              
"""


def _convert_pace_to_speed(pace):
    """Converts pace in MM:SS per mile to MPH."""
    m, s = [int(x) for x in pace.split(':')]
    secs = m * 60 + s
    mph = 1.0 / secs * 60 * 60
    return mph


def _clean_line(line):
    """Converts a line from coolrunning results to a tuple of values."""
    t = line.split()
    if len(t) < 6:
        return None

    place, divtot, div, gun, net, pace = t[0:6]

    if not '/' in divtot:
        return None

    for time in [gun, net, pace]:
        if ':' not in time:
            return None

    return place, divtot, div, gun, net, pace


def _read_results(url=results):
    """Read results from coolrunning and return a list of tuples."""
    results = []
    conn = urllib.urlopen(url)
    for line in conn.fp:
        t = _clean_line(line)
        if t:
            results.append(t)
    return results


def _get_speeds(results, column=5):
    """Extract the pace column and return a list of speeds in MPH."""
    speeds = []
    for t in results:
        pace = t[column]
        speed = _convert_pace_to_speed(pace)
        speeds.append(speed)
    return speeds


def main():
    results = _read_results()
    speeds = _get_speeds(results)
    pmf = _04_Pmf._make_pmf_from_list(speeds, 'speeds')
    _05_myplot._pmf(pmf)
    _05_myplot._show(title='PMF of running speed',
                     xlabel='speed (mph)',
                     ylabel='probability')


if __name__ == '__main__':
    main()
