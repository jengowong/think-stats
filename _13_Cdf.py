"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2008 Allen B. Downey.
Distributed under the GNU General Public License at gnu.org/licenses/gpl.html.

NAME: _13_Cdf.py
"""

import bisect
import math
import random
import _04_Pmf

"""Functions for building CDFs (cumulative distribution functions)."""


class Cdf(object):
    """
    Represents a cumulative distribution function.

    Attributes:
        xs:   sequence of values
        ps:   sequence of probabilities
        name: string used as a graph label.
    """

    def __init__(self, xs=None, ps=None, name=''):
        self.xs = [] if xs is None else xs
        self.ps = [] if ps is None else ps
        self.name = name

    def _values(self):
        """Returns a sorted list of values."""
        return self.xs

    def _items(self):
        """
        Returns a sorted sequence of (value, probability) pairs.

        Note: in Python3, returns an iterator.
        """
        return zip(self.xs, self.ps)

    def _append(self, x, p):
        """
        Add an (x, p) pair to the end of this CDF.

        Note: this us normally used to build a CDF from scratch,
              not to modify existing CDFs.
              It is up to the caller to make sure that the result is a legal CDF.
        """
        self.xs.append(x)
        self.ps.append(p)

    def _prob(self, x):
        """
        Returns CDF(x), the probability that corresponds to value x.

        Args:
            x: number

        Returns:
            float probability
        """
        if x < self.xs[0]: return 0.0
        index = bisect.bisect(self.xs, x)
        p = self.ps[index - 1]
        return p

    def _value(self, p):
        """
        Returns InverseCDF(p), the value that corresponds to probability p.

        Args:
            p: number in the range [0, 1]

        Returns:
            number value
        """
        if p < 0 or p > 1:
            raise ValueError('Probability p must be in range [0, 1]')

        if p == 0: return self.xs[0]
        if p == 1: return self.xs[-1]
        index = bisect.bisect(self.ps, p)
        if p == self.ps[index - 1]:
            return self.xs[index - 1]
        else:
            return self.xs[index]

    def _percentile(self, p):
        """
        Returns the value that corresponds to percentile p.

        Args:
            p: number in the range [0, 100]

        Returns:
            number value
        """
        return self._value(p / 100.0)

    def _random(self):
        """Chooses a random value from this distribution."""
        return self._value(random.random())

    def _sample(self, n):
        """
        Generates a random sample from this distribution.
        
        Args:
            n: int length of the sample
        """
        return [self._random() for i in range(n)]

    def _mean(self):
        """
        Computes the mean of a CDF.

        Returns:
            float mean
        """
        old_p = 0
        total = 0.0
        for x, new_p in zip(self.xs, self.ps):
            p = new_p - old_p
            total += p * x
            old_p = new_p
        return total

    def _round(self, multiplier=1000.0):
        """
        An entry is added to the cdf only if the percentile differs
        from the previous value in a significant digit, where the number
        of significant digits is determined by multiplier.  The
        default is 1000, which keeps log10(1000) = 3 significant digits.
        """
        # TODO(write this method)
        pass

    def _render(self):
        """
        Generates a sequence of points suitable for plotting.

        An empirical CDF is a step function; linear interpolation can be misleading.

        Returns:
            tuple of (xs, ps)
        """
        xs = [self.xs[0]]
        ps = [0.0]
        for i, p in enumerate(self.ps):
            xs.append(self.xs[i])
            ps.append(p)

            try:
                xs.append(self.xs[i + 1])
                ps.append(p)
            except IndexError:
                pass
        return xs, ps


def _make_cdf_from_items(items, name=''):
    """
    Makes a cdf from an unsorted sequence of (value, frequency) pairs.

    Args:
        items: unsorted sequence of (value, frequency) pairs
        name:  string name for this CDF

    Returns:
        cdf: list of (value, fraction) pairs
    """
    runsum = 0
    xs = []
    cs = []

    for value, count in sorted(items):
        runsum += count
        xs.append(value)
        cs.append(runsum)

    total = float(runsum)
    ps = [c / total for c in cs]

    cdf = Cdf(xs, ps, name)
    return cdf


def _make_cdf_from_dict(d, name=''):
    """
    Makes a CDF from a dictionary that maps values to frequencies.

    Args:
       d:    dictionary that maps values to frequencies.
       name: string name for the data.

    Returns:
        Cdf object
    """
    return _make_cdf_from_items(d.iteritems(), name)


def _make_cdf_from_hist(hist, name=''):
    """
    Makes a CDF from a Hist object.

    Args:
       hist: Pmf.Hist object
       name: string name for the data.

    Returns:
        Cdf object
    """
    return _make_cdf_from_items(hist._items(), name)


def _make_cdf_from_pmf(pmf, name=None):
    """
    Makes a CDF from a Pmf object.

    Args:
       pmf:  Pmf.Pmf object
       name: string name for the data.

    Returns:
        Cdf object
    """
    if name == None:
        name = pmf.name
    return _make_cdf_from_items(pmf._items(), name)


def _make_cdf_from_list(seq, name=''):
    """
    Creates a CDF from an unsorted sequence.

    Args:
        seq:  unsorted sequence of sortable values
        name: string name for the cdf

    Returns:
       Cdf object
    """
    hist = _04_Pmf._make_hist_from_list(seq)
    return _make_cdf_from_hist(hist, name)
