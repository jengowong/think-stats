"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from scipy.special import erf, erfinv
import numpy
import math
import _04_Pmf
import _13_Cdf

root2 = math.sqrt(2.0)


def StandardNormalCdf(x):
    return (erf(x / root2) + 1) / 2


def NormalCdf(x, mu=0, sigma=1):
    """
    Evaluates the CDF of the normal distribution.
    
    Args:
        x:     float
        mu:    mean parameter
        sigma: standard deviation parameter
                
    Returns:
        float
    """
    return StandardNormalCdf(float(x - mu) / sigma)


def NormalCdfInverse(p, mu=0, sigma=1):
    """
    Evaluates the inverse CDF of the normal distribution.
    
    Args:
        p:     float
        mu:    mean parameter
        sigma: standard deviation parameter
                
    Returns:
        float
    """
    x = root2 * erfinv(2 * p - 1)
    return mu + x * sigma


spread = 4.0


def MakeNormalCdf(low=-spread, high=spread, digits=2):
    """
    Returns a Cdf object with the standard normal CDF.

    Args:
        low:    how many standard deviations below the mean?
        high:   how many standard deviations above the mean?
        digits:
    """
    n = (high - low) * 10 ** digits + 1
    xs = numpy.linspace(low, high, n)
    ps = (erf(xs / root2) + 1) / 2
    cdf = _13_Cdf.Cdf(xs, ps)
    return cdf


def MakeNormalPmf(low=-spread, high=spread, digits=2):
    """
    Returns a Pmf object with the standard normal CDF.

    Args:
        low:    how many standard deviations below the mean?
        high:   how many standard deviations above the mean?
        digits:
    """
    cdf = MakeNormalCdf(low, high, digits)
    pmf = _04_Pmf._make_pmf_from_cdf(cdf)
    return pmf


class FixedPointNormalPmf(_04_Pmf.Pmf):
    """
    A Pmf that maps from normal scores to probabilities.

    Values are rounded to the given number of digits.
    """

    def __init__(self, spread=4, digits=2, log=False):
        """
        Initalizes a FixedPointNormalPmf.

        Args:
            spread: how many standard deviations in each direction.
            digits: how many digits to round off to
            log:    whether to log-transform the probabilities
        """
        _04_Pmf.Pmf.__init__(self)
        self.spread = spread
        self.digits = digits

        n = 2 * spread * 10 ** digits + 1
        xs = numpy.linspace(-spread, spread, n)
        gap = (xs[1] - xs[0]) / 2

        for x in xs:
            p = StandardNormalCdf(x + gap) - StandardNormalCdf(x - gap)
            self._set(round(x, self.digits), p)

        # save the last (smallest) probability as the default for
        # values beyond the spread
        self.default = p

        self._normalize()
        if log:
            self._log()
            self.default = math.log(self.default)

    def NormalProb(self, x):
        """Looks up the probability for the value closest to x."""
        return self.d.get(round(x, self.digits), self.default)
