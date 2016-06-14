"""
This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: thinkstats2.py
"""

from __future__ import print_function

"""
This file contains class definitions for:
    Hist:         represents a histogram (map from values to integer frequencies).
    Pmf:          represents a probability mass function (map from values to probs).
    _DictWrapper: private parent class for Hist and Pmf.
    Cdf:          represents a discrete cumulative distribution function
    Pdf:          represents a continuous probability density function
"""

import bisect
import copy
import logging
import math
import random
import re
import numpy as np

import thinkplot
import pandas as pd
import scipy.stats
from scipy.special import erf, erfinv

ROOT2 = math.sqrt(2)


def _random_seed(x):
    """
    Initialize the random and np.random generators.

    Args:
        x: int seed
    """
    random.seed(x)
    np.random.seed(x)


def _odds(p):
    """
    Computes odds for a given probability.

    Example: p=0.75 means 75 for and 25 against, or 3:1 odds in favor.

    Note: when p=1, the formula for odds divides by zero, which is
    normally undefined.  But I think it is reasonable to define Odds(1)
    to be infinity, so that's what this function does.

    Args:
        p: float 0-1

    Returns:
        float odds
    """
    if p == 1:
        return float('inf')
    return p / (1 - p)


def _probability(o):
    """
    Computes the probability corresponding to given odds.

    Example: o=2 means 2:1 odds in favor, or 2/3 probability

    Args:
        o: float odds, strictly positive

    Returns:
        float probability
    """
    return o / (o + 1)


def _probability2(yes, no):
    """
    Computes the probability corresponding to given odds.

    Example: yes=2, no=1 means 2:1 odds in favor, or 2/3 probability.

    Args:
        yes, no: int or float odds in favor
    """
    return float(yes) / (yes + no)


class Interpolator(object):
    """
    Represents a mapping between sorted sequences; performs linear interp.

    Attributes:
        xs: sorted list
        ys: sorted list
    """

    def __init__(self, xs, ys):
        self.xs = xs
        self.ys = ys

    def _lookup(self, x):
        """Looks up x and returns the corresponding value of y."""
        return self._bisect(x, self.xs, self.ys)

    def _reverse(self, y):
        """Looks up y and returns the corresponding value of x."""
        return self._bisect(y, self.ys, self.xs)

    @staticmethod
    def _bisect(x, xs, ys):
        """Helper function."""
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        i = bisect.bisect(xs, x)
        frac = 1.0 * (x - xs[i - 1]) / (xs[i] - xs[i - 1])
        y = ys[i - 1] + frac * 1.0 * (ys[i] - ys[i - 1])
        return y


class _DictWrapper(object):
    """An object that contains a dictionary."""

    def __init__(self, values=None, name=''):
        """
        Initializes the distribution.

        hypos: sequence of hypotheses
        """
        self.name = name
        self.d = {}

        # flag whether the distribution is under a log transform
        self.log = False

        if values is None:
            return

        init_methods = [
            self._init_pmf,
            self._init_mapping,
            self._init_sequence,
            self._init_failure,
        ]

        for method in init_methods:
            try:
                method(values)
                break
            except AttributeError:
                continue

        if len(self) > 0 and isinstance(self, Pmf):
            self._normalize()

    def _init_sequence(self, values):
        """
        Initializes with a sequence of equally-likely values.

        values: sequence of values
        """
        for value in values:
            self._set(value, 1)

    def _init_mapping(self, values):
        """
        Initializes with a map from value to probability.

        values: map from value to probability
        """
        for value, prob in values.items():
            self._set(value, prob)

    def _init_pmf(self, values):
        """
        Initializes with a Pmf.

        values: Pmf object
        """
        for value, prob in values._items():
            self._set(value, prob)

    @staticmethod
    def _init_failure(values):
        """Raises an error."""
        raise ValueError('None of the initialization methods worked.')

    def __len__(self):
        return len(self.d)

    def __iter__(self):
        return iter(self.d)

    def iterkeys(self):
        return iter(self.d)

    def __contains__(self, value):
        return value in self.d

    def _copy(self, name=None):
        """
        Returns a copy.

        Make a shallow copy of d.  If you want a deep copy of d,
        use copy.deepcopy on the whole object.

        Args:
            name: string name for the new Hist
        """
        new = copy.copy(self)
        new.d = copy.copy(self.d)
        new.name = name if name is not None else self.name
        return new

    def _scale(self, factor):
        """
        Multiplies the values by a factor.

        Args:
            factor: what to multiply by

        Returns:
            new object
        """
        new = self._copy()
        new.d.clear()

        for val, prob in self._items():
            new._set(val * factor, prob)
        return new

    def _log(self, m=None):
        """
        Log transforms the probabilities.
        
        Removes values with probability 0.

        Normalizes so that the largest logprob is 0.
        """
        if self.log:
            raise ValueError("Pmf/Hist already under a log transform")
        self.log = True

        if m is None:
            m = self._max_like()

        for x, p in self.d.items():
            if p:
                self._set(x, math.log(p / m))
            else:
                self._remove(x)

    def _exp(self, m=None):
        """
        Exponentiates the probabilities.

        Args:
            m: how much to shift the ps before exponentiating

        If m is None, normalizes so that the largest prob is 1.
        """
        if not self.log:
            raise ValueError("Pmf/Hist not under a log transform")
        self.log = False

        if m is None:
            m = self._max_like()

        for x, p in self.d.items():
            self._set(x, math.exp(p - m))

    def _get_dict(self):
        """Gets the dictionary."""
        return self.d

    def _set_dict(self, d):
        """Sets the dictionary."""
        self.d = d

    def _values(self):
        """
        Gets an unsorted sequence of values.

        Note: one source of confusion is that the keys of this
        dictionary are the values of the Hist/Pmf, and the
        values of the dictionary are frequencies/probabilities.
        """
        return self.d.keys()

    def _items(self):
        """Gets an unsorted sequence of (value, freq/prob) pairs."""
        return self.d.items()

    def _render(self):
        """
        Generates a sequence of points suitable for plotting.

        Returns:
            tuple of (sorted value sequence, freq/prob sequence)
        """
        return zip(*sorted(self._items()))

    def _print(self):
        """Prints the values and freqs/probs in ascending order."""
        for val, prob in sorted(self.d.items()):
            print(val, prob)

    def _set(self, x, y=0):
        """
        Sets the freq/prob associated with the value x.

        Args:
            x: number value
            y: number freq or prob
        """
        self.d[x] = y

    def _incr(self, x, term=1):
        """
        Increments the freq/prob associated with the value x.

        Args:
            x: number value
            term: how much to increment by
        """
        self.d[x] = self.d.get(x, 0) + term

    def _mult(self, x, factor):
        """
        Scales the freq/prob associated with the value x.

        Args:
            x: number value
            factor: how much to multiply by
        """
        self.d[x] = self.d.get(x, 0) * factor

    def _remove(self, x):
        """
        Removes a value.

        Throws an exception if the value is not there.

        Args:
            x: value to remove
        """
        del self.d[x]

    def _total(self):
        """Returns the total of the frequencies/probabilities in the map."""
        total = sum(self.d.itervalues())
        return total

    def _max_like(self):
        """Returns the largest frequency/probability in the map."""
        return max(self.d.itervalues())


class Hist(_DictWrapper):
    """
    Represents a histogram, which is a map from values to frequencies.

    Values can be any hashable type; frequencies are integer counters.
    """

    def _freq(self, x):
        """
        Gets the frequency associated with the value x.

        Args:
            x: number value

        Returns:
            int frequency
        """
        return self.d.get(x, 0)

    def _freqs(self, xs):
        """Gets frequencies for a sequence of values."""
        return [self._freq(x) for x in xs]

    def _is_subset(self, other):
        """Checks whether the values in this histogram are a subset of the values in the given histogram."""
        for val, freq in self._items():
            if freq > other._freq(val):
                return False
        return True

    def _subtract(self, other):
        """Subtracts the values in the given histogram from this histogram."""
        for val, freq in other._items():
            self._incr(val, -freq)


class Pmf(_DictWrapper):
    """
    Represents a probability mass function.
    
    Values can be any hashable type; probabilities are floating-point.
    Pmfs are not necessarily normalized.
    """

    def _prob(self, x, default=0):
        """
        Gets the probability associated with the value x.

        Args:
            x: number value
            default: value to return if the key is not there

        Returns:
            float probability
        """
        return self.d.get(x, default)

    def _probs(self, xs):
        """Gets probabilities for a sequence of values."""
        return [self._prob(x) for x in xs]

    def _make_cdf(self, name=None):
        """Makes a Cdf."""
        return _make_cdf_from_pmf(self, name=name)

    def _prob_greater(self, x):
        """
        Probability that a sample from this Pmf exceeds x.

        Args:
            x: number

        returns:
            float probability
        """
        t = [prob for (val, prob) in self.d.items() if val > x]
        return sum(t)

    def _prob_less(self, x):
        """
        Probability that a sample from this Pmf is less than x.

        Args:
            x: number

        returns:
            float probability
        """
        t = [prob for (val, prob) in self.d.items() if val < x]
        return sum(t)

    def __lt__(self, obj):
        """
        Less than.

        Args:
            obj: number or _DictWrapper

        returns:
            float probability
        """
        if isinstance(obj, _DictWrapper):
            return _pmf_prob_less(self, obj)
        else:
            return self._prob_less(obj)

    def __gt__(self, obj):
        """
        Greater than.

        Args:
            obj: number or _DictWrapper

        returns:
            float probability
        """
        if isinstance(obj, _DictWrapper):
            return _pmf_prob_greater(self, obj)
        else:
            return self._prob_greater(obj)

    def __ge__(self, obj):
        """
        Greater than or equal.

        Args:
            obj: number or _DictWrapper

        returns:
            float probability
        """
        return 1 - (self < obj)

    def __le__(self, obj):
        """
        Less than or equal.

        Args:
            obj: number or _DictWrapper

        returns:
            float probability
        """
        return 1 - (self > obj)

    def __eq__(self, obj):
        """
        Less than.

        Args:
            obj: number or _DictWrapper

        returns:
            float probability
        """
        if isinstance(obj, _DictWrapper):
            return _pmf_prob_equal(self, obj)
        else:
            return self._prob(obj)

    def __ne__(self, obj):
        """
        Less than.

        Args:
            obj: number or _DictWrapper

        returns:
            float probability
        """
        return 1 - (self == obj)

    def _normalize(self, fraction=1.0):
        """
        Normalizes this PMF so the sum of all probs is fraction.

        Args:
            fraction: what the total should be after normalization

        Returns:
            the total probability before normalizing
        """
        if self.log:
            raise ValueError("Pmf is under a log transform")

        total = self._total()
        if total == 0.0:
            raise ValueError('total probability is zero.')
            logging.warning('Normalize: total probability is zero.')
            return total

        factor = float(fraction) / total
        for x in self.d:
            self.d[x] *= factor

        return total

    def _random(self):
        """
        Chooses a random element from this PMF.

        Returns:
            float value from the Pmf
        """
        if len(self.d) == 0:
            raise ValueError('Pmf contains no values.')

        target = random.random()
        total = 0.0
        for x, p in self.d.items():
            total += p
            if total >= target:
                return x

        # we shouldn't get here
        assert False

    def _mean(self):
        """
        Computes the mean of a PMF.

        Returns:
            float mean
        """
        mean = 0.0
        for x, p in self.d.items():
            mean += p * x
        return mean

    def _var(self, mu=None):
        """
        Computes the variance of a PMF.

        Args:
            mu: the point around which the variance is computed; if omitted, computes the mean

        Returns:
            float variance
        """
        if mu is None:
            mu = self._mean()

        var = 0.0
        for x, p in self.d.items():
            var += p * (x - mu) ** 2
        return var

    def _maximum_likelihood(self):
        """
        Returns the value with the highest probability.

        Returns:
            float probability
        """
        prob, val = max((prob, val) for val, prob in self._items())
        return val

    def _credible_interval(self, percentage=90):
        """
        Computes the central credible interval.

        If percentage=90, computes the 90% CI.

        Args:
            percentage: float between 0 and 100

        Returns:
            sequence of two floats, low and high
        """
        cdf = self._make_cdf()
        return cdf._credible_interval(percentage)

    def __add__(self, other):
        """
        Computes the Pmf of the sum of values drawn from self and other.

        Args:
            other: another Pmf

        returns:
            new Pmf
        """
        try:
            return self._add_pmf(other)
        except AttributeError:
            return self._add_constant(other)

    def _add_pmf(self, other):
        """
        Computes the Pmf of the sum of values drawn from self and other.

        Args:
            other: another Pmf

        returns:
            new Pmf
        """
        pmf = Pmf()
        for v1, p1 in self._items():
            for v2, p2 in other._items():
                pmf._incr(v1 + v2, p1 * p2)
        return pmf

    def _add_constant(self, other):
        """
        Computes the Pmf of the sum a constant and  values from self.

        Args:
            other: a number

        returns:
            new Pmf
        """
        pmf = Pmf()
        for v1, p1 in self._items():
            pmf._set(v1 + other, p1)
        return pmf

    def __sub__(self, other):
        """
        Computes the Pmf of the diff of values drawn from self and other.

        Args:
            other: another Pmf

        returns:
            new Pmf
        """
        pmf = Pmf()
        for v1, p1 in self._items():
            for v2, p2 in other._items():
                pmf._incr(v1 - v2, p1 * p2)
        return pmf

    def _max(self, k):
        """
        Computes the CDF of the maximum of k selections from this dist.

        Args:
            k: int

        returns:
            new Cdf
        """
        cdf = self._make_cdf()
        cdf.ps = [p ** k for p in cdf.ps]
        return cdf


class Joint(Pmf):
    """
    Represents a joint distribution.

    The values are sequences (usually tuples)
    """

    def _marginal(self, i, name=''):
        """
        Gets the marginal distribution of the indicated variable.

        Args:
            i: index of the variable we want

        Returns:
            Pmf
        """
        pmf = Pmf(name=name)
        for vs, prob in self._items():
            pmf._incr(vs[i], prob)
        return pmf

    def _conditional(self, i, j, val, name=''):
        """
        Gets the conditional distribution of the indicated variable.

        Distribution of vs[i], conditioned on vs[j] = val.

        Args:
            i:   index of the variable we want
            j:   which variable is conditioned on
            val: the value the jth variable has to have

        Returns:
            Pmf
        """
        pmf = Pmf(name=name)
        for vs, prob in self._items():
            if vs[j] != val: continue
            pmf._incr(vs[i], prob)

        pmf._normalize()
        return pmf

    def _max_like_interval(self, percentage=90):
        """
        Returns the maximum-likelihood credible interval.

        If percentage=90, computes a 90% CI containing the values with the highest likelihoods.

        Args:
            percentage: float between 0 and 100

        Returns:
            list of values from the suite
        """
        interval = []
        total = 0

        t = [(prob, val) for val, prob in self._items()]
        t.sort(reverse=True)

        for prob, val in t:
            interval.append(val)
            total += prob
            if total >= percentage / 100.0:
                break

        return interval


def _make_joint(pmf1, pmf2):
    """
    Joint distribution of values from pmf1 and pmf2.

    Args:
        pmf1: Pmf object
        pmf2: Pmf object

    Returns:
        Joint pmf of value pairs
    """
    joint = Joint()
    for v1, p1 in pmf1._items():
        for v2, p2 in pmf2._items():
            joint._set((v1, v2), p1 * p2)
    return joint


def _make_hist_from_list(t, name=''):
    """
    Makes a histogram from an unsorted sequence of values.

    Args:
        t:    sequence of numbers
        name: string name for this histogram

    Returns:
        Hist object
    """
    hist = Hist(name=name)
    [hist._incr(x) for x in t]
    return hist


def _make_hist_from_dict(d, name=''):
    """
    Makes a histogram from a map from values to frequencies.

    Args:
        d:    dictionary that maps values to frequencies
        name: string name for this histogram

    Returns:
        Hist object
    """
    return Hist(d, name)


def _make_pmf_from_list(t, name=''):
    """
    Makes a PMF from an unsorted sequence of values.

    Args:
        t:    sequence of numbers
        name: string name for this PMF

    Returns:
        Pmf object
    """
    hist = _make_hist_from_list(t)
    d = hist._get_dict()
    pmf = Pmf(d, name)
    pmf._normalize()
    return pmf


def _make_pmf_from_dict(d, name=''):
    """
    Makes a PMF from a map from values to probabilities.

    Args:
        d:    dictionary that maps values to probabilities
        name: string name for this PMF

    Returns:
        Pmf object
    """
    pmf = Pmf(d, name)
    pmf._normalize()
    return pmf


def _make_pmf_from_items(t, name=''):
    """
    Makes a PMF from a sequence of value-probability pairs

    Args:
        t:    sequence of value-probability pairs
        name: string name for this PMF

    Returns:
        Pmf object
    """
    pmf = Pmf(dict(t), name)
    pmf._normalize()
    return pmf


def _make_pmf_from_hist(hist, name=None):
    """
    Makes a normalized PMF from a Hist object.

    Args:
        hist: Hist object
        name: string name

    Returns:
        Pmf object
    """
    if name is None:
        name = hist.name

    # make a copy of the dictionary
    d = dict(hist._get_dict())
    pmf = Pmf(d, name)
    pmf._normalize()
    return pmf


def _make_pmf_from_cdf(cdf, name=None):
    """
    Makes a normalized Pmf from a Cdf object.

    Args:
        cdf:  Cdf object
        name: string name for the new Pmf

    Returns:
        Pmf object
    """
    if name is None:
        name = cdf.name

    pmf = Pmf(name=name)

    prev = 0.0
    for val, prob in cdf._items():
        pmf._incr(val, prob - prev)
        prev = prob

    return pmf


def _make_mixture(metapmf, name='mix'):
    """
    Make a mixture distribution.

    Args:
        metapmf: Pmf that maps from Pmfs to probs.
        name:    string name for the new Pmf.

    Returns: Pmf object.
    """
    mix = Pmf(name=name)
    for pmf, p1 in metapmf._items():
        for x, p2 in pmf._items():
            mix._incr(x, p1 * p2)
    return mix


def _make_uniform_pmf(low, high, n):
    """
    Make a uniform Pmf.

    Args:
        low:  lowest value (inclusive)
        high: highest value (inclusize)
        n:    number of values
    """
    pmf = Pmf()
    for x in np.linspace(low, high, n):
        pmf._set(x, 1)
    pmf._normalize()
    return pmf


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

    def _copy(self, name=None):
        """
        Returns a copy of this Cdf.

        Args:
            name: string name for the new Cdf
        """
        if name is None:
            name = self.name
        return Cdf(list(self.xs), list(self.ps), name)

    def _make_pmf(self, name=None):
        """Makes a Pmf."""
        return _make_pmf_from_cdf(self, name=name)

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

        Note: this us normally used to build a CDF from scratch, not
        to modify existing CDFs.  It is up to the caller to make sure
        that the result is a legal CDF.
        """
        self.xs.append(x)
        self.ps.append(p)

    def _shift(self, term):
        """
        Adds a term to the xs.

        Args:
            term: how much to add
        """
        new = self._copy()
        new.xs = [x + term for x in self.xs]
        return new

    def _scale(self, factor):
        """
        Multiplies the xs by a factor.

        Args:
            factor: what to multiply by
        """
        new = self._copy()
        new.xs = [x * factor for x in self.xs]
        return new

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

    def _credible_interval(self, percentage=90):
        """
        Computes the central credible interval.

        If percentage=90, computes the 90% CI.

        Args:
            percentage: float between 0 and 100

        Returns:
            sequence of two floats, low and high
        """
        prob = (1 - percentage / 100.0) / 2
        interval = self._value(prob), self._value(1 - prob)
        return interval

    @staticmethod
    def _round(multiplier=1000.0):
        """
        An entry is added to the cdf only if the percentile differs
        from the previous value in a significant digit, where the number
        of significant digits is determined by multiplier.  The
        default is 1000, which keeps log10(1000) = 3 significant digits.
        """
        # TODO(write this method)
        raise UnimplementedMethodException()

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

    def _max(self, k):
        """
        Computes the CDF of the maximum of k selections from this dist.

        Args:
            k: int

        returns:
            new Cdf
        """
        cdf = self._copy()
        cdf.ps = [p ** k for p in cdf.ps]
        return cdf


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
    return _make_cdf_from_items(d.items(), name)


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
    hist = _make_hist_from_list(seq)
    return _make_cdf_from_hist(hist, name)


class UnimplementedMethodException(Exception):
    """Exception if someone calls a method that should be overridden."""


class Suite(Pmf):
    """Represents a suite of hypotheses and their probabilities."""

    def _update(self, data):
        """
        Updates each hypothesis based on the data.

        Args:
            data: any representation of the data

        returns:
            the normalizing constant
        """
        for hypo in self._values():
            like = self._likelihood(data, hypo)
            self._mult(hypo, like)
        return self._normalize()

    def _log_update(self, data):
        """
        Updates a suite of hypotheses based on new data.

        Modifies the suite directly; if you want to keep the original, make a copy.

        Note: unlike Update, LogUpdate does not normalize.

        Args:
            data: any representation of the data
        """
        for hypo in self._values():
            like = self._log_likelihood(data, hypo)
            self._incr(hypo, like)

    def _update_set(self, dataset):
        """
        Updates each hypothesis based on the dataset.

        This is more efficient than calling Update repeatedly because
        it waits until the end to Normalize.

        Modifies the suite directly; if you want to keep the original, make a copy.

        Args:
            dataset: a sequence of data

        returns:
            the normalizing constant
        """
        for data in dataset:
            for hypo in self._values():
                like = self._likelihood(data, hypo)
                self._mult(hypo, like)
        return self._normalize()

    def _log_update_set(self, dataset):
        """
        Updates each hypothesis based on the dataset.

        Modifies the suite directly; if you want to keep the original, make a copy.

        Args:
            dataset: a sequence of data

        returns:
            None
        """
        for data in dataset:
            self._log_update(data)

    @staticmethod
    def _likelihood(data, hypo):
        """
        Computes the likelihood of the data under the hypothesis.

        Args:
            hypo: some representation of the hypothesis
            data: some representation of the data
        """
        raise UnimplementedMethodException()

    @staticmethod
    def _log_likelihood(data, hypo):
        """
        Computes the log likelihood of the data under the hypothesis.

        Args:
            hypo: some representation of the hypothesis
            data: some representation of the data
        """
        raise UnimplementedMethodException()

    def _print(self):
        """Prints the hypotheses and their probabilities."""
        for hypo, prob in sorted(self._items()):
            print(hypo, prob)

    def _make_odds(self):
        """
        Transforms from probabilities to odds.

        Values with prob=0 are removed.
        """
        for hypo, prob in self._items():
            if prob:
                self._set(hypo, _odds(prob))
            else:
                self._remove(hypo)

    def _make_probs(self):
        """Transforms from odds to probabilities."""
        for hypo, odds in self._items():
            self._set(hypo, _probability(odds))


def _make_suite_from_list(t, name=''):
    """
    Makes a suite from an unsorted sequence of values.

    Args:
        t:    sequence of numbers
        name: string name for this suite

    Returns:
        Suite object
    """
    hist = _make_hist_from_list(t)
    d = hist._get_dict()
    return _make_suite_from_dict(d)


def _make_suite_from_hist(hist, name=None):
    """
    Makes a normalized suite from a Hist object.

    Args:
        hist: Hist object
        name: string name

    Returns:
        Suite object
    """
    if name is None:
        name = hist.name

    # make a copy of the dictionary
    d = dict(hist._get_dict())
    return _make_suite_from_dict(d, name)


def _make_suite_from_dict(d, name=''):
    """
    Makes a suite from a map from values to probabilities.

    Args:
        d:    dictionary that maps values to probabilities
        name: string name for this suite

    Returns:
        Suite object
    """
    suite = Suite(name=name)
    suite._set_dict(d)
    suite._normalize()
    return suite


def _make_suite_from_cdf(cdf, name=None):
    """
    Makes a normalized Suite from a Cdf object.

    Args:
        cdf:  Cdf object
        name: string name for the new Suite

    Returns:
        Suite object
    """
    if name is None:
        name = cdf.name

    suite = Suite(name=name)

    prev = 0.0
    for val, prob in cdf._items():
        suite._incr(val, prob - prev)
        prev = prob

    return suite


class Pdf(object):
    """Represents a probability density function (PDF)."""

    def _density(self, x):
        """
        Evaluates this Pdf at x.

        Returns:
            float probability density
        """
        raise UnimplementedMethodException()

    def _make_pmf(self, xs, name=''):
        """
        Makes a discrete version of this Pdf, evaluated at xs.

        Args:
            xs: equally-spaced sequence of values

        Returns: new Pmf
        """
        pmf = Pmf(name=name)
        for x in xs:
            pmf._set(x, self._density(x))
        pmf._normalize()
        return pmf


class GaussianPdf(Pdf):
    """Represents the PDF of a Gaussian distribution."""

    def __init__(self, mu, sigma):
        """
        Constructs a Gaussian Pdf with given mu and sigma.

        Args:
            mu:    mean
            sigma: standard deviation
        """
        self.mu = mu
        self.sigma = sigma

    def _density(self, x):
        """
        Evaluates this Pdf at x.

        Returns:
            float probability density
        """
        return _eval_gaussian_pdf(x, self.mu, self.sigma)


class EstimatedPdf(Pdf):
    """Represents a PDF estimated by KDE."""

    def __init__(self, sample):
        """
        Estimates the density function based on a sample.

        Args:
            sample: sequence of data
        """
        self.kde = scipy.stats.gaussian_kde(sample)

    def _density(self, x):
        """
        Evaluates this Pdf at x.

        Returns:
            float probability density
        """
        return self.kde.evaluate(x)

    def _make_pmf(self, xs, name=''):
        ps = self.kde.evaluate(xs)
        pmf = _make_pmf_from_items(zip(xs, ps), name=name)
        return pmf


def _percentile(pmf, percentage):
    """
    Computes a percentile of a given Pmf.

    Args:
        percentage: float 0-100
    """
    p = percentage / 100.0
    total = 0
    for val, prob in pmf._items():
        total += prob
        if total >= p:
            return val


def _credible_interval(pmf, percentage=90):
    """
    Computes a credible interval for a given distribution.

    If percentage=90, computes the 90% CI.

    Args:
        pmf:        Pmf object representing a posterior distribution
        percentage: float between 0 and 100

    Returns:
        sequence of two floats, low and high
    """
    cdf = pmf._make_cdf()
    prob = (1 - percentage / 100.0) / 2
    interval = cdf._value(prob), cdf._value(1 - prob)
    return interval


def _pmf_prob_less(pmf1, pmf2):
    """
    Probability that a value from pmf1 is less than a value from pmf2.

    Args:
        pmf1: Pmf object
        pmf2: Pmf object

    Returns:
        float probability
    """
    total = 0.0
    for v1, p1 in pmf1._items():
        for v2, p2 in pmf2._items():
            if v1 < v2:
                total += p1 * p2
    return total


def _pmf_prob_greater(pmf1, pmf2):
    """
    Probability that a value from pmf1 is less than a value from pmf2.

    Args:
        pmf1: Pmf object
        pmf2: Pmf object

    Returns:
        float probability
    """
    total = 0.0
    for v1, p1 in pmf1._items():
        for v2, p2 in pmf2._items():
            if v1 > v2:
                total += p1 * p2
    return total


def _pmf_prob_equal(pmf1, pmf2):
    """
    Probability that a value from pmf1 equals a value from pmf2.

    Args:
        pmf1: Pmf object
        pmf2: Pmf object

    Returns:
        float probability
    """
    total = 0.0
    for v1, p1 in pmf1._items():
        for v2, p2 in pmf2._items():
            if v1 == v2:
                total += p1 * p2
    return total


def _random_sum(dists):
    """
    Chooses a random value from each dist and returns the sum.

    Args:
        dists: sequence of Pmf or Cdf objects

    returns:
        numerical sum
    """
    total = sum(dist._random() for dist in dists)
    return total


def _sample_sum(dists, n):
    """
    Draws a sample of sums from a list of distributions.

    Args:
        dists: sequence of Pmf or Cdf objects
        n:     sample size

    returns:
        new Pmf of sums
    """
    pmf = _make_pmf_from_list(_random_sum(dists) for i in xrange(n))
    return pmf


def _eval_gaussian_pdf(x, mu, sigma):
    """
    Computes the unnormalized PDF of the normal distribution.

    Args:
        x:     value
        mu:    mean
        sigma: standard deviation
    
    returns:
        float probability density
    """
    return scipy.stats.norm.pdf(x, mu, sigma)


def _make_gaussian_pmf(mu, sigma, num_sigmas, n=201):
    """
    Makes a PMF discrete approx to a Gaussian distribution.

    Args:
        mu:         float mean
        sigma:      float standard deviation
        num_sigmas: how many sigmas to extend in each direction
        n:          number of values in the Pmf

    returns:
        normalized Pmf
    """
    pmf = Pmf()
    low = mu - num_sigmas * sigma
    high = mu + num_sigmas * sigma

    for x in np.linspace(low, high, n):
        p = _eval_gaussian_pdf(x, mu, sigma)
        pmf._set(x, p)
    pmf._normalize()
    return pmf


def _eval_binomial_pmf(k, n, p):
    """
    Evaluates the binomial pmf.

    Returns the probabily of k successes in n trials with probability p.
    """
    return scipy.stats.binom.pmf(k, n, p)


def _eval_poisson_pmf(k, lam):
    """
    Computes the Poisson PMF.

    Args:
        k:   number of events
        lam: parameter lambda in events per unit time

    returns:
        float probability
    """
    # don't use the scipy function (yet).  for lam=0 it returns NaN;
    # should be 0.0
    # return scipy.stats.poisson.pmf(k, lam)

    return lam ** k * math.exp(-lam) / math.factorial(k)


def _make_poisson_pmf(lam, high, step=1):
    """
    Makes a PMF discrete approx to a Poisson distribution.

    Args:
        lam:  parameter lambda in events per unit time
        high: upper bound of the Pmf

    returns:
        normalized Pmf
    """
    pmf = Pmf()
    for k in xrange(0, high + 1, step):
        p = _eval_poisson_pmf(k, lam)
        pmf._set(k, p)

    pmf._normalize()
    return pmf


def _eval_exponential_pdf(x, lam):
    """
    Computes the exponential PDF.

    Args:
        x:   value
        lam: parameter lambda in events per unit time

    returns:
        float probability density
    """
    return lam * math.exp(-lam * x)


def _eval_exponential_cdf(x, lam):
    """Evaluates CDF of the exponential distribution with parameter lam."""
    return 1 - math.exp(-lam * x)


def _make_exponential_pmf(lam, high, n=200):
    """
    Makes a PMF discrete approx to an exponential distribution.

    Args:
        lam:  parameter lambda in events per unit time
        high: upper bound
        n:    number of values in the Pmf

    returns:
        normalized Pmf
    """
    pmf = Pmf()
    for x in np.linspace(0, high, n):
        p = _eval_exponential_pdf(x, lam)
        pmf._set(x, p)
    pmf._normalize()
    return pmf


def _standard_gaussian_cdf(x):
    """
    Evaluates the CDF of the standard Gaussian distribution.
    
    See http://en.wikipedia.org/wiki/Normal_distribution#Cumulative_distribution_function

    Args:
        x: float
                
    Returns:
        float
    """
    return (erf(x / ROOT2) + 1) / 2


def _gaussian_cdf(x, mu=0, sigma=1):
    """
    Evaluates the CDF of the gaussian distribution.
    
    Args:
        x:     float
        mu:    mean parameter
        sigma: standard deviation parameter
                
    Returns:
        float
    """
    return _standard_gaussian_cdf(float(x - mu) / sigma)


def _gaussian_cdf_inverse(p, mu=0, sigma=1):
    """
    Evaluates the inverse CDF of the gaussian distribution.

    See http://en.wikipedia.org/wiki/Normal_distribution#Quantile_function  

    Args:
        p:     float
        mu:    mean parameter
        sigma: standard deviation parameter
                
    Returns:
        float
    """
    x = ROOT2 * erfinv(2 * p - 1)
    return mu + x * sigma


class Beta(object):
    """
    Represents a Beta distribution.

    See http://en.wikipedia.org/wiki/Beta_distribution
    """

    def __init__(self, alpha=1, beta=1, name=''):
        """Initializes a Beta distribution."""
        self.alpha = alpha
        self.beta = beta
        self.name = name

    def _update(self, data):
        """
        Updates a Beta distribution.

        Args:
            data: pair of int (heads, tails)
        """
        heads, tails = data
        self.alpha += heads
        self.beta += tails

    def _mean(self):
        """Computes the mean of this distribution."""
        return float(self.alpha) / (self.alpha + self.beta)

    def _random(self):
        """Generates a random variate from this distribution."""
        return random.betavariate(self.alpha, self.beta)

    def _sample(self, n):
        """
        Generates a random sample from this distribution.

        Args:
            n: int sample size
        """
        size = n,
        return np.random.beta(self.alpha, self.beta, size)

    def _eval_pdf(self, x):
        """Evaluates the PDF at x."""
        return x ** (self.alpha - 1) * (1 - x) ** (self.beta - 1)

    def _make_pmf(self, steps=101, name=''):
        """
        Returns a Pmf of this distribution.

        Note: Normally, we just evaluate the PDF at a sequence
        of points and treat the probability density as a probability
        mass.

        But if alpha or beta is less than one, we have to be
        more careful because the PDF goes to infinity at x=0
        and x=1.  In that case we evaluate the CDF and compute
        differences.
        """
        if self.alpha < 1 or self.beta < 1:
            cdf = self._make_cdf()
            pmf = cdf._make_pmf()
            return pmf

        xs = [i / (steps - 1.0) for i in xrange(steps)]
        probs = [self._eval_pdf(x) for x in xs]
        pmf = _make_pmf_from_dict(dict(zip(xs, probs)), name)
        return pmf

    def _make_cdf(self, steps=101):
        """Returns the CDF of this distribution."""
        xs = [i / (steps - 1.0) for i in xrange(steps)]
        ps = [scipy.special.betainc(self.alpha, self.beta, x) for x in xs]
        cdf = Cdf(xs, ps)
        return cdf


class Dirichlet(object):
    """
    Represents a Dirichlet distribution.

    See http://en.wikipedia.org/wiki/Dirichlet_distribution
    """

    def __init__(self, n, conc=1, name=''):
        """
        Initializes a Dirichlet distribution.

        Args:
            n:    number of dimensions
            conc: concentration parameter (smaller yields more concentration)
            name: string name
        """
        if n < 2:
            raise ValueError('A Dirichlet distribution with '
                             'n<2 makes no sense')

        self.n = n
        self.params = np.ones(n, dtype=np.float) * conc
        self.name = name

    def _update(self, data):
        """
        Updates a Dirichlet distribution.

        Args:
            data: sequence of observations, in order corresponding to params
        """
        m = len(data)
        self.params[:m] += data

    def _random(self):
        """
        Generates a random variate from this distribution.

        Returns:
            normalized vector of fractions
        """
        p = np.random.gamma(self.params)
        return p / p.sum()

    def _likelihood(self, data):
        """
        Computes the likelihood of the data.

        Selects a random vector of probabilities from this distribution.

        Returns:
            float probability
        """
        m = len(data)
        if self.n < m:
            return 0

        x = data
        p = self._random()
        q = p[:m] ** x
        return q.prod()

    def _log_likelihood(self, data):
        """
        Computes the log likelihood of the data.

        Selects a random vector of probabilities from this distribution.

        Returns:
            float log probability
        """
        m = len(data)
        if self.n < m:
            return float('-inf')

        x = self._random()
        y = np.log(x[:m]) * data
        return y.sum()

    def _marginal_beta(self, i):
        """
        Computes the marginal distribution of the ith element.

        See http://en.wikipedia.org/wiki/Dirichlet_distribution#Marginal_distributions

        Args:
            i: int

        Returns:
            Beta object
        """
        alpha0 = self.params.sum()
        alpha = self.params[i]
        return Beta(alpha, alpha0 - alpha)

    def _predictive_pmf(self, xs, name=''):
        """
        Makes a predictive distribution.

        Args:
            xs: values to go into the Pmf

        Returns:
            Pmf that maps from x to the mean prevalence of x
        """
        alpha0 = self.params.sum()
        ps = self.params / alpha0
        return _make_pmf_from_items(zip(xs, ps), name=name)


def _binomial_coef(n, k):
    """
    Compute the binomial coefficient "n choose k".

    Args:
        n: number of trials
        k: number of successes

    Returns:
        float
    """
    return scipy.misc.comb(n, k)


def _log_binomial_coef(n, k):
    """
    Computes the log of the binomial coefficient.

    http://math.stackexchange.com/questions/64716/approximating-the-logarithm-of-the-binomial-coefficient

    Args:
        n: number of trials
        k: number of successes

    Returns:
        float
    """
    return n * math.log(n) - k * math.log(k) - (n - k) * math.log(n - k)


def _normal_probability(ys, jitter=0.0):
    """
    Generates data for a normal probability plot.

    Args:
        ys:     sequence of values
        jitter: float magnitude of jitter added to the ys

    returns:
        xs, ys
    """
    n = len(ys)
    xs = np.random.normal(0, 1, n)
    xs.sort()

    if jitter:
        ys = np.random.uniform(-jitter, +jitter, n) + ys
    ys.sort()

    return xs, ys


def _jitter(values, jitter=0.5):
    """Jitters the values by adding a uniform variate in (-jitter, jitter)."""
    return np.random.uniform(-jitter, +jitter, n) + values


def _fit_line(xs, inter, slope):
    """
    Fits a line to the given data.

    Args:
        xs: sequence of x

    returns:
        tuple of numpy arrays (sorted xs, fit ys)
    """
    fit_xs = np.sort(xs)
    fit_ys = inter + slope * fit_xs
    return fit_xs, fit_ys


def _normal_probability_plot(sample, label, data_color='blue', fit_color='gray'):
    """
    Makes a normal probability plot with a fitted line.

    Args:
        sample:     sequence of numbers
        label:      string
        data_color: color string for the data
        fit_color:  color string for the fitted line
    """
    data = _normal_probability(sample)
    fit = _fit_line(*data)

    thinkplot.plot(*fit, color=fit_color, alpha=0.5)
    thinkplot.plot(*data,
                   label=label,
                   color=data_color,
                   marker='.',
                   markersize=5,
                   alpha=0.5)


def _cov(xs, ys, mux=None, muy=None):
    """
    Computes Cov(X, Y).

    Args:
        xs:  sequence of values
        ys:  sequence of values
        mux: optional float mean of xs
        muy: optional float mean of ys

    Returns:
        Cov(X, Y)
    """
    if mux is None:
        mux = np.mean(xs)
    if muy is None:
        muy = np.mean(ys)

    total = 0.0
    for x, y in zip(xs, ys):
        total += (x - mux) * (y - muy)

    return total / len(xs)


def _mean(xs):
    """
    Computes mean.

    Args:
        xs: sequence of values

    returns:
        float mean
    """
    return np.mean(xs)


def _var(xs, ddof=None):
    """
    Computes variance.

    Args:
        xs: sequence of values

    returns:
        float
    """
    return np.var(xs, ddof=ddof)


def _mean_var(xs):
    """
    Computes mean and variance.

    Args:
        xs: sequence of values

    returns:
        pair of float, mean and var
    """
    return np.mean(xs), np.var(xs)


def _trim(t, p=0.01):
    """
    Trims the largest and smallest elements of t.

    Args:
        t: sequence of numbers
        p: fraction of values to trim off each end

    Returns:
        sequence of values
    """
    n = int(p * len(t))
    t = sorted(t)[n:-n]
    return t


def _trimmed_mean(t, p=0.01):
    """
    Computes the trimmed mean of a sequence of numbers.

    Side effect: sorts the list.

    Args:
        t: sequence of numbers
        p: fraction of values to trim off each end

    Returns:
        float
    """
    t = _trim(t, p)
    return _mean(t)


def _trimmed_mean_var(t, p=0.01):
    """
    Computes the trimmed mean and variance of a sequence of numbers.

    Side effect: sorts the list.

    Args:
        t: sequence of numbers
        p: fraction of values to trim off each end

    Returns:
        float
    """
    t = _trim(t, p)
    mu, var = _mean_var(t)
    return mu, var


def _corr(xs, ys):
    """
    Computes Corr(X, Y).

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        Corr(X, Y)
    """
    xbar, varx = _mean_var(xs)
    ybar, vary = _mean_var(ys)

    corr = _cov(xs, ys, xbar, ybar) / math.sqrt(varx * vary)

    return corr


def _serial_corr(xs):
    """
    Computes the serial correlation of a sequence.

    Args:
        xs: sequence of numbers

    returns:
        float correlation coefficient
    """
    return _corr(xs[:-1], xs[1:])


def _spearman_corr(xs, ys):
    """
    Computes Spearman's rank correlation.

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        float Spearman's correlation
    """
    xranks = _map_to_ranks(xs)
    yranks = _map_to_ranks(ys)
    return _corr(xranks, yranks)


def _least_squares(xs, ys):
    """
    Computes a linear least squares fit for ys as a function of xs.

    Args:
        xs: sequence of values
        ys: sequence of values

    Returns:
        tuple of (intercept, slope)
    """
    xbar, varx = _mean_var(xs)
    ybar, vary = _mean_var(ys)

    slope = _cov(xs, ys, xbar, ybar) / varx
    inter = ybar - slope * xbar

    return inter, slope


def _residuals(xs, ys, inter, slope):
    """
    Computes residuals for a linear fit with parameters inter and slope.

    Args:
        xs:    independent variable
        ys:    dependent variable
        inter: float intercept
        slope: float slope

    Returns:
        list of residuals
    """
    res = [y - inter - slope * x for x, y in zip(xs, ys)]
    return res


def _coef_determination(ys, res):
    """
    Computes the coefficient of determination (R^2) for given residuals.

    Args:
        ys:  dependent variable
        res: residuals
        
    Returns:
        float coefficient of determination
    """
    ybar, vary = _mean_var(ys)
    resbar, varres = _mean_var(res)
    return 1 - varres / vary


def _map_to_ranks(t):
    """
    Returns a list of ranks corresponding to the elements in t.

    Args:
        t: sequence of numbers
    
    Returns:
        list of integer ranks, starting at 1
    """
    # pair up each value with its index
    pairs = enumerate(t)

    # sort by value
    sorted_pairs = sorted(pairs, key=lambda pair: pair[1])

    # pair up each pair with its rank
    ranked = enumerate(sorted_pairs)

    # sort by index
    resorted = sorted(ranked, key=lambda trip: trip[1][0])

    # extract the ranks
    ranks = [trip[0] + 1 for trip in resorted]
    return ranks


def _correlated_generator(rho):
    """
    Generates standard normal variates with serial correlation.

    Args:
        rho: target coefficient of correlation

    Returns:
        iterable
    """
    x = random.gauss(0, 1)
    yield x

    sigma = math.sqrt(1 - rho ** 2);
    while True:
        x = random.gauss(x * rho, sigma)
        yield x


def _correlated_gaussian_generator(mu, sigma, rho):
    """
    Generates normal variates with serial correlation.

    Args:
        mu:    mean of variate
        sigma: standard deviation of variate
        rho:   target coefficient of correlation

    Returns:
        iterable
    """
    for x in _correlated_generator(rho):
        yield x * sigma + mu


def _raw_moment(xs, k):
    """Computes the kth raw moment of xs."""
    return sum(x ** k for x in xs) / float(len(xs))


def _central_moment(xs, k):
    """Computes the kth central moment of xs."""
    xbar = _raw_moment(xs, 1)
    return sum((x - xbar) ** k for x in xs) / len(xs)


def _standardized_moment(xs, k):
    """Computes the kth standardized moment of xs."""
    var = _central_moment(xs, 2)
    sigma = math.sqrt(var)
    return _central_moment(xs, k) / sigma ** k


def _skewness(xs):
    """Computes skewness."""
    return _standardized_moment(xs, 3)


def _median(xs):
    """Computes the median (50th percentile) of a sequence."""
    cdf = _make_cdf_from_list(xs)
    return cdf._value(0.5)


def _pearson_median_skewness(xs):
    """Computes the Pearson median skewness."""
    median = _median(xs)
    mean = _raw_moment(xs, 1)
    var = _central_moment(xs, 2)
    std = math.sqrt(var)
    gp = 3 * (mean - median) / std
    return gp


class Dictionary(object):
    """Represents a set of variables in a fixed width file."""

    def __init__(self, variables, colspecs, names):
        """
        Initializes.

        Args:
            variables: list of (start, vtype, name, fstring, long_desc) tuples
            colspecs:  list of (start, end) index tuples
            names:     list of string variable names
        """
        self.variables = variables
        self.colspecs = colspecs
        self.names = names

    def _read_fixed_width(self, dat_file, compression='gzip'):
        """
        Reads a fixed width ASCII file.

        Args:
            dat_file:    string filename
            compression: string

        returns:
            DataFrame
        """
        frame = pd.read_fwf(dat_file,
                            compression=compression,
                            colspecs=self.colspecs,
                            names=self.names,
                            header=None)
        return frame


def _read_stata_dct(dct_file):
    """
    Reads a Stata dictionary file.

    returns:
        Dictionary object
    """
    type_map = dict(byte=int, int=int, float=float, double=float)

    variables = []
    for line in open(dct_file):
        match = re.search(r'_column\(([^)]*)\)', line)
        if match:
            start = int(match.group(1))
            t = line.split()
            vtype, name, fstring = t[1:4]
            if vtype.startswith('str'):
                vtype = str
            else:
                vtype = type_map[vtype]
            long_desc = ' '.join(t[4:]).strip('"')
            variables.append((start, vtype, name, fstring, long_desc))

    colspecs = []
    names = []
    for i in range(len(variables)):
        start, vtype, name, fstring, long_desc = variables[i]
        names.append(name)
        try:
            end = variables[i + 1][0]
            colspecs.append((start - 1, end - 1))
        except IndexError:
            # Note: this won't work properly until Pandas Issue 7079 is
            # resolved so pd.read_fwf accepts None as a colspec

            # In the meantime, it lops one character off the end of the
            # last field.

            # TODO: replace -1 with None (see DocString above)
            colspecs.append((start - 1, -1))

    return Dictionary(variables, colspecs, names)


def main():
    pass


if __name__ == '__main__':
    main()
