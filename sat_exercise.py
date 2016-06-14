"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: sat_exercise.py
"""

import csv
import sys
import _03_thinkstats
import _04_Pmf
import _05_myplot


def _read_ranks(filename='sat_ranks.csv'):
    """
    Reads a CSV file of SAT scores.

    Args:
        filename: string filename

    Returns:
        list of (score, number) pairs
    """
    fp = open(filename)
    reader = csv.reader(fp)
    res = []

    for t in reader:
        try:
            score = int(t[0])
            number = int(t[1])
            res.append((score, number))
        except ValueError:
            pass

    return res


def _read_scale(filename='sat_scale.csv', col=2):
    """
    Reads a CSV file of SAT scales (maps from raw score to standard score).

    Args:
        filename: string filename
        col:      which column to start with (0=Reading, 2=Math, 4=Writing)

    Returns:
        list of (raw score, standardize score) pairs
    """

    def ParseRange(s):
        t = [int(x) for x in s.split('-')]
        return 1.0 * sum(t) / len(t)

    fp = open(filename)
    reader = csv.reader(fp)
    raws = []
    scores = []

    for t in reader:
        try:
            raw = int(t[col])
            raws.append(raw)
            score = ParseRange(t[col + 1])
            scores.append(score)
        except:
            pass

    raws.sort()
    scores.sort()
    return _03_thinkstats.Interpolator(raws, scores)


def _reverse_scale(pmf, scale):
    """
    Applies the reverse scale to the values of a PMF.

    Args:
        pmf:   Pmf of scaled scores
        scale: Interpolator object

    Returns:
        Pmf of raw scores
    """
    new = _04_Pmf.Pmf()
    for val, prob in pmf._items():
        raw = scale._reverse(val)
        new._incr(raw, prob)
    return new


def _divide_values(pmf, denom):
    """Divides the values in a PMF by denom.  Returns a new PMF."""
    new = _04_Pmf.Pmf()
    for val, prob in pmf._items():
        if val >= 0:
            x = 1.0 * val / denom
            new._incr(x, prob)
    return new


class Exam:
    """
    Encapsulates information about an exam.

    Contains the distribution of scaled scores and an
    Interpolator that maps between scaled and raw scores.
    """

    def __init__(self):
        # scores is a list of (scaled score, number pairs)
        scores = _read_ranks()

        # hist is the histogram of scaled scores
        hist = _04_Pmf._make_hist_from_dict(dict(scores))

        # scaled is the PMF of scaled scores
        self.scaled = _04_Pmf._make_pmf_from_hist(hist)

        # scale is an Interpolator from raw scores to scaled scores
        self.scale = _read_scale()

        # raw is the PMF of raw scores
        self.raw = _reverse_scale(self.scaled, self.scale)

        # max_score is the highest raw score
        self.max_score = max(self.raw._values())

    def _get_raw_score(self, scaled_score):
        """Looks up a scaled score and returns a raw score."""
        return self.scale._reverse(scaled_score)

    def _get_prior(self):
        """Returns a new PMF of p, which is (raw_score / max_score)."""
        prior = _divide_values(self.raw, denom=self.max_score)
        return prior

    def _get_max_score(self):
        """Returns the highest raw score, presumed to the the max possible."""
        return self.max_score


def main(script):
    # make an exam object with data from the 2010 SAT
    exam = Exam()

    # look up Alice's raw score
    alice = 780
    alice_correct = exam._get_raw_score(alice)
    print('Alice raw score', alice_correct)

    # display the distribution of raw scores for the population
    prior = exam._get_prior()
    _05_myplot._pmf(prior, show=True)


if __name__ == '__main__':
    main(*sys.argv)
