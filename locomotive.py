"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

This file contains a solution to the locomotive problem adapted from 
Mosteller, Fifty Challenging Problems in Probability:

"A railroad numbers its locomotives in order 1..N.  One day you see a 
locomotive with the number 60.  Estimate how many locomotives the 
railroad has."

NAME: locomotive.py
"""

import matplotlib.pyplot as pyplot
from math import pow
import _04_Pmf
import _05_myplot
import _13_Cdf


def _make_uniform_suite(low, high, steps):
    """
    Makes a PMF that represents a suite of hypotheses with equal p.
    
    Args:
        low:   low end of range
        high:  high end of range
        steps: number of values

    Returns:
        Pmf object
    """
    hypos = [low + (high - low) * i / (steps - 1.0) for i in range(steps)]
    pmf = _04_Pmf._make_pmf_from_list(hypos)
    return pmf


def _update(suite, evidence):
    """
    Updates a suite of hypotheses based on new evidence.

    Modifies the suite directly; if you want to keep the original, make a copy.

    Args:
        suite:    Pmf object
        evidence: whatever kind of object Likelihood expects
    """
    for hypo in suite._values():
        likelihood = _likelihood(evidence, hypo)
        suite._mult(hypo, likelihood)
    suite._normalize()


def _likelihood(evidence, hypo):
    """
    Computes the likelihood of the evidence assuming the hypothesis is true.

    Args:
        evidence: the serial number of the observed train
        hypo:     int hypothetical number of trains

    Returns:
        probability of seeing a given train assuming that there are N trains
    """
    train_seen = evidence
    num_trains = hypo
    if train_seen > num_trains:
        return 0.0
    else:
        return 1.0 / num_trains


def _credible_interval(pmf, percentage):
    """
    Computes a credible interval for a given distribution.

    If percentage=90, computes the 90% CI.

    Args:
        pmf:        Pmf object representing a posterior distribution
        percentage: float between 0 and 100

    Returns:
        sequence of two floats, low and high
    """
    cdf = _13_Cdf._make_cdf_from_dict(pmf._get_dict())
    prob = (1 - percentage / 100.0) / 2
    interval = [cdf._value(p) for p in [prob, 1 - prob]]
    return interval


def main():
    upper_bound = 200
    prior = _make_uniform_suite(1, upper_bound, upper_bound)
    prior.name = 'prior'

    evidence = 60
    posterior = prior._copy()
    _update(posterior, evidence)
    posterior.name = 'posterior'

    print(_credible_interval(posterior, 90))

    # plot the posterior distribution
    pyplot.subplots_adjust(wspace=0.4, left=0.15)
    plot_options = dict(linewidth=2)

    _05_myplot._pmf(posterior, **plot_options)
    _05_myplot._save(root='locomotive',
                     title='Locomotive problem',
                     xlabel='Number of trains',
                     ylabel='Posterior probability')


if __name__ == '__main__':
    main()
