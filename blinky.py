"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

This file contains a solution to the Blink Monty Problem, by
Allen Downey:

Suppose you are on Let's Make a Deal and you are playing the Monty
Hall Game, with one difference: before you went on the show you
analyzed tapes of previous shows and discovered that Monty has a tell:
when the contestant picks the correct door, Monty is more likely to
blink.

Specifically, of the 15 shows you watched, the contestant chose the
correct door 5 times, and Monty blinked three of those times.  Of the
other 10 times, Monty blinked three times.

Assume that you choose Door A.  Monty opens door B and blinks.  What
should you do, and what is your chance of winning?


You can read a discussion of this problem at XXX

NAME: blinky.py
"""

import _04_Pmf
import _05_myplot


def _make_uniform_suite(low, high, steps, name=''):
    """
    Makes a PMF that represents a suite of hypotheses with equal p.
    
    Args:
        low:   low end of range
        high:  high end of range
        steps: number of values
        name:  string name for the Pmf

    Returns:
        Pmf object
    """
    hypos = [low + (high - low) * i / (steps - 1.0) for i in range(steps)]
    pmf = _04_Pmf._make_pmf_from_list(hypos, name=name)
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
        evidence: a tuple of (number of successes, number of failures)
        hypo:     float probability of success

    Returns:
        unnormalized likelihood of getting the given number of successes
        and failures if the probability of success is p
    """
    heads, tails = evidence
    p = hypo
    return pow(p, heads) * pow(1 - p, tails)


def _total_probability(pmf1, pmf2, func):
    """
    Enumerates pairs from the Pmfs, calls the func, and returns
    the total probability.

    Args:
        pmf1: Pmf object
        pmf2: Pmf object
        func: a callable that takes a value from each Pmf and returns probability.
    """
    total = 0.0
    for x, px in pmf1._items():
        for y, py in pmf2._items():
            if px and py:
                total += px * py * func(x, y)
    return total


def _prob_winning(pbA, pbC):
    """
    Computes the probability that the car is behind door A:

    Args:
        pbA: probability that Monty blinks if the car is behind A
        pbC: probability that Monty blinks if the car is behind C
    """
    pea = 0.5 * pbA
    pec = pbC

    pae = pea / (pea + pec)
    return pae


def main():
    print('pae', 0.3 / (0.3 + 3.0 / 13))

    doorA = _make_uniform_suite(0.0, 1.0, 101, name='Door A')
    evidence = 3, 2
    _update(doorA, evidence)

    doorC = _make_uniform_suite(0.0, 1.0, 101, name='Door C')
    evidence = 3, 10
    _update(doorC, evidence)

    print(_total_probability(doorA, doorC, _prob_winning))

    # plot the posterior distributions
    _05_myplot._pmfs([doorA, doorC])
    _05_myplot._save(root='blinky',
                     formats=['pdf', 'png'],
                     title='Probability of blinking',
                     xlabel='P(blink)',
                     ylabel='Posterior probability')


if __name__ == '__main__':
    main()
