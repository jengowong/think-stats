"""
This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import _04_Pmf
import _05_myplot


def BiasPmf(pmf, name, invert=False):
    """
    Returns the Pmf with oversampling proportional to value.

    If pmf is the distribution of true values, the result is the
    distribution that would be seen if values are oversampled in
    proportion to their values; for example, if you ask students
    how big their classes are, large classes are oversampled in
    proportion to their size.

    If invert=True, computes in inverse operation;
    for example, unbiasing a sample collected from students.

    Args:
      pmf:    Pmf object.
      name:   string name for the new Pmf.
      invert: boolean

     Returns:
       Pmf object
    """
    new_pmf = pmf._copy()
    new_pmf.name = name

    for x, p in pmf._items():
        if invert:
            new_pmf._mult(x, 1.0 / x)
        else:
            new_pmf._mult(x, x)

    new_pmf._normalize()
    return new_pmf


def UnbiasPmf(pmf, name):
    """
    Returns the Pmf with oversampling proportional to 1/value.

    Args:
      pmf:  Pmf object.
      name: string name for the new Pmf.

     Returns:
       Pmf object
    """
    return BiasPmf(pmf, name, invert=True)


def ClassSizes():
    # start with the actual distribution of class sizes from the book
    d = {
        7: 8,
        12: 8,
        17: 14,
        22: 4,
        27: 6,
        32: 12,
        37: 8,
        42: 3,
        47: 2,
    }

    # form the pmf
    pmf = _04_Pmf._make_pmf_from_dict(d, 'actual')
    print('mean', pmf._mean())
    print('var', pmf._var())

    # compute the biased pmf
    biased_pmf = BiasPmf(pmf, 'observed')
    print('mean', biased_pmf._mean())
    print('var', biased_pmf._var())

    # unbias the biased pmf
    unbiased_pmf = UnbiasPmf(biased_pmf, 'unbiased')
    print('mean', unbiased_pmf._mean())
    print('var', unbiased_pmf._var())

    # plot the Pmfs
    _05_myplot.Pmfs([pmf, biased_pmf])
    _05_myplot.Show(xlabel='Class size', ylabel='PMF')


def main():
    ClassSizes()


if __name__ == '__main__':
    main()
