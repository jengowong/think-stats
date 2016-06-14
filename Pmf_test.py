"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: Pmf_test.py
"""

import unittest
import _03_thinkstats
import _04_Pmf
import _13_Cdf


class Test(unittest.TestCase):
    def testHist(self):
        t = [1, 2, 2, 3, 5]
        hist = _04_Pmf._make_hist_from_list(t)

        self.assertEquals(hist._freq(1), 1)
        self.assertEquals(hist._freq(2), 2)
        self.assertEquals(hist._freq(3), 1)
        self.assertEquals(hist._freq(4), 0)
        self.assertEquals(hist._freq(5), 1)

        pmf = _04_Pmf._make_pmf_from_hist(hist)

        pmf._print()
        self.checkPmf(pmf)

    def checkPmf(self, pmf):
        self.assertAlmostEquals(pmf._prob(1), 0.2)
        self.assertAlmostEquals(pmf._prob(2), 0.4)
        self.assertAlmostEquals(pmf._prob(3), 0.2)
        self.assertAlmostEquals(pmf._prob(4), 0.0)
        self.assertAlmostEquals(pmf._prob(5), 0.2)

    def testMakePmf(self):
        t = [1, 2, 2, 3, 5]
        pmf = _04_Pmf._make_pmf_from_list(t)
        self.checkPmf(pmf)

        d = pmf._get_dict()
        self.assertAlmostEquals(d[2], 0.4)

        vals = pmf._values()
        self.assertEquals(sorted(vals), [1, 2, 3, 5])

        items = pmf._items()
        d = dict(items)
        new_pmf = _04_Pmf._make_pmf_from_dict(d)
        self.checkPmf(new_pmf)

    def testSetAndNormalize(self):
        pmf = _04_Pmf.Pmf()
        t = [1, 2, 2, 3, 5]
        for x in t:
            pmf._set(x, 1)
        pmf._incr(2)
        pmf._normalize()
        self.checkPmf(pmf)

    def testIncrAndNormalize(self):
        pmf = _04_Pmf.Pmf()
        t = [1, 2, 2, 3, 5]
        for x in t:
            pmf._incr(x)
        pmf._normalize()
        self.checkPmf(pmf)

    def testMultAndNormalize(self):
        t = [1, 2, 3, 5]
        pmf = _04_Pmf._make_pmf_from_list(t)
        pmf._mult(2, 2)
        pmf._normalize()
        self.checkPmf(pmf)

    def testRender(self):
        t = [1, 2, 2, 3, 5]
        pmf = _04_Pmf._make_pmf_from_list(t)
        xs, ps = pmf._render()

        d = dict(zip(xs, ps))
        new_pmf = _04_Pmf._make_pmf_from_dict(d)
        self.checkPmf(new_pmf)

    def testMeanAndVar(self):
        t = [1, 2, 2, 3, 5]
        mu = _03_thinkstats._mean(t)
        var = _03_thinkstats._var(t, mu)

        pmf = _04_Pmf._make_pmf_from_list(t)
        mu2 = pmf._mean()
        var2 = pmf._var()
        var3 = pmf._var(mu2)

        self.assertAlmostEquals(mu, mu2)
        self.assertAlmostEquals(var, var2)
        self.assertAlmostEquals(var, var3)

    def testMakePmfFromCdf(self):
        t = [1, 2, 2, 3, 5]
        pmf = _04_Pmf._make_pmf_from_list(t)
        self.checkPmf(pmf)

        cdf = _13_Cdf._make_cdf_from_pmf(pmf)
        pmf2 = _04_Pmf._make_pmf_from_cdf(cdf)
        self.checkPmf(pmf2)


if __name__ == "__main__":
    unittest.main()
