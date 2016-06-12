"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import cPickle
import continuous
import math
import matplotlib
import matplotlib.pyplot as pyplot
import random
import sys
import _01_survey
import _03_thinkstats
import _05_myplot
import _13_Cdf
import _17_rankit
import _19_brfss


class Respondents(_19_brfss.Respondents):
    """Represents the respondent table."""

    def MakeNormalModel(self,
                        weights,
                        root,
                        xmax=175,
                        xlabel='adult weight (kg)',
                        axis=None):
        cdf = _13_Cdf.MakeCdfFromList(weights)

        pyplot.clf()

        t = weights[:]
        t.sort()
        mu, var = _03_thinkstats.TrimmedMeanVar(t)
        print('n, Mean, Var', len(weights), mu, var)

        sigma = math.sqrt(var)
        print('Sigma', sigma)

        xs, ps = continuous.RenderNormalCdf(mu, sigma, xmax)
        pyplot.plot(xs, ps, label='model', linewidth=4, color='0.7')

        xs, ps = cdf.Render()
        pyplot.plot(xs, ps, label='data', linewidth=2, color='blue')

        _05_myplot.Save(root,
                        title='Adult weight',
                        xlabel=xlabel,
                        ylabel='CDF',
                        axis=axis or [0, xmax, 0, 1])

    def MakeFigures(self):
        """Generates CDFs and normal prob plots for weights and log weights."""
        weights = [record.wtkg2 for record in self.records if record.wtkg2 != 'NA']
        self.MakeNormalModel(weights, root='brfss_weight_model')
        _17_rankit.MakeNormalPlot(weights,
                                  root='brfss_weight_normal',
                                  title='Adult weight',
                                  ylabel='Weight (kg)')

        log_weights = [math.log(weight) for weight in weights]
        xmax = math.log(175.0)
        axis = [3.5, 5.2, 0, 1]
        self.MakeNormalModel(log_weights,
                             root='brfss_weight_log',
                             xmax=xmax,
                             xlabel='adult weight (log kg)',
                             axis=axis)
        _17_rankit.MakeNormalPlot(log_weights,
                                  root='brfss_weight_lognormal',
                                  title='Adult weight',
                                  ylabel='Weight (log kg)')


def main(name):
    resp = Respondents()
    resp.ReadRecords()
    resp.MakeFigures()


if __name__ == '__main__':
    main(*sys.argv)
