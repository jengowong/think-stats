"""
This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

NAME: brfss_scatter.py
"""

import math
import matplotlib
import matplotlib.pyplot as pyplot
import random
import sys
import _05_myplot
import _19_brfss


class Respondents(_19_brfss.Respondents):
    """Represents the respondent table."""

    def _get_height_weight(self, jitter=0.0):
        """
        Get sequences of height and weight.

        Args:
            jitter: float magnitude of random noise added to heights

        Returns:
            tuple of sequences (heights, weights)
        """
        heights = []
        weights = []
        for r in self.records:
            if r.wtkg2 == 'NA' or r.htm3 == 'NA':
                continue

            height = r.htm3 + random.uniform(-jitter, jitter)

            heights.append(height)
            weights.append(r.wtkg2)

        return heights, weights

    @staticmethod
    def _scatter_plot(root, heights, weights, alpha=1.0):
        pyplot.scatter(heights, weights, alpha=alpha, edgecolors='none')
        _05_myplot._save(root=root,
                         xlabel='Height (cm)',
                         ylabel='Weight (kg)',
                         axis=[140, 210, 20, 200],
                         legend=False)

    @staticmethod
    def _hex_bin(root, heights, weights, cmap=matplotlib.cm.Blues):
        pyplot.hexbin(heights, weights, cmap=cmap)
        _05_myplot._save(root=root,
                         xlabel='Height (cm)',
                         ylabel='Weight (kg)',
                         axis=[140, 210, 20, 200],
                         legend=False)


def _make_figures():
    resp = Respondents()
    resp._read_records(n=1000)

    heights, weights = resp._get_height_weight(jitter=0.0)
    pyplot.clf()
    resp._scatter_plot('scatter1', heights, weights)

    heights, weights = resp._get_height_weight(jitter=1.3)
    pyplot.clf()
    resp._scatter_plot('scatter2', heights, weights)

    pyplot.clf()
    resp._scatter_plot('scatter3', heights, weights, alpha=0.2)

    # read more respondents for the hexplot
    resp = Respondents()
    resp._read_records(n=10000)
    heights, weights = resp._get_height_weight(jitter=1.3)

    pyplot.clf()
    resp._hex_bin('scatter4', heights, weights)


def main(name):
    _make_figures()


if __name__ == '__main__':
    main(*sys.argv)
