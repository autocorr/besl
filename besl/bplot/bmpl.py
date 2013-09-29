"""
===============
BFigure & BAxes
===============

Baseclass for extensions of the Matplotlib Figure and Axes classes.

"""

import numpy as np
from matplotlib.figure import (Figure, Bbox)
from matplotlib.axes import Axes


class BFigure(Figure):
    def __init__(self, **kwargs):
        super(BFigure, self).__init__(**kwargs)


class BAxes(Axes):
    def __init__(self, fig, rect, **kwargs):
        bbox = Bbox(rect)
        super(BAxes, self).__init__(fig, bbox, **kwargs)


def one_plot(fig=None, rect=None, **kwargs):
    """
    Construct a BAxes instance

    Parameters
    ----------
    fig : matplotlib.figure.Figure, besl.bplot.baxes.BFigure
    rect : numpy.array

    Returns : besl.bplot.baxes.BAxes
    """
    if fig is None:
        fig = BFigure()
    if rect is None:
        rect = np.array([[0.125, 0.1], [0.9, 0.9]])
    return BAxes(fig, rect, **kwargs)
