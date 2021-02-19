"""Module that allows the selection of data points interactively.

Classes:
    Highlighter
"""

from matplotlib.widgets import RectangleSelector
import numpy as np


class Highlighter(object):
    """Class that highlights and returns selected objects in plot.

    Attributes:
        ax: (matplotlib.axes._subplots.AxesSubplot).
        x: (numpy.ndarray) X values.
        y: (numpy.ndarray) Y values.
        canvas: (backend_interagg.FigureCanvasInterAgg)
        mask: (numpy.ndarray of booleans) Array that indicates which coordinates are highlighted.
        selector: (matplotlib.widgets.RectangleSelector)
    """

    def __init__(self, ax, x, y):
        """

        :param ax: (matplotlib.axes._subplots.AxesSubplot).
        :param x: (numpy.ndarray) X values.
        :param y: (numpy.ndarray) X values.
        """

        self.ax = ax
        self.canvas = ax.figure.canvas
        self.x, self.y = x, y
        self.mask = np.zeros(x.shape, dtype=bool)

        self._highlight = ax.scatter([], [], s=50, color='yellow', zorder=10)
        self.selector = RectangleSelector(ax, self, useblit=True)

    def __call__(self, event1, event2):
        """The __call__ method enables Python programmers to write classes where the instances behave like functions
        and can be called like a function.
        Draws selected points captures by mouse event to plot.


        :param event1:
        :param event2:
        :return:
        """

        # |= allow points to be selected in turn
        # ^ -> XOR operator
        self.mask ^= self.inside(event1, event2)

        xy = np.column_stack([self.x[self.mask], self.y[self.mask]])

        self._highlight.set_offsets(xy)
        self.canvas.draw()

    def inside(self, event1, event2):
        """Returns a boolean mask of the points inside the rectangle defined by
        event1 and event2.

        :param event1: (matplotlib.backend_bases.MouseEvent)
        :param event2: (matplotlib.backend_bases.MouseEvent)
        :return: (array) Boolean mask of the points inside rectangle.
        """

        # Note: Could use points_inside_poly, as well
        x0, x1 = sorted([event1.xdata, event2.xdata])
        y0, y1 = sorted([event1.ydata, event2.ydata])
        mask = ((self.x > x0) & (self.x < x1) & (self.y > y0) & (self.y < y1))

        # return not selected points
        # return ~mask

        # return selected points
        return mask
