"""Various curves."""

from numpy import arange
from ._compatibility import interp
from .descriptorlib import reify
from .functions import _to_array

class Splines(object):
    def __init__(self, table):
        table = _to_array(table)
        if len(table.shape) < 2:  # Let caller provide a single row.
            table = table.reshape(table.shape + (1,))
        self.table = table
        self.lower = lower = table[0]
        self.upper = upper = table[1]
        self._width = upper - lower
        self._n = arange(len(lower))
        self.coefficients = table[2:]

    def __call__(self, x):
        i = interp(x, self.lower, self._n)
        i = i.astype(int)
        t = (x - self.lower[i]) / self._width[i]
        coefficients = iter(self.coefficients)
        value = next(coefficients)[i]
        for c in coefficients:
            value *= t
            value += c[i]
        return value

    @reify
    def derivative(self):
        columns = [self.table[0], self.table[1]]
        coefficients = self.table[2:-1]
        for i, c in enumerate(coefficients):
            n = len(coefficients) - i
            columns.append(n * c / self._width)
        return Splines(columns)

def build_spline_given_ends(x0, y0, slope0, x1, y1, slope1):
    width = x1 - x0
    slope0 = slope0 * width
    slope1 = slope1 * width
    a0 = y0
    a1 = slope0
    a2 = -2*slope0 - slope1 - 3*y0 + 3*y1
    a3 = slope0 + slope1 + 2*y0 - 2*y1
    return x0, x1, a3, a2, a1, a0
