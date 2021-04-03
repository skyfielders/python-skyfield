import numpy as np

# Old versions of interp() sometimes return a plain Python float,
# instead of always returning at least a NumPy float64 object.

interp = np.interp

if isinstance(interp(1.0, [0.0, 2.0], [10.0, 20.0]), float):
    def interp(*args):
        return np.array(np.interp(*args))
