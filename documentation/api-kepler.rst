
===============================
 API Reference — Kepler Orbits
===============================

.. currentmodule:: skyfield.data.mpc

Skyfield includes some contributed code
that can compute predicted positions for minor planets and comets
given simple Kepler orbital elements.
The computations are not yet very fast,
and can’t yet compute multiple body positions efficiently
using NumPy array operations—if you ask for ten comet positions,
then Skyfield has to run an internal loop
computing each of the ten positions separately.

The library is not yet fully documented;
see :doc:`kepler-orbits` for an introduction
to the features it already supports.

Here, we have only yet documented the routines for loading data,
as they are the most stable.

.. autofunction:: load_mpcorb_dataframe
.. autofunction:: load_comets_dataframe
.. autofunction:: load_comets_dataframe_slow
