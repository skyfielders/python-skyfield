
==============
 Installation
==============

Skyfield has only a single binary dependency,
the `NumPy <http://www.numpy.org/>`_ vector library,
and is designed to install cleanly with a single invocation
of the standard Python package installation tool::

    pip install skyfield

Note that Skyfield has not yet reached 1.0
so tweaks to the API are still possible.
While everything you read in the documentation should keep working,
you can protect your project from any abrupt changes
by specifying its Skyfield dependency with an explicit version
in your ``requirements.txt`` or ``setup.py`` or install instructions::

    skyfield==0.1

This lets you choose when you opt-in to future versions
and not be blindsided by improvements that take place
between now and the eventual 1.0 version.
If you find any problems or would like to suggest an improvement,
simply create an issue on the project’s GitHub page:

    https://github.com/brandon-rhodes/python-skyfield

If you try to install Skyfield but get errors during NumPy installation,
remember that there are several ways to get NumPy for your system
before you then install Skyfield as a second step.
Here are your options, in order from easiest to hardest:

* It is best to simply use
  `a science-ready Python distribution
  <http://www.scipy.org/install.html#scientific-python-distributions>`_
  like `Anaconda <http://docs.continuum.io/anaconda/install.html>`_
  or `Python(x,y) <https://code.google.com/p/pythonxy/>`_
  because they come with NumPy built in.

* Or you can follow the standard install instructions for the
  `full SciPy stack <http://www.scipy.org/install.html>`_
  or `directly download NumPy
  <https://sourceforge.net/projects/numpy/files/NumPy/>`_
  for your normal install of Python.

* If you are going to let ``pip`` try to install NumPy,
  then make sure that your system has a functioning C compiler
  and that it matches the compiler used on your version of Python.
  Windows users, remember that the Python you download
  from the official web site is compiled with Visual Studio 2008
  and that you can download and install
  `Visual C++ Express 2008 <http://go.microsoft.com/?linkid=7729279>`_
  from Microsoft for free.
  Mac users should install the “Xcode Command Line Tools”
  to give ``pip`` the superpower of being able to build and install
  binary Python dependencies like NumPy.

Good luck!
