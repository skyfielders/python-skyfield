===================
Hacking on Skyfield
===================

If you would like to get started hacking on the codebase, we highly suggest that you use Virtualenv_. We also suggest checking out VirtualenvWrapper_.

=================
Running The Tests
=================

First install the developer dependencies into your current environment with::
    pip install -e .
    pip install -r requirements.txt

Now run the tests against your current environment with::
    py.test

To run the tests against all of the platforms that Skyfield supports, use tox_::
    pip install tox
    tox

.. _Virtualenv: http://www.virtualenv.org/en/latest/
.. _VirtualenvWrapper: http://doughellmann.com/2008/05/virtualenvwrapper.html
.. _tox: http://tox.readthedocs.org/en/latest/
