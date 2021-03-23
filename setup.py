import os
from distutils.core import setup
from distutils.command.sdist import sdist

import skyfield  # safe, because __init__.py contains no import statements

class my_sdist(sdist):
    def make_distribution(self):
        # See https://github.com/skyfielders/python-skyfield/issues/378
        for path in self.filelist.files:
            os.chmod(path, 0o644)
        sdist.make_distribution(self)

setup(
    cmdclass={'sdist': my_sdist},
    name='skyfield',
    version=skyfield.__version__,
    description=skyfield.__doc__.split('\n', 1)[0],
    long_description=open('README.rst', 'rb').read().decode('utf-8'),
    license='MIT',
    author='Brandon Rhodes',
    author_email='brandon@rhodesmill.org',
    url='http://github.com/brandon-rhodes/python-skyfield/',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
    packages=[
        'skyfield',
        'skyfield.data',
        'skyfield.tests',
        ],
    package_data = {
        'skyfield': ['documentation/*.rst'],
        'skyfield.data': ['*.gz', '*.npy', '*.npz'],
        'skyfield.tests': ['data/*'],
        },
    install_requires=[
        'certifi>=2017.4.17',  # version constraint copied from Requests
        'jplephem>=2.13',
        'numpy',
        'sgp4>=2.2',
        ],
)
