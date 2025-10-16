from distutils.core import setup

import skyfield  # safe, because __init__.py contains no import statements

setup(
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
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    packages=[
        'skyfield',
        'skyfield.data',
        'skyfield.tests',
    ],
    package_data = {
        'skyfield.data': ['*.gz', '*.npy', '*.npz'],
        'skyfield.tests': ['data/*'],
    },
    install_requires=[
        'certifi>=2017.4.17',  # version constraint copied from Requests
        'jplephem>=2.13',
        'numpy',
        'sgp4>=2.13',
    ],
)
