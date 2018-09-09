from distutils.core import setup
import skyfield  # safe, because __init__.py contains no import statements

extras = {
    'tests': [
        # TODO: add assay
        'pytz',
    ],
}

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
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
    packages=[
        'skyfield',
        'skyfield.data',
        'skyfield.tests',
        ],
    package_data = {
        'skyfield': ['documentation/*.rst'],
        'skyfield.data': ['*.npy', '*.npz'],
        'skyfield.tests': ['data/*'],
        },
    install_requires=[
        'jplephem>=2.3',
        'numpy',
        'sgp4>=1.4',
        ],
    extras_require=extras,          # support "pip install skyfield[tests]"
)
