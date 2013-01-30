from distutils.core import setup
import skyfield  # to learn the version

setup(
    name='skyfield',
    version=skyfield.__version__,
    description=skyfield.__doc__,
    long_description=open('README.rst').read(),
    license='MIT',
    author='Brandon Rhodes',
    author_email='brandon@rhodesmill.org',
    url='http://github.com/brandon-rhodes/skyfield/',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
    packages=[ 'skyfield', 'skyfield.tests' ],
    install_requires=['jplephem', 'numpy', 'sgp4'],
    )
