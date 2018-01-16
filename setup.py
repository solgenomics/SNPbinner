#!/usr/bin/env python

from setuptools import setup

setup(name='snpbinner',
    version='0.2.3',
    description='',
    author='David Lyon',
    author_email='dlyon@fandm.edu',
    url='https://github.com/solgenomics/snpbinner',
    install_requires=['pillow'],
    packages=['snpbinner'],
    include_package_data=True,
    package_data={'snpbinner':['package_data/Tuffy.ttf']},
    entry_points = {
        'console_scripts': ['snpbinner = snpbinner.__main__:main'],
    }
)
