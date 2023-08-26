#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

requirements = ['redis', 'pika']

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Johannes Griss",
    author_email='jgriss@ebi.ac.uk',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Basic utilities used by the REACTOME Analysis system",
    entry_points={
        'console_scripts': [
            'reactome_analysis_utils=reactome_analysis_utils.cli:main',
        ],
    },
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description="Shared util package of the Reactome Analysis System",
    include_package_data=True,
    keywords='reactome_analysis_utils',
    name='reactome_analysis_utils',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/jgriss/reactome_analysis_utils',
    version='0.1.0',
    zip_safe=False,
)
