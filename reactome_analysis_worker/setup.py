#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

requirements = ["reactome-analysis-utils", "rpy2 > 2.9.0, < 3.0.0", "statsmodels", "scipy", "prometheus_client"]

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
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="The worker process of the REACTOME analysis system performing the actual GSA.",
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description="The worker process of the REACTOME analysis system performing the actual GSA.",
    include_package_data=True,
    keywords='reactome_analysis_worker',
    name='reactome_analysis_worker',
    packages=find_packages(),
    package_data={'': ['*.R', 'IntAct_Static.txt']},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/jgriss/reactome_analysis_worker',
    version='0.1.0',
    zip_safe=False,
)
