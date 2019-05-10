# This Python file uses the following encoding: utf-8
from setuptools import setup, find_packages

setup(
    name='aa2atom',
    packages=find_packages(),
    version='1.2',
    description='Change amino-acidic sequence into counts of atoms.',
    long_description='Change amino-acidic sequence into counts of atoms.',
    author='Mateusz Krzysztof Łącki',
    author_email='matteo.lacki@gmail.com',
    keywords=['Mass Spectrometry', 'fasta', 'Analitical Chemistry', 'Amino acids'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.5'],
    install_requires=[],
    scripts = [
        "bin/aa2atom"
    ]
)
