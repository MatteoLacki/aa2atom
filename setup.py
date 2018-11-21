# This Python file uses the following encoding: utf-8
from setuptools import setup, find_packages

setup(
    name='aa2atom',
    packages=find_packages(),
    version='0.1',
    description='Change amino-acidic sequence into counts of atoms.',
    long_description='Change amino-acidic sequence into counts of atoms.',
    author='Mateusz Krzysztof Łącki',
    author_email='matteo.lacki@gmail.com',
    url='',
    # download_url='https://github.com/MatteoLacki/MassTodonPy/tree/GutenTag',
    keywords=[
        'Mass Spectrometry',
        'fasta'
        'reverse fastas'],
    classifiers=[
        'Development Status :: 1 - Planning',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3.6'],
    install_requires=[
        'linearCounter',
    ],
    # scripts = [
    #     "bin/update_tenzer"
    # ]
)
