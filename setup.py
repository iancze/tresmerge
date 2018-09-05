import setuptools
from codecs import open
from os import path

import numpy as np

# set the command line scripts
entry_points = {'console_scripts': [
    'tresmerge-process = tresmerge.process:main'
    ]}


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tresmerge",
    version="0.0.1",
    author="Ian Czekala",
    author_email="iancze@gmail.com",
    description="Flatten and merge echelle orders of the TRES echelle spectrograph.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/iancze/TRES-merge",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy','scipy','astropy','specutils','matplotlib'],
    entry_points=entry_points,
    package_data={'tresmerge':['data/sensfuncs.npy']} # the sensitivity curves for fluxcal
)
