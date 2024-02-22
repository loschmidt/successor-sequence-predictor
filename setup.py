"""
Successor reconstruction setup script
Robust prediction of successor based on multiple phylogeny trees for inference
"""

from setuptools import setup, find_packages

short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name='SuccessorRobustReconstruction',
    author='Pavel Kohout, Loschmidt Laboratories',
    author_email='xkohou15@vutbr.cz',
    version="1.0",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.loschmidt.cz/llws/llws-runners",
    license='MIT license',

    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Environment :: Console",
    ],
    python_requires=">=3.8",
    install_requires=[
        "click~=8.0.4",
        "pyyaml~=6.0",
        "numpy==1.26.4",
        "scipy==1.10.1",
        "scikit-learn~=1.3.0",
        "BioPython",
        "seaborn~=0.13",
        "setuptools~=67.8.0"
    ],
    entry_points={
        'console_scripts': [
            'successor=successors.cli.main:cli'
        ],
    }
)
