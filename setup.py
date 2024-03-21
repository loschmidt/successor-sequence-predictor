from setuptools import setup, find_packages

setup(
    # Self-descriptive entries which should always be present
    name='loschmidt-ssp',
    author='Pavel Kohout, Loschmidt Laboratories',
    author_email='xkohou15@vutbr.cz',
    version="1.0",
    description="Successor Sequence Predictor",
    url="https://github.com/loschmidt/successor-sequence-predictor",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: LGPL License",
        "Operating System :: OS Independent",
        "Environment :: Console",
    ],
    python_requires=">=3.9",
    install_requires=[
        "click~=8.1.7",
        "pandas~=2.2.0",
        "pyyaml~=6.0",
        "numpy==1.26.4",
        "scipy==1.12.0",
        "scikit-learn~=1.3.0",
        "BioPython",
        "seaborn~=0.13"
    ],
    entry_points={
        'console_scripts': [
            'loschmidt-ssp=loschmidt.ssp.cli.main:cli'
        ],
    }
)
