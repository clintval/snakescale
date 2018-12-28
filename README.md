# snakescale

[![Testing Status](https://travis-ci.org/clintval/snakescale.svg?branch=master)](https://travis-ci.org/clintval/snakescale)
[![codecov](https://codecov.io/gh/clintval/snakescale/branch/master/graph/badge.svg)](https://codecov.io/gh/clintval/snakescale)
[![Documentation Build Status](https://readthedocs.org/projects/snakescale/badge/?version=latest)](https://snakescale.readthedocs.io/en/latest/?badge=latest)
[![PyPi Release](https://badge.fury.io/py/snakescale.svg)](https://badge.fury.io/py/snakescale)
[![Anaconda-Server Badge](https://anaconda.org/clintval/snakescale/badges/version.svg)](https://anaconda.org/clintval/snakescale)
[![Python Versions](https://img.shields.io/pypi/pyversions/snakescale.svg)](https://pypi.python.org/pypi/snakescale/)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Non-strict wrappers for the data pipelining language Snakemake.

```bash
❯ pip install 'snakescale[full]'
```

Features:

- Do the wrappers in the [official wrapper repository](https://bitbucket.org/snakemake/snakemake-wrappers) get you half of the way to writing rules in only Python syntax?
- Do want your rules fully parameterized with the `input`, `output`, `resources`, and `params` keys only?
- Do you want to use the builtin Python types as values to a rule?
- Do you want to use the Snakemake resource system for JVM resources?
- Do you want a Snakemake wrapper which hard-codes as little as possible besides the **style** of the CLI it's wrapping?
Read the documentation at: [snakescale.readthedocs.io](http://snakescale.readthedocs.io/)
- Snakescale does all this and is conda environment compatible!

This project aims to wrap bioinformatics utilities with style and variable converters instead of strict, inflexible shell templates. The wrappers in this project are unaware of the command line flags of the tool the wrapper is wrapping!

## Example

```python
from snakescale import scale

rule bedtools_subtract:
    input:
        a='data/a.bed',
        b='data/b.bed'
    output: 'data/result.bed'
    params:
        no_name_check=True,
        g='data/ref.genome'
    wrapper: scale('bedtools', 'subtract')
```

Which executes this under the hood:

```bash
❯ bedtools subtract -a data/a.bed -b data/b.bed -nonamecheck -g data/ref.genome > data/result.bed
```

By invoking the following:

```bash
❯ snakemake -F --use-conda

Building DAG of jobs...

Creating conda environment .../bedtools/subtract/environment.yaml...
Downloading remote packages.
Environment for .../bedtools/subtract/environment.yaml created (location: .snakemake/conda/32f9fcde)
Using shell: /usr/local/bin/bash
Provided cores: 1

Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	bedtools_subtract
	1

[Fri Dec 28 13:13:47 2018]
rule bedtools_subtract:
    input: data/a.bed, data/b.bed
    output: data/result.bed
    jobid: 0

Activating conda environment: .snakemake/conda/32f9fcde

[Fri Dec 28 13:13:47 2018]
Finished job 0.

1 of 1 steps (100%) done
Complete log: .snakemake/log/2018-12-28T131312.471617.snakemake.log
```