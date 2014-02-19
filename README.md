# Source size calculator
======================
## Introduction

The following package contains the source size calculation tools introduced by Lovric _et al._ [Opt. Express 22, 2745 (2014), http://doi.org/q9m].
Feel free to report any kind of bugs or suggestions for improvement: [here](https://github.com/gnudo/source-size-calculator/issues)

## Structure

- `classes` folder: contains all methods for calculating source sizes
- `scripts` folder: contains both the source size calculation and the multilayer characterization implementations
- `examples` folder: contains examples that can be run to see how the tools work

## Tutorial

Best way to start using the tool is to run `run_all_tests.m` from `examples` folder. It currently consists of 2 parts: (1) it creates a subfolder called `exp_data` with "fake" experimental data; (2) it performs the fitting algorithm on this data, prints and plots the results.

## Roadmap

There are still some steps left that will follow in the near future:

- automatic uncertaintly analysis
- rebase `scripts` and `examples` folders to reduce redundancies
- remove unused methods and complete documentation
- remove *hard-coded* variables, etc.
