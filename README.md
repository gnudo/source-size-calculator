# Source size calculator
======================
## Introduction

The following scripts include the source size calculation tools introduced by Lovric _et al._ [Opt. Express 22, 2745 (2014), http://doi.org/q9m].
Feel free to report anz kind of bugs or suggestions: [\[Issues\]](https://github.com/gnudo/source-size-calculator/issues)

## Structure

- `classes` folder: contains all methods for calculating source sizes
- `scripts` folder: contains the implementation for the source size calculation as well for the multilayer characterization
- `examples` folder: contains examples that can be run to see how the tools work

## Tutorial

Best way to start using the tool is to run `run_all_tests.m` from `examples` folder. I currently consists of 2 parts: (1) it creates a subfolder called `exp_data` with "fake" experimental data; (2) it performs the fitting algorithm on this data, prints and plots the results.

## Roadmap

There are still some steps left to make the tools fully user-friendly and usable. Following TODO-s will be resolved in the near future:

- automatic uncertaintly analysis
- rebase `scripts` and `examples` folders to reduce redundancies
- remove unused methods and complete documentation
- remove *hard-coaded* variables etc.