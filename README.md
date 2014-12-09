# Source size calculator

## Introduction
This repository contains the source size calculation tools introduced by Lovric _et al._ [Opt. Express 22, 2745 (2014), http://doi.org/q9m].
Feel free to report any kind of bugs or suggestions for improvement as [GitHub issues here](https://github.com/gnudo/source-size-calculator/issues).

## Structure
- [`classes`](classes) contains all methods for calculating source sizes
- [`scripts`](scripts) contains both the source size calculation and the multilayer characterization implementations
- [`examples`](examples) contains examples that can be run to see how the tools work

## Tutorial
The best way to start using the tool is to run `run_all_tests.m` from the [`examples`]().
The test currently consists of 2 parts:

1. it creates a subfolder called `exp_data` with *fake* experimental data
2. it performs the fitting algorithm on this data, prints and plots the results.

## Roadmap
There are still some steps left that will follow in the near future:

- [ ] automatic uncertainty analysis
- [ ] rebase [`scripts`](scripts) and [`examples`](examples) folders to reduce redundancies
- [ ] remove unused methods and complete documentation
- [ ] remove hard-coded variables
