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

Best way to start using the tools is to run `demo.m` from `examples` folder. It consists of 4 parts and basically includes all steps that are presented in the paper:

1. Creation of a subfolder called `exp_data` with "fake" experimental data.
2. Run fitting algorithm on this data and save all results in a file called "results.mat"
3. Calculate uncertainties for the energy and the source sizes
4. Print and plot all the results.

## Roadmap

Following TODO-s are currently being conducted:

- rebase `scripts` and `examples` folders to reduce redundancies
- remove unused methods and complete documentation
- remove *hard-coded* variables, etc.
- apply the tools without knowing the grating parameters (i.e. one has a pure phase grating etc.)