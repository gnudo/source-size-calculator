# Source size calculator
======================
## Introduction

The following package contains the source size calculation tools introduced by Lovric _et al._ [Opt. Express 22, 2745 (2014), http://doi.org/q9m].
Feel free to report any kind of [bugs or suggestions for improvement](https://github.com/gnudo/source-size-calculator/issues) or contact the authors for further questions.

## Structure

The tools consist of the following 3 layers (folders), listed from bottom to top:

- `classes`: contains all methods for calculating source sizes, interval nesting etc.
- `functions`: contains implementations for both the source size calculation and the multilayer characterization as well as the calculations for energy and source size uncertainties. They are implemented as Matlab functions. 
- `examples`: contains examples that can be run to see how the tools work. All methods are called from scripts located in this folder and all constants/parameters are set in scripts within this folder.

## Tutorial

Best way to start using the tools is to run `demo.m` from root folder. It consists of 4 parts and basically includes all steps that are presented in the paper:

1. Creation of a subfolder called `exp_data` with "fake" experimental data.
2. Run fitting algorithm on these data and save all results in a file called "results.mat"
3. Calculate uncertainties for the energy and the source sizes
4. Print and plot all the results.

## Roadmap

With version 1.1.0, all functions have been implemented that were originally planned. As we are continuing to work with the tools, small fixes etc. will be applied when necessary.

Remaing functions to be added in an eventual future 1.2.0 release are the following:
- determine magnification from experimental data automatically.