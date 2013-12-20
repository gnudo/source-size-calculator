# Source size calculator
======================
## Introduction

The following tools provide the source size calculation introduced by G. Lovric _et al._ (submitted to Opt. Exp.). Please be aware that the current version is in **beta state** and not (yet) very user friendly. However we are currently refactoring the code to facilitate simple usage.

## Structure

- `classes` folder: contains all methods for calculating source sizes
- `scripts` folder: will contain all higher-level methods (for fitting etc.) that can be called directly with the needed parameters (all files in this folder are currently being rewritten)
- `examples` folder: contains examples that can be run to see how the tools work (to be populated further)

## Tutorial

Best way to start using the tool is to run `run_all_tests.m` from `examples` folder. I currently consists of 2 parts: (1) it creates a subfolder called `exp_data` with "fake" experimental data; (2) it performs the fitting algorithm on this data, prints and plots the results.

## Roadmap

There are still some steps left to make the tools fully user-friendly and usable. A tentative list includes:

- consistency in variable names etc. according to paper
- automatic uncertaintly analysis
- rebase `scripts` and `examples` folders to reduce redundancies
- remove unused methods and complete documentation
- remove *hard-coaded* variables etc.
- ...