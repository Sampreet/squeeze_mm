# Robust Mechanical Squeezing beyond 3 dB in a Quadratically Coupled Optomechanical System

[![Manuscript Version](https://img.shields.io/badge/version-2.2-red?style=for-the-badge)](https://doi.org/10.1364/JOSAB.483944)
[![Toolbox Version](https://img.shields.io/badge/qom-v1.0.1-red?style=for-the-badge)](https://sampreet.github.io/qom-docs/v1.0.1)
[![Last Updated](https://img.shields.io/github/last-commit/sampreet/squeeze_mm?style=for-the-badge)](https://github.com/sampreet/squeeze_mm/blob/master/CHANGELOG.md)

> J. Opt. Soc. Am. B 40, 1398 (2023)

Author | Affiliation
------------ | -------------
[Priyankar Banerjee](https://www.iitg.ac.in/stud/p.banerjee/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India
[Sampreet Kalita](https://www.iitg.ac.in/stud/sampreet/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India
[Amarendra Kumar Sarma](https://www.iitg.ac.in/aksarma/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India

Contributing Part | PB | SK
------------ | ------------ | -------------
Literature review | 70% | 30%
Idea and formulation | 50% | 50%
Derivations of expressions | 80% | 20%
Parameter sweeping | 30% | 70%
Illustrations and plots | 30% | 70%
Results and discussion | 60% | 40%
Manuscript preparation | 60% | 40%


## About the Work

We demonstrate the dissipation-enabled generation of a strong mechanical squeezing in a cavity optomechanical system by periodically modulating the amplitude of a single-tone laser driving the system. 
The Bogoliubov mode of the quadratically-coupled mechanical oscillator cools down to its ground state due to optomechanical sideband cooling, which contributes to a strong squeezing exceeding the $3$ dB standard quantum limit.
This sideband cooling mechanism is further optimized by numerically maximizing the ratio of the coupling sidebands.
Then we look at the crucial role of the cavity mode dissipation in inducing an enhanced squeezing.
We also verify our results with the analytical solution (under adiabatic approximation) and the exact numerical solution.
Compared with previous setups, the quadratic coupling between the mechanical oscillator and the optical mode gives rise to a robust mechanical squeezing and a strong optomechanical entanglement even for a large thermal occupancy of the mechanical mode.

## Notebooks

* [Plots in the Latest Version of the Manuscript](notebooks/v2.2_qom-v1.0.1/plots.ipynb)

## Structure of the Repository

```
ROOT_DIR/
|
├───data/
│   ├───bar/
│   │   ├───baz_xyz.npz
│   │   └───...
│   └───...
|
├───notebooks/
│   ├───bar/
│   │   ├───baz.ipynb
│   │   └───...
│   │
│   ├───foo_baz.ipynb
│   └───...
|
│───scripts/
│   ├───bar/
│   │   ├───baz.py
│   │   └───...
│   └───...
|
├───systems/
│   ├───__init__.py
│   ├───Foo.py
│   └───...
│
├───.gitignore
├───CHANGELOG.md
└───README.md
```

Here, `foo` represents the module or system and `bar` represents the version.

## Installing Dependencies

All numerical data and plots are obtained using the [Quantum Optomechanics Toolbox](https://github.com/sampreet/qom), an open-source Python framework to simulate optomechanical systems.
Refer to the [QOM toolbox documentation](https://sampreet.github.io/qom-docs/v1.0.1) for the steps to install this libary.

## Running the Scripts

To run the scripts, navigate *inside* the top-level directory, and execute:

```bash
python scripts/bar/baz.py
```

Here, `bar` is the name of the folder (containing the version information) inside `scripts` and `baz.py` is the name of the script (refer to the repository structure).