# Mechanical Squeezing in Quadratically-coupled Optomechanical Systems

[![Version](https://img.shields.io/badge/manuscript-v1.4-red?style=for-the-badge)](#)
[![Version](https://img.shields.io/badge/qom-v0.9.0-red?style=for-the-badge)](#)

> A collection of all data and scripts for the work.

Author | Affiliation
------------ | -------------
[Priyankar Banerjee](https://www.iitg.ac.in/stud/p.banerjee/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India
[Sampreet Kalita](https://www.iitg.ac.in/stud/sampreet/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India
[Amarendra Kumar Sarma](https://www.iitg.ac.in/aksarma/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India

## About the Work

We demonstrate the generation of a strong mechanical squeezing in a dissipative optomechanical system by introducing a periodic modulation in the amplitude of a single-tone laser driving the system.
The mechanical oscillator is quadratically coupled to the optical mode, which contributes to a strong squeezing exceeding the $3$ dB standard quantum limit.
The Bogoliubov mode of the mechanical oscillator also cools down to its ground state due to sideband cooling.
We further optimize this ratio of sideband strengths to introduce enhanced squeezing.
We also compare our results with the analytical (under adiabatic approximation) and the exact numerical solution.
Even for a thermal occupancy of $10^{4}$ phonons, mechanical squeezing beyond $3$ dB and a strong optomechanical entanglement is observed.

## Structure of the Repository

```
ROOT_DIR/
|
├───data/
│   ├───foo-bar/
│   │   ├───baz_xyz.npz
│   │   └───...
│   └───...
|
│───scripts/
│   ├───foo-bar/
│   │   ├───baz.py
│   │   └───...
│   └───...
│
├───.gitignore
├───CHANGELOG.md
└───README.md
```

## Execution

### Installing Dependencies

The project requires `Python 3.8+` installed via the [Anaconda distribution](https://www.anaconda.com/products/individual). 
An extensive guide to set up your python environment same can be found [here](https://sampreet.github.io/python-for-physicists/modules/m01-getting-started/m01t01-setting-up-python.html).

Once the installation is complete and `conda` is configured, it is preferable to create a new conda environment (say `qom`) and activate it using:

```bash
conda create -n qom python=3
conda activate qom
```

This project uses [The Quantum Optomechanics Toolbox](https://github.com/Sampreet/qom) via Python Package Index using `pip`:

```bash
pip install -i https://test.pypi.org/simple/ qom
```

Alternatively, [clone the repository](https://github.com/Sampreet/qom) or [download the sources](https://github.com/Sampreet/qom/archive/refs/heads/master.zip) as `.zip` and extract the contents.
Now, execute the following from *outside* the top-level directory, `ROOT_DIR`, inside which `setup.py` is located:

```bash
pip install -e ROOT_DIR
```

### Running the Scripts

To run the scripts, navigate *inside* the top-level directory, `ROOT_DIR`, and execute:

```bash
python scripts/foo-bar/baz.py
```

Here, `foo-bar` is the name of the folder and  `baz.py` is the name of the file.