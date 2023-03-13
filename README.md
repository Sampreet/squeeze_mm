# Mechanical Squeezing in Quadratically-coupled Optomechanical Systems

[![Version](https://img.shields.io/badge/manuscript-v2.2-red?style=for-the-badge)](#)
[![Version](https://img.shields.io/badge/qom-v0.9.0-red?style=for-the-badge)](#)

> A collection of all data and scripts for the work.

Author | Affiliation
------------ | -------------
[Priyankar Banerjee](https://www.iitg.ac.in/stud/p.banerjee/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India
[Sampreet Kalita](https://www.iitg.ac.in/stud/sampreet/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India
[Amarendra Kumar Sarma](https://www.iitg.ac.in/aksarma/) | Department of Physics, Indian Institute of Technology Guwahati, Assam 781039, India

## About the Work

We demonstrate the dissipation-enabled generation of a strong mechanical squeezing in a cavity optomechanical system by periodically modulating the amplitude of a single-tone laser driving the system. 
The Bogoliubov mode of the quadratically-coupled mechanical oscillator cools down to its ground state due to optomechanical sideband cooling, which contributes to a strong squeezing exceeding the $3$ dB standard quantum limit.
This sideband cooling mechanism is further optimized by numerically maximizing the ratio of the coupling sidebands.
Then we look at the crucial role of the cavity mode dissipation in inducing an enhanced squeezing.
We also verify our results with the analytical solution (under adiabatic approximation) and the exact numerical solution.
Compared with previous setups, the quadratic coupling between the mechanical oscillator and the optical mode gives rise to a robust mechanical squeezing and a strong optomechanical entanglement even for a large thermal occupancy of the mechanical mode.

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