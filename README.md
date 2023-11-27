##### Table of Contents  

- [LRG-Signed Documentation](#lrg-signed-documentation)
  - [Installation and Activation](#installation-and-activation)
  - [How to use](#how-to-use)
    - [Signed 2D lattices](#signed-2d-lattices)


# LRG-Signed Documentation
## Installation and Activation
1. Download the project from GitHub:
```
git clone https://github.com/pvgongora/LRG-Signed
```
2. Move to the main folder and create the `LRG-Signed` anaconda enviroment by running:
```
conda env create -f LRG-Signed_env.yml
```
3. Verify that installation was succesfull by means of `conda list`. Activate the `conda` environment before executing scripts:
```
conda activate LRG-Signed
``` 
4. For installing in [development mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) the package just move to the package directory (with `cd src` from the project directory) and then run
```
pip install --editable .
```
## How to use
### Signed 2D lattices
A full description of these scripts can be read off by executing them with the `-h` flag. Note that the script must run from the parent folder (the root folder of the project).
```
$ python3 src/signed_lattice2dsq_input.py -h
```
Basic usage is given by, e.g.
```
python3 src/signed_lattice2dsq_input.py -nA 20 10 0.1
```
where optional parameter `-nA` set the number of averages to do at 20 (1000 if not specified), positional arguments 10 and 0.1 are, respectively `L` (the side of the lattice) and `p` (the fraction of flipped edges). 