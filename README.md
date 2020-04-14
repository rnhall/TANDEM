# TANDEM
A tool for finding tandem repeats in genomic assembles and reads.

## Installation Instructions:

Install the latest version of Miniconda (or Anaconda):
```
$wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

$chmod +x Miniconda3-latest-Linux-x86_64.sh

$./Miniconda3-latest-Linux-x86_64.sh
```

Restart your terminal to auto-enter your base conda environment.

Install git if it's not installed:
```
$sudo apt-get install git
```

Clone this repository into a directory of your choice:
```
$git clone https://github.com/rnhall/TANDEM.git

$cd TANDEM
```

In the main directory, you'll see environment.yml, this tells conda which dependencies to install.

We'll make a new conda environment and set it up with environment.yml.
```
$conda config --add channels conda-forge

$conda config --add channels bioconda

$conda create --name tandem

$conda activate tandem

$conda env update -f environment.yml
```

Once the environment is done being set up, try running TANDEM on a test file:

(You may have to enable the file to be run)
```
$chmod +x tandem.py
$./tandem.py -genome test_sequence
```
You may also need to install libgtk-3-0:
```
$sudo apt install libgtk-3-0
```

## How to use TANDEM:

//IN PROGRESS



