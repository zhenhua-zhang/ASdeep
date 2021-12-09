asdeep: Deep-learning tool to interpret variatn function by allelic imbalance.


## Introduction

![Pipeline overview](./figures/overview-of-the-asedlp-pipeline.png)


## Dependencies

The tool was written in Python 3 and depends on a set of well-developed
packages, mainly including the followings:

- numpy == 1.21.4
- h5py == 3.6.0
- pymc3 == 3.11.4
- pysam == 0.18.0
- torch == 1.10.0
- captum == 0.4.1
- tensorboard == 2.7.0
- torchvision == 0.11.1
- scikit-learn == 1.0.1


## Installation

We highly recommend to use a Python virtual environment.
Here is an example to create one by the `venv` module shiped with Python 3.

```
$> python -m venv .env
$> source .env/bin/activate
```

Then, clone the repository and install it using `pip` module.

```
(.env) $> git clone https://github.com/zhenhua-zhang/asdeep
(.env) $> python -m pip install -e .
```

Or you can install it from PyPi.

```
(.env) $> python -m pip install asdeep
```

## Usage

### Example dataset

We also provide a small dataset for (potential) users to play with. The dataset
is permantly available by the [link](https://zenodo.org/xxx) at Zenodo

### Step-by-step usage example

#### Create a working space


#### Prepare input files


#### Prepare the allelic read counts


#### Estimate the ASE effects for an isoform


#### Create databases for train/test/prediction


#### Train a model and test it


#### Predict the unseen samples


### Outputs

- Model state
- Evaluation matrix
- Prediction results
- Feature attribution plot


## Known issues

1. `No section: 'blas'` problem raised by `pymc3`

The `pymc3` package depends on `Theano-PyMC` which is a special `Theano` branch.
However, one knownd issues when importing pymc3 is the `No section: 'blas'`.
if this is the issues you come across, please check if the BLAS installed
correctly in your system.

## Contact

If you have issue or question, please create a new issue or send an email to
the following address.

Zhenhua Zhang <zhenhua.zhang217@gmail.com

