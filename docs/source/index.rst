.. ASdeep documentation master file, created by
   sphinx-quickstart on Mon Nov  8 11:36:29 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======
ASdeep
======
``ASdeep`` is designed to attribute genetic variants to allelic imbalance.

.. image:: overview.png
   :alt:

============
Introduction
============

The allelic imbalance of the heritable architecture of human genome, which is introduced by sexual reproduction, results in phenotypic variantion within and among populations.
Many studies have suggested that it plays roles in complex human phenotypes by differentially regulating genes expression.
In routine genomic diagnosis researchers usually start with genetic variants (e.g., single nucleotide variants, SNV) and infer the causal variants to explain the abnomality.
However, the practice exclusively takes genetic variants in coding regions, which usually have large effect-size, into account but not the regulation affected by non-coding variants, which are missing pieces of the puzzle.
Therefore, to fill the missing part, a handy tool is required to study the regulation potential of non-coding variants.


---------------------
Counting allelic read
---------------------

Omics profiling methods can measure millions of molecules time- and cost-efficiently.
Large-scale profiling of allelic imbalance is usually accomplished using RNA sequencing technology.
We quantified allelic read counts from RNA-seq results including xxx samples in total.
In brief, the quality control and preprocessing for raw reads in :token:`FASTQ` format were done by ``fastp``.
Then, the clean reads were aligned to reference genome (GRCh37) by ``STAR``.
To remove mapping bias to reference allele, ``WASP`` algorithm was applied using individual genotypes as known genetic variants (i.e., SNPs).
From the obtained :token:`BAM` files, ``GATK/ASEReadCounter`` were used to count allelic read.


--------------------
Infering ASE effects
--------------------

After obtaining allelic read counts, we firstly identified parental allele read counts using the phasing information in the genotype data set.
Transcripts
Then, supposing that the read counts from one parent are respect to Beta-Binomial distribution, the P of success and the corresponding high density intervals were estimated by Markov chain Monte-Carlo (MCMC) methods.
The inference were implemented into the sub-command ``inferai`` with a few commandline options to tweak the hyper-parameters.

.. tip::
   Tip



------------------------------
Model inputs and Hilbert curve
------------------------------

Theoretically, one organism's genome contains all heritable information that controls its biological structures, developments, and behaviours.


--------------------------------------------
Attribution methods for model interpretation
--------------------------------------------

The complex of NN is one of the main barriers that hang back its biological application.
Thanks to the advances of attribution methods (AM) to dissect the NN.
``ASdeep`` incorporates recently published AM, including xxx, xxx, and xxx.


========================
Dependencies and install
========================

------------
Dependencies
------------

:command:`ASdeep` was built upon a few well-designed frameworks, for instance
the widely-used deep-learning framework ``PyTorch``.

.. code-block:: text

   PyTorch
   PyMC3
   pysam
   h5py


------------
Installation
------------

:command:`ASdeep` is available at PyPi (here).

.. code:: bash

   $> python -m venv .env
   $> source .env/bin/activate
   (.env) $> pip install asdeep


====
Help
====

To show the help from CLI, use the following command:

.. code:: bash

  $> asdeep --help

The following outputs will be diplayed on the console:

.. code-block:: text

  usage: asdeep [-h] [-l LOGFILE] [-v] [-d] {inferai,makedb,train,predict} ...

  A tool to interpret variants function by deep-learning

  positional arguments:
    {inferai,makedb,train,predict}
      inferai             Infer allelic imbalance using allelic read counts.
      makedb              Make a HDF5 database for train, validation, and
                          prediction.
      train               Train the model on quantified ASE effects
      predict             Predict based on given net work state

  optional arguments:
    -h, --help            show this help message and exit
    -l LOGFILE            The file into which logging write.
    -v                    Verbose level. A counting keyword argument.
    -d, --debug           Enter debug mode using pdb.


``ASdeep`` has four sub-commands, for each sub-command, its short help text is
available by flagging ``-h/--help``:

.. code:: bash

   $> asdeep sub-command --help


=====
Usage
=====

.. tip::
   This is a test tip.


====
FAQs
====

.. DANGER::
   Please do not load models you don't trust.
