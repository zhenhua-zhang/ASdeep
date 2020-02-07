#!/bin/bash

#
## GCC: 
#
gcc_version=
module load ${gcc_version}

#
## fastp: download and compile
#
# FAQ
# fastp only relies on zlib, which is already available on most Linux-like systems. If you get an
# error like undefined reference to gzbuffer when compiling fastp, you may have to update the zlib
# (http://zlib.net)

zlib_version=zlib/1.2.11-foss-2015b
module load ${zlib_version}


#
## WASP: Python virtual environment
#
python_version=

#
## WASP: `snp2h5`
#
hdf5_version=
