#!/bin/bash

#
## A short script to convert README.md to README.pdf
#

# NOTE: other formats are also possible.

set -Eeu -o pipefail

mypath=$(dirname $(realpath $0))
input_md=$mypath/../README.md
output_pdf=$mypath/../README.pdf
pandoc \
    --pdf-engine xelatex \
    --highlight-style tango \
    --output $output_pdf \
    $input_md

tput bold; tput setaf 6; echo Check $(readlink -f $output_pdf) for the PDF file!

