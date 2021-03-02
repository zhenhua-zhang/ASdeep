#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import gzip
import logging
from argparse import ArgumentParser

def check_io_opts(args):
    '''Check input files and output dir.'''
    output_dir = args.output_dir
    input_files = args.input_files

    if os.path.isdir(output_dir):
        logging.warning(output_dir + ' FOUND.')
    else:
        logging.info(output_dir + ' NOT FOUND! Create it.')
        os.makedirs(output_dir)

    for fpath in input_files[:]:
        if not (os.path.exists(fpath) and os.access(fpath, os.R_OK)):
            logging.error(fpath + ' NOT FOUND! Please ensure the file exists and you have access.')
            input_files.remove(fpath)

    if len(input_files):
        logging.info('In total, ' + str(len(input_files)) + ' input files available.')
    else:
        logging.error('No input available! Exit...')
        os.sys.exit()

def cond_decode(_s):
    '''The compressed line is bytes'''
    if isinstance(_s, bytes):
        return _s.decode()

    if isinstance(_s, str):
        return _s

    raise TypeError("The _s is neither a string nor a byte.")

def fetch_write_seq(geneid, fpath, optdir, beta=5e-6):
    if fpath.endswith('.gz'):
        myopen = gzip.open
    else:
        myopen = open

    p_val = None
    sample_id = None
    direction = None

    with myopen(fpath, 'r') as fhandle:
        _line = ''
        sequence = ''
        while _line != 'EOF':
            _line = cond_decode(next(fhandle, 'EOF'))
            if _line.startswith(">") and geneid in _line:
                sample_id, gene_id, p_val, direction = _line.split('|')
                sequence += _line

                _line = cond_decode(next(fhandle, 'EOF'))
                while not _line.startswith('>'):
                    sequence += _line
                    _line = cond_decode(next(fhandle, 'EOF'))
                break

    if sample_id is None or p_val is None or direction is None:
        logging.warning("Did not fetch sequence for " + geneid + " from " + fpath)
    else:
        direction = int(direction) + 1
        if float(p_val) > beta:
            direction = 1

        optdir = os.path.join(optdir, str(direction))
        if not os.path.isdir(optdir):
            logging.info(optdir + ' NOT FOUND! Create it.')
            os.makedirs(optdir)

        ofpath = os.path.join(optdir, str(sample_id.replace('>', '')) + '-' + str(geneid) + '.fa')
        with open(ofpath, 'w') as fh:
            fh.write(sequence)

def main():
    # Set up a logger
    logging.basicConfig(format='{levelname: ^8}| {asctime} | {name} | {message}', style='{', level=logging.INFO)

    parser = ArgumentParser('A utility script to split samples by ASE effects')
    # Input options
    parser.add_argument('-g', '--gene-id', dest='gene_id', required=True, help='The gene ID to split the sequence. Required.')
    parser.add_argument('-i', '--input-files', dest='input_files', nargs='*', required=True, help='Input files. Required')
    # Output options
    parser.add_argument('-a', '--alpha', dest='alpha', default=5e-6, help='The threhold p-value of ASE gene. Default: %(default)s')
    parser.add_argument('-o', '--output-dir', dest='output_dir', default='./output_dir', help='The output directory. Default: %(default)s')
    # Parse options
    args = parser.parse_args()

    # Check the I/O options
    check_io_opts(args)

    # Options
    gene_id = args.gene_id
    input_files = args.input_files
    alpha = args.alpha
    output_dir = args.output_dir

    # Fetch sequence and write to disk per effect: ASE-to-A(0), ASE(1), ASE-to-B(2)
    for ipfile in input_files:
        fetch_write_seq(gene_id, ipfile, output_dir, alpha)

if __name__ == '__main__':
    main()
