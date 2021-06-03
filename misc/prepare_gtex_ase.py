#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2021-05-08
# Updated: 2021-05-08

import argparse
import pandas as pd

DEBUG = False


def getargs():
    '''Get commandline arguments.'''
    parser = argparse.ArgumentParser(description='Prepare GTEx ASE data for asedlp/quantify.')
    parser.add_argument('-l', '--lookup-table', dest='lookup_table', required=True,
                        help='The lookup table from GTEx.')
    parser.add_argument('-a', '--ase-table', dest='ase_table', required=True,
                        help='The ASE table (tsv/csv/ssv) from GTEx.')
    parser.add_argument('-t', '--tissue', dest='tissue', default='whlbld', required=True,
                        help='The tissue to be abstract.')
    parser.add_argument('-o', '--output-table', dest='output_table', default='GTEx_ASE_table.csv',
                        help='The output file.')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug',
                        help='Debug')

    return parser.parse_args()


def liftover(item):
    '''Liftover genetic coordination of variants based on the GTEx lookup table.'''
    global DEBUG
    if DEBUG: import pdb; pdb.set_trace()

    _, _, b38_ref, b38_alt, _ = item['variant_id'].split('_')
    b37_variantid = item['variant_id_b37']

    # Some b38 variants do not have corresponding b37 version.
    if isinstance(b37_variantid, str):
        b37_chrom, b37_pos, b37_ref, b37_alt, _ = b37_variantid.split('_')
    else:
        return None

    # Refuse variants of incorrect ref / alt alleles between two builds (b37 and b38).
    if b38_ref != b37_ref or b38_alt != b37_alt:
        return None

    item['CHR'] = b37_chrom
    item['POS'] = b37_pos

    return item


def main():
    '''The main entry of the script.'''
    global DEBUG
    if DEBUG: import pdb; pdb.set_trace()  # Enter debug if -d/--debug is on

    args = getargs()
    # Inputs
    lookup_table = args.lookup_table
    ase_table = args.ase_table

    # Behaviors
    tissue = args.tissue
    DEBUG = args.debug

    # Output
    output_table = args.output_table

    nrows = 10000 if DEBUG else None  # Only read 10000 rows if DEBUG is on.
    lookup_dataframe = pd.read_csv(lookup_table, sep='\t', nrows=nrows)
    ase_dataframe = pd.read_csv(ase_table, sep='\t', nrows=nrows)

    kept_cols = ['CHR', 'POS', 'REF_ALLELE', 'ALT_ALLELE', 'REF_COUNT', 'ALT_COUNT',
                 'TOTAL_COUNT', 'TISSUE_ID', 'variant_id_b37', 'variant_id']

    output_cols = ['CHR', 'POS', 'variant_id_b37', 'REF_ALLELE', 'ALT_ALLELE', 'REF_COUNT',
                   'ALT_COUNT', 'TOTAL_COUNT']

    new_name = {'CHR': 'contig', 'POS': 'position', 'variant_id_b37': 'variantID',
                'REF_ALLELE': 'refAllele', 'ALT_ALLELE': 'altAllele', 'REF_COUNT': 'refCount',
                'ALT_COUNT': 'altCount', 'TOTAL_COUNT': 'totalCount'}

    # Using Pandas to process the data stream.
    _ = (pd.merge(lookup_dataframe, ase_dataframe, left_on='variant_id',
                  right_on='VARIANT_ID')
         .loc[:, kept_cols]
         .query('TISSUE_ID == "{}"'.format(tissue))
         .apply(liftover, axis='columns', result_type='broadcast')
         .loc[:, output_cols]
         .dropna(axis=0)
         .rename(columns=new_name)
         .to_csv(output_table, index=False))


if __name__ == '__main__':
    main()
