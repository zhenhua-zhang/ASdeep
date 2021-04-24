#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''A simple implementation for Hilbert curve.
'''
import logging
import itertools as it

import numpy as np
import matplotlib.pyplot as plt

from .zutils import B2M

class HelbertCurve:
    def __init__(self, seq, kmers=4, biallele=True, seqorder='sm'):
        super(HelbertCurve, self).__init__()

        self.seq = seq              # Input sequence
        self.kmers = kmers          # K-mers
        self.biallele = biallele    # Biallelic or monoallelic seqeunce?

        self.seqorder = seqorder    # Way to concatate biallelic sequence
        self.seqlen = -1            # Class wide variable raw sequence length

        self.order = -1             # The order for Helbert curve (HC)
        self.x_pool = []            # Y coord for each mer, paired with self.x_pool
        self.y_pool = []            # X coord for each mer, paired with self.y_pool
        self.hcurve_matrix = []     # The elements are index of kmers by one-hot encoding

        self.kmers_idx_pool = []    # Index for each mer in HC
        self.onehot_enc_list = None # Onehot encoded list

        self._init(biallele)

    def _init(self, biallele=True):
        self.seqlen = len(self.seq)

        if self.kmers > self.seqlen:
            raise ValueError('Sequence length should be greater equal to k')

        if biallele:
            self.seq = self._bcode2mcode()

        self.onehot_enc_list = [''.join(x) for x in it.product('ACGT', repeat=self.kmers)]

    def _bcode2mcode(self): # Convert biallelic seq into two monoallelic seq
        if isinstance(self.seq, (list, tuple)) and len(self.seq) == 2:
            return self.seq
        elif isinstance(self.seq, str):
            biseq = ''.join([B2M[base] for base in self.seq])
            return biseq[0::2], biseq[1::2]
        else:
            raise TypeError('Unsupport seq format: {}'.format(type(self.seq)))

    def _make_kmers_base(self, seq):
        k = self.kmers
        seqlen = len(seq)

        if seqlen < k:
            raise ValueError('Sequence length should be greater equal to k')

        # FIXME: mer could be missing in the mer pool if there is N in sequence.
        return [self.onehot_enc_list.index(seq[i:i+k]) for i in range(seqlen-k+1)]

    def _make_kmers(self):
        if self.biallele:
            kmers_allele_1 = self._make_kmers_base(self.seq[0])
            kmers_allele_2 = self._make_kmers_base(self.seq[1])

            if self.seqorder == 'ct':
                self.kmers_idx_pool = kmers_allele_1 + kmers_allele_2
            else:
                self.kmers_idx_pool = kmers_allele_1 + kmers_allele_2[::-1]
                if self.seqorder != 'sm':
                    logging.error('Unknown way to order biallelic seq,'
                                 ' try default: sm')
        else:
            self.kmers_idx_pool = self._make_kmers_base(self.seq)

    def _calc_order_for_seq(self):
        order = -1
        exp_order = 1
        kmers_idx_pool_len = len(self.kmers_idx_pool)

        while exp_order < kmers_idx_pool_len:
            order += 1
            exp_order = 4 ** order

        return order

    def _hilbert_base(self, x0, y0, xi, xj, yi, yj, n):
        '''Make n-order Hilbert curve.

        Param:
            x0, y0: Start point
            xi, xj: Vector to indicate X-vector
            yi, yj: Vector to indicate Y-vector
            n: Number of order

        Refer:
            http://www.fundza.com/algorithmic/space_filling/hilbert/basics/
        '''
        if n <= 0:
            X = x0 + (xi + yi) / 2
            Y = y0 + (xj + yj) / 2
            self.x_pool.append(X)
            self.y_pool.append(Y)
        else:
            self._hilbert_base(x0,           y0,            yi/2,  yj/2,  xi/2,  xj/2, n-1)
            self._hilbert_base(x0+xi/2,      y0+xj/2,       xi/2,  xj/2,  yi/2,  yj/2, n-1)
            self._hilbert_base(x0+xi/2+yi/2, y0+xj/2+yj/2,  xi/2,  xj/2,  yi/2,  yj/2, n-1)
            self._hilbert_base(x0+xi/2+yi,   y0+xj/2+yj,   -yi/2, -yj/2, -xi/2, -xj/2, n-1)

    def _helbert(self, start_point=None, x_vector=None, y_vector=None,
                 order=None):
        '''Generate coordinations for Hilbert curve.
        '''
        if order is None:
            order = self._calc_order_for_seq()
            self.order = order

        x0, y0 = start_point if start_point else -0.5, -0.5
        xi, xj = x_vector if x_vector else 0, 2 ** order
        yi, yj = y_vector if y_vector else 2 ** order, 0

        self._hilbert_base(x0, y0, xi, xj, yi, yj, order)

    def _fit_kmers_to_hcurve(self):
        kmers_idx_pool_len = len(self.kmers_idx_pool)
        coord_pool_len = len(self.x_pool)

        if kmers_idx_pool_len != coord_pool_len:
            filler = [-1] * int((coord_pool_len - kmers_idx_pool_len) / 2)
            self.kmers_idx_pool = filler + self.kmers_idx_pool + filler

    def _seq_to_hcurve(self):
        matrix_length = 2 ** self.order

        self.hcurve_matrix = np.zeros((matrix_length, matrix_length))
        for base, i, j in zip(self.kmers_idx_pool, self.x_pool, self.y_pool):
            self.hcurve_matrix[int(j)][int(i)] = base

    def _crop_blank(self):
        logging.error('Not implementated yet'.format(self.__str__))

    def seq_to_hcurve(self, crop_blank=False):
        '''Convert sequence into a Helbert curve.
        '''
        self._make_kmers()
        self._helbert()
        self._fit_kmers_to_hcurve()
        self._seq_to_hcurve()

        if crop_blank:
            self._crop_blank()

        return self

    def hcurve_to_img(self, output_prefix='./', cmap="summer", overlay=None, overlay_cmap='Reds',
                      overlay_alpha=0.5, scatter=True, connect=True, color=True, kmer_text=True,
                      figsize=-1, fmt='pdf'):
        '''Plot the Helbert curve.'''
        fig, ax = plt.subplots()

        if scatter:
            ax.scatter(self.x_pool, self.y_pool, c='0', marker='.', alpha=0.5)

        if connect:
            ax.plot(self.x_pool, self.y_pool, color='0', linewidth=0.5, alpha=0.5)

        if color:
            ax.imshow(self.hcurve_matrix, cmap=cmap)

            if isinstance(overlay, np.ndarray):
                ax.imshow(overlay, alpha=overlay_alpha, cmap=overlay_cmap)

        if kmer_text:
            for i, j in zip(self.x_pool, self.y_pool):
                i, j = int(i), int(j)
                if self.hcurve_matrix[j][i] != -1:
                    text = self.onehot_enc_list[int(self.hcurve_matrix[j][i])]
                else:
                    text = 'NULL'
                ax.text(i, j, text, ha='center', va='center', fontsize='x-small', rotation=45)

        if figsize < 0:
            figsize = max(self.hcurve_matrix.shape) / 4

        ax.set_xticks([])
        ax.set_yticks([])

        fig.set_figwidth(figsize)     # Set size by k-mers
        fig.set_figheight(figsize)
        fig.set_tight_layout(True)

        save_path = '{}helbert_curve.{}'.format(output_prefix, fmt)
        fig.savefig(save_path)

        plt.close(fig)

        return self

    def get_hcurve(self, onehot=False):
        matrix_depth = int(4 ** self.kmers) if onehot else 1
        matrix_length = int(2 ** self.order)

        _hcurve = np.zeros((matrix_depth, matrix_length, matrix_length))
        for base, i, j in zip(self.kmers_idx_pool, self.x_pool, self.y_pool):
            idx, value = [base, 1.0] if onehot else [0, base]
            if base != -1:
                _hcurve[idx][int(j)][int(i)] = value

        return _hcurve


if __name__ == '__main__':
    logging.warning('This module should not be executed directly.')