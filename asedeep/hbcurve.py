#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''A simple implementation for Hilbert curve.
'''
import logging
import itertools as it

import numpy as np
import matplotlib.pyplot as plt

from .zutils import B2M

class HilbertCurve:
    '''Hilbert curve'''
    def __init__(self, seq, kmers=4, biallele=True, seqorder='sm', mask_homo=True, mask_flank=25):
        super(HilbertCurve, self).__init__()

        self.seq = seq               # Input sequence
        self.kmers = kmers           # K-mers
        self.biallele = biallele     # Bi-allelic or monoallelic sequence?
        self.mask_homo = mask_homo   # Whether mask homogeneous sites.
        self.mask_flank = mask_flank # Whether keep a flank unmasked to the heterogeneous site
        self.seqorder = seqorder     # Way to concatenate bi-allelic sequence

        self.seqlen: int = -1
        self.order: int = -1
        self.kmers_idx_pool: list = []
        self.kmers_pool: list = []
        self.x_pool: list = []
        self.y_pool: list = []
        self.hcurve_matrix: np.ndarray = np.array([])

        self._init(biallele)

    def _init(self, biallele=True):
        self.seqlen = len(self.seq)

        if self.kmers > self.seqlen:
            raise ValueError('Sequence length should be greater equal to k')

        if biallele:
            self.seq = self._bcode2mcode()

        self.kmers_pool = [''.join(x) for x in it.product('ACGT', repeat=self.kmers)]

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
        return [self.kmers_pool.index(seq[i:i+k]) for i in range(seqlen-k+1)]

    def _mask_homo_base(self, kmers_1, kmers_2):
        kmers_1 = np.array(kmers_1)
        kmers_2 = np.array(kmers_2)
        hete_sites = (kmers_1 != kmers_2).nonzero()[0]
        n_hete_sites = hete_sites.size
        n_kmers = len(kmers_1)

        if n_hete_sites > 0:
            mask_bins = []
            for idx, val in enumerate(hete_sites):
                win_start = max(0, val - self.mask_flank)
                win_stop = min(n_kmers, val + self.mask_flank + 1)

                if idx == 0:
                    mask_bins.append((win_start, win_stop))
                else:
                    last_win_start, last_win_stop = mask_bins[-1]
                    if win_start > last_win_stop:
                        mask_bins.append((win_start, win_stop))
                    else:
                        mask_bins.pop()
                        mask_bins.append((last_win_start, win_stop))

            n_kmers = kmers_1.shape[0]
            if mask_bins:
                masked_kmers_1 = np.ones(n_kmers) * -1
                masked_kmers_2 = np.ones(n_kmers) * -1

                for start, stop in mask_bins:
                    masked_kmers_1[start:stop] = kmers_1[start:stop]
                    masked_kmers_2[start:stop] = kmers_2[start:stop]

                return list(masked_kmers_1), list(masked_kmers_2)

        return list(kmers_1), list(kmers_2)

    def _make_kmers(self):
        if self.biallele:
            kmers_a1 = self._make_kmers_base(self.seq[0])
            kmers_a2 = self._make_kmers_base(self.seq[1])

            if self.mask_homo:
                kmers_a1, kmers_a2 = self._mask_homo_base(kmers_a1, kmers_a2)

            if self.seqorder == 'ct':
                self.kmers_idx_pool = kmers_a1 + kmers_a2
            else: # FIXME: should I care about the strand?
                self.kmers_idx_pool = kmers_a1[::-1] + kmers_a2
                if self.seqorder != 'sm':
                    logging.error('Unknown way to order biallelic seq, try default: sm')
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
            x = x0 + (xi + yi) / 2
            y = y0 + (xj + yj) / 2
            self.x_pool.append(x)
            self.y_pool.append(y)
        else:
            self._hilbert_base(x0,           y0,            yi/2,  yj/2,  xi/2,  xj/2, n-1)
            self._hilbert_base(x0+xi/2,      y0+xj/2,       xi/2,  xj/2,  yi/2,  yj/2, n-1)
            self._hilbert_base(x0+xi/2+yi/2, y0+xj/2+yj/2,  xi/2,  xj/2,  yi/2,  yj/2, n-1)
            self._hilbert_base(x0+xi/2+yi,   y0+xj/2+yj,   -yi/2, -yj/2, -xi/2, -xj/2, n-1)

    def _hilbert(self, start_point=None, x_vector=None, y_vector=None,
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
            i, j = int(i), int(j)
            self.hcurve_matrix[j, i] = base

    def _crop_blank(self):
        logging.error('Not implementated yet {}'.format(__name__))

    def seq_to_hcurve(self, crop_blank=False):
        '''Convert sequence into a Hilbert curve.
        '''
        self._make_kmers()
        self._hilbert()
        self._fit_kmers_to_hcurve()
        self._seq_to_hcurve()

        if crop_blank:
            self._crop_blank()

        return self

    def hcurve_to_img(self, output_prefix='./', cmap="summer", overlay=None, overlay_cmap='Reds',
                      overlay_alpha=0.5, scatter=True, connect=True, color=True, kmer_text=True,
                      figsize=(0, 0), fmt='pdf'):
        '''Plot the Hilbert curve.'''
        if sum(figsize) <= 0:
            _size = max(self.hcurve_matrix.shape) / 4
            figsize = (_size, _size)

        fig, ax = plt.subplots(figsize=figsize)

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
                _kmer_idx = int(self.hcurve_matrix[j, i])
                if _kmer_idx != -1:
                    text = self.kmers_pool[_kmer_idx]
                else:
                    text = 'NULL'

                ax.text(i, j, text, ha='center', va='center', fontsize='x-small', rotation=45)

        ax.set_axis_off()
        fig.set_tight_layout(True)

        save_path = '{}hilbert_curve.{}'.format(output_prefix, fmt)
        fig.savefig(save_path)

        plt.close(fig)

        return self

    def get_hcurve(self, onehot=False):
        '''Generate a numpy.ndarray for the Hilbert curve of the input sequence.'''
        if onehot:
            matrix_depth = int(4 ** self.kmers)
            matrix_length = int(2 ** self.order)
            _hcurve = np.zeros((matrix_depth, matrix_length, matrix_length))
            for base, i, j in zip(self.kmers_idx_pool, self.x_pool, self.y_pool):
                i, j = int(i), int(j)  # i is column index, j is row index
                if base != -1:
                    _hcurve[base][j, i] = 1.0

            return _hcurve

        return np.expand_dims(self.hcurve_matrix, axis=0)

    def get_kmerseq(self):
        '''Get the kmer sequence of the input sequence.'''
        return self.kmers_idx_pool


if __name__ == '__main__':
    logging.warning('This module should not be executed directly.')
