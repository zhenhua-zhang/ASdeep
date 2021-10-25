# Author : Zhenhua Zhang
# Email  : zhenhua.zhang217@gmail.com
# Created: Mar 12, 2020
# Updated: Oct 06, 2021
"""Database"""

import csv
import copy
import traceback
from argparse import Namespace
from collections import OrderedDict

import h5py as h5
import numpy as np
from hilbert import decode
from pysam.libcbcf import VariantFile
from pysam.libctabix import TabixFile
from pysam.libcfaidx import FastaFile

try:
    from .zutils import B2M
    from .zutils import M2B
    from .zutils import calc_bits
    from .zutils import LogManager
    from .zutils import make_all_mers
except ImportError as e:
    from zutils import B2M
    from zutils import M2B
    from zutils import calc_bits
    from zutils import LogManager
    from zutils import make_all_mers


class AllelicDiffLenError(Exception):
    pass


class HbmatrixNotSquareErr(Exception):
    pass


class HDF5Database(h5.File):
    """The HDF5 database for Hilbert Curve."""
    def __init__(self, matrix: np.ndarray = None, key: str = None,
                 attrs: dict = None, dbpath: str = "./db.h5",
                 logman: LogManager = LogManager("HDF5Database"), **kwargs):
        """Create a new HDF5 database to store the Hilbert matrix."""
        if kwargs.get("mode", None) is None:
            kwargs["mode"] = "a"

        super(HDF5Database, self).__init__(dbpath, **kwargs)

        self._logman = logman
        if (isinstance(matrix, np.ndarray) and isinstance(key, str)
                and isinstance(attrs, dict)):
            self.add_matrix(matrix, key, attrs)

    @property
    def dbpath(self):
        return self.filename

    def add_matrix(self, matrix: np.ndarray, key: str, attrs: dict):
        """Add a matrix to the the database."""
        _mg = self.get(f"/{key}")
        if _mg is None:
            self.create_dataset(key, matrix.shape, "i", matrix)
        else:
            self._logman.warning(f"{key} exists, use update_matrix instead.")
        self.update_matrix_attrs(key, attrs)

    def del_matrix(self, key):
        """Delete a matrix by the key."""
        if self.get(f"/{key}", None) is None:
            self._logman.error(f"No {key} was found in current database.")
        else:
            try:
                del self[f"/{key}"]
            except Exception as e:
                self._logman.error(e)

    def update_matrix(self, matrix: np.ndarray, key: str, attrs: dict = {},
                      del_old: bool = False):
        """Update a matrix after removing the old one if already exists."""
        if del_old:
            self.del_matrix(key)
        self.add_matrix(matrix, key, attrs)

    def update_matrix_attrs(self, key: str, attrs: dict):
        """Update the attributes of a matrix by the key."""
        for attr, val in attrs.items():
            self.get(f"/{key}").attrs.create(attr, val)

    def get_matrix(self, key):
        """Return the Hilbert curve matrix for given key."""
        return self.get(f"/{key}")

    def get_attrs(self, key):
        """Return attributes of the Hilbert curve matrix under given key."""
        return dict(self.get(f"/{key}").attrs)


class HilbertCurve:
    """Hilbert curve.

    NOTE:
        `mask_homo` only works for string DNA sequence for now (Oct 11, 2021).
    """
    def __init__(self, source=None, kmer=4, mask_homo=False, flank=25, pads=-1,
                 logman: LogManager = LogManager("HilbertCurve")):
        self._logman = logman
        self._mer2idx, self._idx2mer = make_all_mers(kmer)

        if isinstance(source, str):
            self.from_sequence(source, kmer, mask_homo, flank, pads)
        elif isinstance(source, np.ndarray):
            self.from_hbcmat(source)
        else:
            self._bits = 0
            self._hbcurve = np.array([])
            self._sequence = ""
            self._kmered_biseq = np.array([])

    def __repr__(self):
        trunc_seq = self.sequence[:10]
        trunc_seq = trunc_seq + " ..." if len(trunc_seq) else ""
        hbc_size = self.hbcurve.shape
        hbc_bits = self.bits
        return f"HilbertCurve:\n" + \
               f"Sequence: {trunc_seq}\n" + \
               f"HBC size: {hbc_size}\n" + \
               f"HBC bits: {hbc_bits}\n"

    def _mkmers(self, sequence, kmer):
        return np.array([
            self._mer2idx.get(sequence[i:i + kmer], 0)
            for i in range(len(sequence) - kmer + 1)
        ])

    def from_sequence(self, seq: str, kmer: int = 4, mask_homo: bool = False,
                      flank: int = 25, pads: int = -1):
        self._sequence = copy.deepcopy(seq)

        biseq = "".join([B2M[base] for base in seq])
        a1_seq, a2_seq = biseq[0::2], biseq[1::2]
        self._a1_idx = self._mkmers(a1_seq, kmer)
        self._a2_idx = self._mkmers(a2_seq, kmer)

        biseq_len = len(self._a1_idx) + len(self._a2_idx)
        self._bits = calc_bits(biseq_len)

        if mask_homo:
            if self._a1_idx.size != self._a2_idx.size:
                raise AllelicDiffLenError("Alleles should have equal length!")

            het_sites = (self._a1_idx != self._a2_idx).nonzero()[0]
            if het_sites.size:
                max_idx = self._a1_idx.size
                masks = [
                    np.arange(max(x - flank, 0), min(x + flank, max_idx))
                    for x in het_sites ]
                masks = np.concatenate(masks, axis=None)
                trans = np.zeros_like(self._a1_idx)
                trans[masks] = 1
                self._a1_idx, self._a2_idx = self._a1_idx * trans, self._a2_idx * trans
            else:
                self._logman.info("No heterozygous sites was found")

        hcwidth = 2 ** self._bits
        min_hbc_size = int(hcwidth ** 2)
        if biseq_len != min_hbc_size:
            pads_array = np.zeros(int((min_hbc_size - biseq_len)/2)) + pads
            self._a1_idx = np.append(self._a1_idx, pads_array)
            self._a2_idx = np.append(self._a2_idx, pads_array)

        self._kmered_biseq = np.append(self._a1_idx[::-1], self._a2_idx)

        self._locus = decode(np.arange(min_hbc_size), 2, self._bits)
        self._hbcurve = np.zeros(min_hbc_size).reshape((1, hcwidth, hcwidth))
        for idx in np.arange(min_hbc_size):
            x, y = self._locus[idx]
            self._hbcurve[0, y, x] = self._kmered_biseq[idx]

        return self

    def from_hbcmat(self, hbmatrix: np.ndarray, make_seq=True):
        _, hcwidth, hcheight = hbmatrix.shape
        if hcwidth != hcheight:
            raise HbmatrixNotSquareErr("The HilbertCurve should be a square.")

        min_hbc_size = hcwidth * hcheight
        self._bits = calc_bits(min_hbc_size)
        self._locus = decode(np.arange(min_hbc_size), 2, self._bits)

        self._hbcurve = copy.deepcopy(hbmatrix)
        self._kmered_biseq = np.zeros(min_hbc_size)
        for i, (x, y) in enumerate(self._locus):
            self._kmered_biseq[i] = self._hbcurve[0, y, x]

        n_kmers = int(min_hbc_size/2)
        self._a1_idx = self._kmered_biseq[:n_kmers][::-1]
        self._a2_idx = self._kmered_biseq[n_kmers:]

        self._sequence = ""
        if make_seq:
            for p_a1_idx, p_a2_idx in zip(self._a1_idx, self._a2_idx):
                if p_a1_idx in self._idx2mer and p_a2_idx in self._idx2mer:
                    p_a1_mer = self._idx2mer[p_a1_idx]
                    p_a2_mer = self._idx2mer[p_a2_idx]
                    self._sequence += M2B[p_a1_mer[0] + p_a2_mer[0]]
                else:
                    for a1_base, a2_base in zip(p_a1_mer[1:], p_a2_mer[1:]):
                        self._sequence += M2B[a1_base + a2_base]
                    break

        return self
    
    def subset(self, length: int = 1024, strand = 1, **kwargs):
        if strand in [-1, "-1", "-"]:
            return HilbertCurve(self._sequence[:length], **kwargs)

        return HilbertCurve(self._sequence[length:], **kwargs)

    @property
    def bits(self):
        return self._bits

    @property
    def hbcurve(self):
        return self._hbcurve

    @property
    def kmered_seq(self):
        return self._kmered_biseq

    @property
    def sequence(self):
        return self._sequence

    @property
    def allelic_attrs(self):
        return self._a1_idx, self._a2_idx


class BEDDict(csv.DictReader):
    def __init__(self, csvfile, idx_col=0, row_key_tran=None, **kwargs):
        self._row_key_tran = row_key_tran

        if csvfile is None:
            self._csv_handler = None
            self._metadict = None
            self._colkeys = None
            self._rowkeys = None
        else:
            self._csv_handler = open(csvfile, "r")
            super(BEDDict, self).__init__(self._csv_handler, **kwargs)
            self._mk_metadict(idx_col)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)

        if hasattr(self._csv_handler, "close"):
            self._csv_handler.close()

    def __getitem__(self, key):
        if self._metadict:
            if key in self._metadict:
                return self._metadict[key]
        return {}

    def __contains__(self, key):
        return key in self._metadict

    def __len__(self):
        return 0 if self._rowkeys is None else len(self._rowkeys)

    def _mk_metadict(self, idx_col, discard_idx_col=False):
        _metadict = OrderedDict()
        _colkeys, _rowkeys = [], []
        for idx, record in enumerate(self):
            row_key, col_key = idx, "Index"
            if isinstance(idx_col, int):
                if idx_col >= 0:
                    col_key, row_key = list(record.items())[idx_col]
            elif isinstance(idx_col, (list, tuple)):
                rc_key = [r for i, r in enumerate(record.items())
                          if i in idx_col]
                col_key = tuple([key[0] for key in rc_key])
                row_key = tuple([key[1] for key in rc_key])

                if isinstance(self._row_key_tran, (tuple, list)):
                    if len(self._row_key_tran) != len(row_key):
                        raise ValueError("row_key_tran functions should match"
                                         " row_key one by one")
                    row_key = [t(k) if callable(t) else k
                               for t, k in zip(self._row_key_tran, row_key)]
                    row_key = tuple(row_key)
            else:
                raise TypeError(f"Unsupported key type: {type(idx_col)}")

            if row_key in _metadict:
                raise KeyError("Duplicated key entry")

            if discard_idx_col:
                if isinstance(idx_col, int):
                    record.pop(col_key)
                elif isinstance(idx_col, (list, tuple)):
                    for key in col_key:
                        record.pop(key)
                else:
                    raise TypeError(f"Unsupported key type: {type(idx_col)}")

            _metadict[row_key] = record
            if idx == 0:
                _colkeys = list(record.keys())
            _rowkeys.append(row_key)

        self._metadict = _metadict
        self._colkeys = _colkeys
        self._rowkeys = _rowkeys

    @property
    def csv_handler(self):
        return self._csv_handler

    @property
    def meta_dict(self):
        return self._metadict

    @property
    def row_keys(self):
        return self._rowkeys

    @property
    def col_keys(self):
        return self._colkeys

    def close(self):
        if not self._csv_handler.closed:
            self._csv_handler.close()


def _makedb(vcf, bed, fasta, db_path, shift, meta=None, process=1):
    """Make HDF5 database.

    Make a HDF5 database from given genetic variants (SNP) for specific genome
    regions according to given reference genome
    """
    def _parse_bed(itvl, sep="\t"):
        itvls = itvl.split(sep)
        if len(itvls) < 6:
            return None
        chrom, start, end, name, _, strand, *_ = itvls
        return chrom, start, end, name, strand

    samples = None
    with VariantFile(vcf) as variants:
        samples = variants.header.samples

    with TabixFile(bed) as intervals, FastaFile(fasta) as sequence, \
            BEDDict(meta) as metadata:
        for per_itvl in intervals.fetch():
            itvl_lst = _parse_bed(per_itvl)
            if itvl_lst is None:
                continue

            chrom, start, end, name, strand = itvl_lst
            start, end = int(start), int(end)
            attr_lst = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "strand": strand
            }

            if strand in ["-", "-1", -1]:
                start, end = end, end + shift - 1
            else:
                end, start = start, (start - shift if start > shift else 1)
            region = f"{chrom}:{start}-{end}"

            with HDF5Database(dbpath=f"{db_path}/{name}.h5") as hdf5db:
                for per_sample in samples:
                    sub_seqs = sequence.fetch(region=region)

                    with VariantFile(vcf) as variants:
                        variants.subset_samples([per_sample])
                        sub_vars = variants.fetch(region=region)

                        for per_var in sub_vars:
                            a1, a2 = per_var.samples[per_sample].alleles
                            if a1 == a2:
                                continue
                            istpos = per_var.pos - start
                            sub_seqs = "".join([
                                sub_seqs[:istpos - 1], M2B[f"{a1}{a2}"],
                                sub_seqs[istpos:]
                            ])

                        meta_info = metadata[per_sample]
                        if meta_info:
                            attr_lst.update(meta_info)
                            hbcurve = HilbertCurve(sub_seqs).hbcurve
                            hdf5db.add_matrix(hbcurve, per_sample, attr_lst)


def makedb(args: Namespace, logman: LogManager = LogManager("Make DB")):
    """Make a database for training."""
    reference_genome = args.reference_genome
    genetic_variants = args.genetic_variants
    genome_intervals = args.genome_intervals
    metadata_table = args.metadata_table
    n_base_pairs = args.n_base_pairs
    output_dir = args.output_dir

    logman.info(f"Reference genome: {reference_genome}")
    logman.info(f"Genetic variants: {genetic_variants}")
    logman.info(f"Genome intervals: {genome_intervals}")
    logman.info(f"Metadata table:   {metadata_table}")
    logman.info(f"No. base pairs:   {n_base_pairs}")
    logman.info(f"Output directory: {output_dir}")

    _makedb(vcf=genetic_variants, bed=genome_intervals, fasta=reference_genome,
            db_path=output_dir, shift=n_base_pairs, meta=metadata_table)
