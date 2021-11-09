# Author : Zhenhua Zhang
# Email  : zhenhua.zhang217@gmail.com
# Created: Mar 12, 2020
# Updated: Oct 06, 2021
"""Database"""

import copy
from argparse import Namespace

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
    from .tabdict import CSVDict
except ImportError as e:
    from zutils import B2M
    from zutils import M2B
    from zutils import calc_bits
    from zutils import LogManager
    from zutils import make_all_mers
    from tabdict import CSVDict


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
        1. `mask_homo` only works for string DNA sequence (Oct 11, 2021).
        2. The allelic sequences are connected at the TSS end, which take the
        strand in consideration.
    """
    def __init__(self, source=None, kmer: int = 4, strand=1, dtype=np.uint8,
                 logman: LogManager = LogManager("HilbertCurve")):
        self._logman = logman
        self._mer2idx, self._idx2mer = make_all_mers(kmer)

        self._kmer = kmer
        self._strand = strand
        self._dtype = dtype

        # Footprint tracker
        self._is_hbcmat = False
        self._is_masked = False
        self._is_subset = False

        if isinstance(source, str) and len(source) != 0:
            self._from_dnaseq(source)
        elif isinstance(source, np.ndarray) and source.size != 0:
            self._from_hbcmat(source)
        else:
            self._bits = -1
            self._hbc_size = 0
            self._dnaseq = ""
            self._hbcurve = np.array([])
            self._kmer_biseq = np.array([])

            self._a1_val = np.array([])
            self._a2_val = np.array([])

    def __repr__(self):
        trunc_seq = ""
        if len(self.dnaseq):
            trunc_seq = self._dnaseq[:10]
            trunc_seq = trunc_seq + " ..." if len(trunc_seq) else ""

        hbc_size = "None"
        if self.hbcmat.size:
            hbc_size = self.hbcmat.shape

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

    def _from_dnaseq(self, dnaseq: str, kmer=None, strand=None):
        kmer = self._kmer if kmer is None else kmer
        strand = self._strand if strand is None else strand

        self._dnaseq = copy.deepcopy(dnaseq)

        biseq = "".join([B2M[base] for base in dnaseq])
        self._a1_val = self._mkmers(biseq[0::2], kmer)
        self._a2_val = self._mkmers(biseq[1::2], kmer)

        biseq_len = len(self._a1_val) + len(self._a2_val)
        self._bits = calc_bits(biseq_len)

        self._hbc_size = int((2 ** self._bits) ** 2)
        # Allele 1 should be always on the left-hand side.
        if strand in [1, "1", "+"]:
            self._kmer_biseq = np.append(self._a1_val, self._a2_val[::-1])
        else:
            self._kmer_biseq = np.append(self._a1_val[::-1], self._a2_val)

        self._locus = decode(np.arange(self._hbc_size), 2, self._bits)

        return self

    def _from_hbcmat(self, hbcmat: np.ndarray, strand=1, dtype=None):
        """Construct an instance from existing Hilbert curve matrix"""
        if strand is None:
            strand = self._strand

        if dtype is None:
            dtype = self._dtype

        self._hbcurve = copy.deepcopy(hbcmat)

        _, hcwidth, hcheight = hbcmat.shape
        if hcwidth != hcheight:
            raise HbmatrixNotSquareErr("The HilbertCurve should be a square.")

        self._hbc_size = hbcmat.size
        self._bits = calc_bits(self._hbc_size)

        self._locus = decode(np.arange(self._hbc_size), 2, self._bits)
        self._kmer_biseq = np.zeros(self._hbc_size, dtype=dtype)
        for i, (x, y) in enumerate(self._locus):
            self._kmer_biseq[i] = hbcmat[0, y, x]

        a12_len = int(self._hbc_size / 2)
        assert a12_len * 2 == self._hbc_size, "HBcurve size should be even."

        self._a1_val = self._kmer_biseq[:a12_len]
        self._a2_val = self._kmer_biseq[a12_len:]
        if strand in [1, "1", "+"]:
            self._a2_val = self._a2_val[::-1]
        else:
            self._a1_val = self._a1_val[::-1]

        self._dnaseq = ""
        self._is_hbcmat = True # Footprint tracker

        return self

    def subset(self, length: int = 1024):
        lf_bound = int((self._hbc_size - 2 * length) / 2)
        rt_bound = lf_bound + 2 * length

        self._kmer_biseq = self._kmer_biseq[lf_bound:rt_bound]
        biseq_len = self._kmer_biseq.size # Update class-wide variables.

        # If we subset the kmer_biseq, the self._locus should be re-calculated
        self._bits = calc_bits(biseq_len)
        self._hbc_size = int((2 ** self._bits) ** 2)
        self._locus = decode(np.arange(self._hbc_size), 2, self._bits)

        if self._dnaseq:
            if self._strand in [1, "1", "+"]:
                self._dnaseq = self._dnaseq[-length:]
            else:
                self._dnaseq = self._dnaseq[:length]

        # TODO: dulicated code, try to use a method.
        a12_len = int(self._hbc_size / 2)
        assert a12_len * 2 == self._hbc_size, "HBcurve size should be even."

        self._a1_val = self._kmer_biseq[:a12_len]
        self._a2_val = self._kmer_biseq[a12_len:]
        if self._strand in [1, "1", "+"]:
            self._a2_val = self._a2_val[::-1]
        else:
            self._a1_val = self._a1_val[::-1]

        self._is_subset = True  # Footprint tracker

        return self
    
    def mask_homo(self, flank: int = 25, fills: int = -1):
        het_sites = (self._kmer_biseq == self._kmer_biseq[::-1]).nonzero()[0]
        if het_sites.size:
            masks = [ np.arange(max(x-flank, 0), min(x+flank, self._hbc_size))
                      for x in het_sites ]
            masks = np.concatenate(masks, axis=None)
            self._kmer_biseq[masks] = fills

            self._is_masked = True # Footprint tracker

        return self

    def get_hbcmat(self, fills=-1, dtype=None):
        if dtype is None:
            dtype = self._dtype

        hcwidth = int(2 ** self._bits)
        hbcurve = np.full(self._hbc_size, fills, dtype=dtype)
        hbcurve = hbcurve.reshape((1, hcwidth, hcwidth))

        start_idx = int((self._hbc_size - self._kmer_biseq.size) / 2)
        end_idx = self._hbc_size - start_idx

        idx = start_idx
        while start_idx <= idx < end_idx:
            x, y = self._locus[idx]
            hbcurve[0, y, x] = self._kmer_biseq[idx - start_idx]
            idx += 1

        return hbcurve

    def get_dnaseq(self, strand="1"):
        if self._is_masked:
            self._logman.warning("The instance was masked, no DNA sequence was"
                                 " available.")
            return ""

        if self._is_hbcmat:
            self._logman.warning("The source is a HilbertCurve matrix, no DNA"
                                 " sequence was available.")
            return ""

        return self._dnaseq

    @property
    def bits(self):
        return self._bits

    @property
    def hbcmat(self):
        return self.get_hbcmat()

    @property
    def dnaseq(self):
        return self.get_dnaseq()

    @property
    def kmered_seq(self):
        return self._kmer_biseq

    @property
    def allelic_attrs(self):
        return self._a1_val, self._a2_val


def create_database(vcf, bed, fasta, db_path, shift, meta=None, kmer: int = 4):
    """Make HDF5 database.

    Make a HDF5 database from given genetic variants (SNP) for specific genome
    regions according to given reference genome
    """
    def _parse_bed(itvl, sep="\t"):
        itvls = itvl.split(sep)
        if len(itvls) < 6:
            raise ValueError("Require the BED file with at least 6 fields.")
        chrom, start, end, name, _, strand, *_ = itvls
        return chrom, start, end, name, strand

    samples = None
    with VariantFile(vcf) as variants:
        samples = variants.header.samples

    with TabixFile(bed) as intervals, FastaFile(fasta) as sequence, \
            CSVDict(meta) as metadata:
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
                            hbc = HilbertCurve(sub_seqs, kmer, strand).hbcmat
                            hdf5db.add_matrix(hbc, per_sample, attr_lst)


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

    create_database(vcf=genetic_variants, bed=genome_intervals,
                    fasta=reference_genome, db_path=output_dir,
                    shift=n_base_pairs, meta=metadata_table)
