"""Infer allelic imbalance effects from allelic read counts using Bayesian inference."""

import os
import shutil
import logging
import tempfile
import traceback
from argparse import Namespace
from collections import OrderedDict

import arviz as az
import pymc3 as pm
from pysam.libctabix import asBed, asTuple

from .zutils import LogManager
from .tabdict import BEDDict
from .tabdict import GTFDict
from .tabdict import VCFDict

# Suppress the loging of pymc3
logger = logging.getLogger("pymc3")
logger.propagate = False
logger.setLevel(logging.ERROR)


class AllelicCounts:
    """A class to handling allelic read counts."""
    def __init__(self, sample_id: str, vcf_path: str, gtf_path: str,
                 bed_path: str, threads: int = 4, tar_feature: str = "exon",
                 hdi_prob: float = 0.90,
                 logman: LogManager = LogManager("RCPool")):
        self._logman = logman
        self._sample_id = sample_id
        self._tar_feature = tar_feature
        self._readcounts: OrderedDict = OrderedDict()
        self._ai_summary: list = []
        self._mrna_id: list = []

        self._vcf_recs = VCFDict(vcf_path, mode="r", sample_id=sample_id,
                                 threads=threads)
        self._gtf_recs = GTFDict(gtf_path, mode="r", parser=asTuple(),
                                 threads=threads)
        self._bed_recs = BEDDict(bed_path, mode="r", parser=asBed(),
                                 threads=threads)

        self._hdi_prob = hdi_prob

        self._trace = None

    def __bool__(self):
        return len(self._readcounts) > 0

    def __getitem__(self, idx):
        if isinstance(idx, str):
            return self._readcounts[idx]
        elif isinstance(idx, tuple) and len(idx) == 2:
            gene_id, mrna_id = idx
            return self._readcounts[gene_id][mrna_id]

        raise KeyError("Unsupported way to index.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)

        if self._vcf_recs.is_open():
            self._vcf_recs.close()

        if self._gtf_recs.is_open():
            self._gtf_recs.close()

        if self._bed_recs.is_open():
            self._bed_recs.close()

    @property
    def model(self):
        return self._model

    @property
    def trace(self):
        return self._trace

    @property
    def results(self):
        return self._ai_summary

    def fetch(self, **kwargs):
        gtf_recs = self._gtf_recs.subset(**kwargs)
        for gene_id, per_gene_rec in gtf_recs.tabdict.items():
            if gene_id in self._readcounts:
                raise KeyError(f"Duplicated entry: {gene_id}")
            self._readcounts[gene_id] = OrderedDict()

            for mrna_id, per_mrna_rec in per_gene_rec.items():
                if mrna_id in self._readcounts[gene_id]:
                    raise KeyError(f"Duplicated entry: {mrna_id}")

                self._readcounts[gene_id][mrna_id] = OrderedDict()
                a12 = []
                for per_exon_rec in per_mrna_rec:
                    chrom, _, feature, start, end, *_ = per_exon_rec

                    if feature != self._tar_feature:
                        continue

                    region = f"{chrom}:{start}-{end}"
                    rcpool = self._bed_recs.subset(region=region).tabdict
                    variants = self._vcf_recs.subset(region=region).tabdict

                    for key, per_var in variants.items():
                        if key not in rcpool:
                            continue

                        per_rc = rcpool[key]
                        chrom, pos, rsid, ref, alt, refrc, altrc, *_ = per_rc
                        _, _, _, _, _, is_phase, (a1_idx, a2_idx) = per_var

                        if a1_idx == a2_idx:
                            continue

                        if is_phase:
                            birc = (refrc, altrc)
                            a1_rc, a2_rc = birc[a1_idx], birc[a2_idx]
                            phase = f"{a1_idx}|{a2_idx}"
                        else:
                            a1_rc, a2_rc = altrc, refrc
                            phase = f"{a1_idx}/{a2_idx}"

                        a12.append((chrom, pos, ref, alt, rsid, phase, a1_rc,
                                    a2_rc))

                self._mrna_id.append((gene_id, mrna_id))
                self._readcounts.update({gene_id: {mrna_id: a12}})

        return self

    def inferai(self, gene_id: list = None, mrna_id: list = None,
                ab_sigma: float = 10, hdi_prob: float = None, **kwargs):
        """Infer allelic difference using MCMC."""

        if "tune" not in kwargs:
            kwargs["tune"] = 500

        if "draws" not in kwargs:
            kwargs["draws"] = 500

        if "chains" not in kwargs:
            kwargs["chains"] = 2

        if "cores" not in kwargs:
            kwargs["cores"] = kwargs["chains"]

        if "progressbar" not in kwargs:
            kwargs["progressbar"] = False

        if "random_seed" not in kwargs:
            kwargs["random_seed"] = 42

        if "return_inferencedata" not in kwargs:
            kwargs["return_inferencedata"] = True

        if isinstance(hdi_prob, float) and 0 < hdi_prob < 1:
            self._hdi_prob = hdi_prob

        if gene_id is None:
            gene_id = [x[0] for x in self._mrna_id]

        if mrna_id is None:
            mrna_id = [x[1] for x in self._mrna_id]

        for per_gene_id, per_mrna_id in self._mrna_id:
            if per_gene_id not in gene_id or per_mrna_id not in mrna_id:
                continue

            a12 = self._readcounts[per_gene_id][per_mrna_id]
            a1_rc, a2_rc = [x[-2] for x in a12], [x[-1] for x in a12]

            if not a1_rc or not a2_rc:
                self._logman.warning(f"No enough reads for {per_mrna_id} of "
                                     f"{per_gene_id}")
                self._ai_summary.append((per_gene_id, per_mrna_id, "", "",
                                     "", "", ""))
                continue

            n = sum(a1_rc) + sum(a2_rc)
            k = sum(a1_rc)

            self._model = pm.Model()
            with self._model:
                alpha = pm.HalfNormal(name="alpha", sigma=ab_sigma)
                beta = pm.HalfNormal(name="beta", sigma=ab_sigma)
                theta = pm.Beta(name="theta", alpha=alpha, beta=beta)
                _ = pm.Binomial(name="exp", p=theta, n=n, observed=k)

                self._trace = pm.sample(**kwargs)

            obs_data = ";".join([":".join([str(i) for i in x]) for x in a12])

            ai_summary = az.summary(self._trace, hdi_prob=self._hdi_prob)
            mean, sd, hdi_lower, hdi_upper = ai_summary.iloc[-1, :4]
            self._ai_summary.append((per_gene_id, per_mrna_id, mean, sd,
                                     hdi_lower, hdi_upper, obs_data))

        return self

    def save_to_dist(self, fpath: str, sep=","):
        lower_bound = (1 - self._hdi_prob) / 2
        upper_bound = 1 - lower_bound
        header = sep.join(["gene_id", "mrna_id", "mean", "sd",
                           f"hdi_{lower_bound:.3}",
                           f"hdi_{upper_bound:.3}", "evidence"])

        with open(fpath, "w") as fhandle:
            fhandle.write(header + "\n")
            for rec in self._ai_summary:
                line = sep.join([str(x) for x in rec]) + "\n"
                fhandle.write(line)

        return self


def inferai(args: Namespace, logman: LogManager = LogManager("InferAI")):
    """Infer the allelic imbalance using Bayesian inference."""
    n_cpu = args.n_cpu
    n_draw = args.n_draw
    n_tune = args.n_tune
    n_chain = args.n_chain
    hdi_prob = args.hdi_prob
    bed_path = args.readcounts_table
    gtf_path = args.genome_intervals
    vcf_path = args.genetic_variants
    sample_id = args.sample_id
    tar_feature = args.feature
    out_file = args.out_file

    logman.info(f"HDI              : {hdi_prob}")
    logman.info(f"N cpu            : {n_cpu}")
    logman.info(f"N draws          : {n_draw}")
    logman.info(f"N tunes          : {n_tune}")
    logman.info(f"N chains         : {n_chain}")
    logman.info(f"Variants         : {vcf_path}")
    logman.info(f"Sample ID        : {sample_id}")
    logman.info(f"Output file      : {out_file}")
    logman.info(f"Genome region    : {gtf_path}")
    logman.info(f"Target features  : {tar_feature}")
    logman.info(f"Rreadcounts table: {bed_path}")

    out_dir, _ = os.path.split(out_file)
    out_dir = os.path.realpath(out_dir)
    cache_path = tempfile.mkdtemp(prefix="theano-", dir=out_dir)
    os.environ["THEANO_FLAGS"] = f"base_compiledir={cache_path}"

    with AllelicCounts(sample_id=sample_id,
                       vcf_path=vcf_path,
                       gtf_path=gtf_path,
                       bed_path=bed_path,
                       tar_feature=tar_feature,
                       hdi_prob=hdi_prob) as allelic_counts:
        (allelic_counts
         .fetch()
         .inferai(draws=n_draw, chains=n_chain, tune=n_tune, cores=n_cpu)
         .save_to_dist(out_file))

    # Clean-up the cache path
    shutil.rmtree(cache_path, ignore_errors=True)
