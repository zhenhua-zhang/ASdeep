import csv
import traceback
from collections import OrderedDict

from pysam.libcbcf import VariantFile
from pysam.libctabix import TabixFile, asBed, asGTF, asGFF3

class CSVDict(csv.DictReader):
    def __init__(self, csvfile, idx_col=0, row_key_tran=None, **kwargs):
        self._row_key_tran = row_key_tran

        if csvfile is None:
            self._csv_handler = None
            self._metadict = None
            self._colkeys = None
            self._rowkeys = None
        else:
            self._csv_handler = open(csvfile, "r")
            super(CSVDict, self).__init__(self._csv_handler, **kwargs)
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


class TabDict(TabixFile):
    def __init__(self, *args, **kwargs):
        self._rec_iters = self.fetch()
        self._args = args
        self._kwargs = kwargs

        if "parser" in kwargs:
            self._parser = kwargs["parser"]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)

        if hasattr(self, "close"):
            self.close()

    def _make_dict(self):
        raise NotImplementedError("The hiden method should be implemented "
                                  "accordingly.")

    @property
    def is_bed(self):
        return isinstance(self._parser, asBed)

    @property
    def is_gtf(self):
        return isinstance(self._parser, asGTF)

    @property
    def is_gff3(self):
        return isinstance(self._parser, asGFF3)

    @property
    def tabdict(self):
        return self._make_dict()

    def subset(self, **kwargs):
        if "parser" in kwargs:
            self._parser = kwargs["parser"]

        self._rec_iters = self.fetch(**kwargs)
        return self


class BEDDict(TabDict):
    def __init__(self, *args, **kwargs):
        super(BEDDict, self).__init__(*args, **kwargs)

    def _make_dict(self):
        bed_dict = OrderedDict()
        for per_rec in self._rec_iters:
            chrom, _, pos = per_rec.contig, per_rec.start, per_rec.end
            rsid, ref, alt, refrc, altrc, allrc, *oth = per_rec.name.split(";")
            key = (chrom, pos, ref, alt)
            bed_dict[key] = (chrom, pos, rsid, ref, alt, int(refrc),
                             int(altrc), int(allrc), oth)

        return bed_dict


class GTFDict(TabDict):
    def __init__(self, *args, **kwargs):
        self._gene_id = "gene_id"
        if "gene_id" in kwargs:
            self._gene_id = kwargs.pop("gene_id")

        self._transcript_id = "transcript_id"
        if "transcript_id" in kwargs:
            self._transcript_id = kwargs.pop("transcript_id")

        super(GTFDict, self).__init__(*args, **kwargs)

    @staticmethod
    def _parse_attrs(attr):
        fields = [x.strip() for x in attr.strip().split(";")]
        for per_field in fields:
            if not per_field:
                continue

            if per_field.endswith(";"):
                per_field = per_field[:-1]

            if "=" in per_field:
                key, val = [x.strip() for x in per_field.split("=", 1)]
            else:
                key, val = [x.strip() for x in per_field.split(" ", 1)]

            if val[0] == '"' and val[1] == '"':
                val = val[1:-1]
            else:
                try:
                    val = float(val)
                    val = int(val)
                except ValueError:
                    pass
                except TypeError:
                    pass

            yield key, val

    def _make_dict(self):
        gtf_dict = OrderedDict()
        for per_rec in self._rec_iters:
            if len(per_rec) != 9:
                raise ValueError("The GTF/GFF record should have 9 fields")
            
            seqname, source, feature, start, end, score, strand, frame, attr = per_rec

            attr_dict = dict(self._parse_attrs(attr))
            gene_id = attr_dict.get(self._gene_id, None)
            mrna_id = attr_dict.get(self._transcript_id, None)

            if mrna_id is None: continue

            rec_tuple = (seqname,
                         self._try_to_dot(source),
                         self._try_to_dot(feature),
                         int(start),
                         int(end),
                         self._try_to_dot(score),
                         self._try_to_dot(strand),
                         self._try_to_dot(frame),
                         attr_dict)

            if gene_id in gtf_dict:
                if mrna_id in gtf_dict[gene_id]:
                    gtf_dict[gene_id][mrna_id].append(rec_tuple)
                else:
                    gtf_dict[gene_id][mrna_id] = [rec_tuple]
            else:
                gtf_dict[gene_id] = {mrna_id: [rec_tuple]}

        return gtf_dict

    @staticmethod
    def _try_to_dot(x):
        if x is None: return "."
        else:
            try:
                return int(x)
            except ValueError:
                try:
                    return float(x)
                except ValueError:
                    return str(x)


class VCFDict:
    def __init__(self, *args, **kwargs):
        self._sample_id = kwargs.get("sample_id", None)
        if self._sample_id is not None:
            kwargs.pop("sample_id")

        self._vcf = VariantFile(*args, **kwargs)

        if self._sample_id:
            self._vcf.subset_samples([self._sample_id])

        self._rec_iters = self._vcf.fetch()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)

        if self._vcf.is_open:
            self._vcf.close()

    def _mk_vcfdict(self):
        vcf_dict = OrderedDict()
        sample_id = self._sample_id
        for per_rec in self._rec_iters:
            chrom, pos, rsid = per_rec.chrom, per_rec.pos, per_rec.id
            ref, alts = per_rec.ref, per_rec.alts

            if len(alts) != 1: continue
            alt = alts[0]

            sample = per_rec.samples.get(sample_id)
            phased = True if sample.phased else False
            genotype = sample.get("GT", None)

            key = (chrom, pos, ref, alt)
            if key in vcf_dict:
                raise KeyError(f"Duplicated key: {key}")

            vcf_dict[key] = (chrom, pos, rsid, ref, alt, phased, genotype)

        return vcf_dict

    @property
    def tabdict(self):
        return self._mk_vcfdict()

    def is_open(self):
        return self._vcf.is_open

    def close(self):
        if self._vcf.is_open:
            self._vcf.close()

    def subset(self, **kwargs):
        self._rec_iters = self._vcf.fetch(**kwargs)
        return self

