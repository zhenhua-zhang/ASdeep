"""ASE report.
"""

import os
import math
import logging

import pandas
import seaborn as sbn
import gffutils as gut
import matplotlib.pyplot as plt


class ASEReport:
    """A class to process ASE report.
    """

    def __init__(self, file_path_list, save_prefix, **kwargs):
        self.file_path_list = file_path_list
        self.save_prefix = save_prefix

        self.pvm_raw = None
        self.pvm_by_coord = None

        self.genome_annot = kwargs.get("genome_annot", None)
        self.ant_db_name = None
        self.ant_sql = None

    def _pr_join_dtfm(self, **kwargs):
        """Merge multiple data frame of ASE quantification report.

        Parameters
        ------
        kwargs (optional): Parameters(i.e. sep, usecols) used in pandas.read_csv()
        """
        sep = kwargs.get("sep", "\t")
        use_cols = kwargs.get("usecols", None)
        dtfm_list = [pandas.read_csv(file_path, sep=sep, usecols=use_cols)
                     for file_path in self.file_path_list]

        if dtfm_list:
            return pandas.concat(dtfm_list, ignore_index=True)

        logging.warning("The dataframe list is empty.")
        return None

    @staticmethod
    def _pr_check_read_count(ase_rec: pandas.Series, gross_count=10, allelic_count=3):
        allele_count_str = ase_rec["allele_counts"]
        allele_count_lst = allele_count_str.split(":")[1].split(";")

        a1_sum, a2_sum = 0, 0
        for count_pair in allele_count_lst:
            a1_count, a2_count = count_pair.split("|")
            a1_sum += int(a1_count)
            a2_sum += int(a2_count)

        return (a1_sum >= allelic_count and a2_sum >= allelic_count) or ((a1_sum + a2_sum) > gross_count)

    def _pr_fetch_gene_coord(self, gene_id):
        gene = self.ant_sql[gene_id]    # Feature attributes, accessed by gene.attributes["gene_name"]
        return gene_id, gene.chrom, gene.start, gene.end

    def _pr_p_val_matrix(self, max_na_per_gene=1500):
        # Get P values.
        merged_dtfm = self._pr_join_dtfm()
        if merged_dtfm is not None:
            pvm_raw = (merged_dtfm
                       .loc[:, ["gene_id", "sample_id", "p_val"]]
                       .pivot(index="gene_id", columns="sample_id"))

            # Remove genes without herterozygous locus in `max_na_per_gene` individuals at maximum
            pvm_by_coord = (pvm_raw.loc[pvm_raw.isna().sum(axis=1) < max_na_per_gene, :].fillna(1))

            gene_coord = [self._pr_fetch_gene_coord(gene_id) for gene_id in pvm_by_coord.index]
            gene_coord = sorted(gene_coord, key=lambda x: (int(x[1]), int(x[2])))

            pvm_by_coord = pvm_by_coord.loc[[coord[0] for coord in gene_coord], :]
            pvm_by_coord["chrom"] = [coord[1] for coord in gene_coord]
            pvm_by_coord["pos"] = [coord[2] for coord in gene_coord]

            return pvm_raw, pvm_by_coord

        return None, None

    def _pr_draw_p_val_htmp(self):
        # Transform p-values by -log10()
        p_val_matrix = (self.pvm_raw.apply(lambda x: [-math.log10(e) if e>5e-6 else 6 for e in x]))

        # Heatmap
        fig_size = [x / 2 if x < 200 else 100 for x in p_val_matrix.shape]
        fig_size[0], fig_size[1] = fig_size[1], fig_size[0]
        try:
            grid = sbn.clustermap(p_val_matrix, figsize=fig_size, cmap="Greens")
        except:
            try:
                logging.warning("Failed to cluster rows or columns, try again without cluster.")
                grid = sbn.clustermap(p_val_matrix, figsize=fig_size, cmap="Greens",
                                      row_cluster=False, col_cluster=False)
            except:
                logging.error("Failed to do cluster, skip heatmap.")
                grid = None

        return grid

    # Show the number of individuals is ASE for the gene.
    def _pr_draw_ase_gene_count_across_genome(self):
        mtrx = self.pvm_by_coord

        chr_len_pl = mtrx["chrom"].value_counts()
        chr_len_pl.index = [int(x) for x in chr_len_pl.index]
        chr_len_pl = chr_len_pl.sort_index()

        chrom_vspan = [
            [0, xmax] if idx == 0 else [sum(chr_len_pl[:idx]), sum(chr_len_pl[:idx+1])]
            for idx, xmax in enumerate(chr_len_pl)
        ]

        mtrx = mtrx.loc[:, [x for x in mtrx.columns if x not in [("chrom", ""), ("pos", "")]]]
        ase_gene_count = mtrx[mtrx < 5e-6].count(axis=1) / (mtrx.shape[1] - 1) * 100

        fig, axes = plt.subplots()
        axes.set_title('Percentage of ASE genes across the genome')
        axes.plot(ase_gene_count.index, ase_gene_count, linewidth=0.5)

        for idx, (xmin, xmax) in enumerate(chrom_vspan):
            color = "g" if idx % 2 == 0 else "r"
            axes.axvspan(xmin=xmin, xmax=xmax, ymin=0, facecolor=color, alpha=0.1)

        chr_xticks = [
            xmax / 2 if idx == 0 else (sum(chr_len_pl[:idx]) + sum(chr_len_pl[:idx+1])) / 2
            for idx, xmax in enumerate(chr_len_pl)
        ]
        axes.set_xticks(chr_xticks)

        chr_labels = ["chr" + str(x) for x in chr_len_pl.index]
        axes.set_xticklabels(chr_labels, rotation=45, rotation_mode='anchor', ha='right')

        axes.set_xlabel("Genome coordination (gene binned by position)")
        axes.set_ylabel("Percentage (out of " + str(len(ase_gene_count.index)) + ')')

        return fig

    # Generate descriptive statistic for the merged ASE report.
    def _pr_desc_stat(self, alpha=5e-6):
        if self.pvm_raw is not None:

            ase_gene_per_individual = self.pvm_raw[self.pvm_raw < alpha].count(axis=0)
            individual_per_ase_gene = self.pvm_raw[self.pvm_raw < alpha].count(axis=1)

            return ase_gene_per_individual, individual_per_ase_gene

        return None, None

    def init(self, **kwargs):
        """Initialize the instance.
        """
        logging.info("Init")
        # Create a GFF/GTF database to fetch coordination of given gene.
        ant_db_name = kwargs.get("ant_db_name", None)
        max_na_per_gene = kwargs.get("max_na_per_gene", 1500)

        if ant_db_name is None:
            ant_db_name = os.path.splitext(self.genome_annot)[0] + ".db"

        if not os.path.exists(ant_db_name):
            gut.create_db(self.genome_annot, ant_db_name, disable_infer_transcripts=True,
                          disable_infer_genes=True)

        self.ant_sql = gut.FeatureDB(ant_db_name)
        self.pvm_raw, self.pvm_by_coord = self._pr_p_val_matrix(max_na_per_gene)

        return self

    def report(self, report_fmt="csv"):
        """Generate files to show quantification results.
        """
        logging.info("Report")
        save_prefix = self.save_prefix
        if self.pvm_raw is not None:
            if report_fmt == "csv":
                sep = ","
            elif report_fmt == "tsv":
                sep = "\t"
            else:
                logging.warning("Unknown format: {}, using csv as default".format(report_fmt))
                sep = ","

            self.pvm_raw.to_csv(save_prefix + "-pval_matrix_raw." + report_fmt, sep=sep)
            self.pvm_by_coord.to_csv(save_prefix + "-pval_matrix_sorted_filtered." + report_fmt, sep=sep)

        return self

    def visualize(self, fig_fmt: str = "pdf"):
        """Draw figures to show the result.
        """
        logging.info("Visualize")
        save_prefix = self.save_prefix

        if self.pvm_raw is not None:
            # grid = self._pr_draw_p_val_htmp()
            # if grid:
            #     grid.fig.savefig(save_prefix + "-heatmap." + fig_fmt)
            #     plt.close("all")

            fig = self._pr_draw_ase_gene_count_across_genome()
            fig.set_figheight(9)
            fig.set_figwidth(16)
            fig.savefig(save_prefix + "-ase_gene_count_across_genome." + fig_fmt)
            plt.close("all")

            ase_gene_per_individual, individual_per_ase_gene = self._pr_desc_stat()
            axes = ase_gene_per_individual.plot(kind="hist")
            axes.set_title("ASE genes per individual (p<5e-6)")
            fig = axes.get_figure()
            fig.savefig(save_prefix + "-ase_gene_per_individual." + fig_fmt)
            plt.close("all")

            axes = individual_per_ase_gene.plot(kind="hist")
            axes.set_title("Individuals per ASE gene (p<5e-6)")
            fig = axes.get_figure()
            fig.savefig(save_prefix + "-individual_per_ase_gene." + fig_fmt)
            plt.close("all")

        return self
