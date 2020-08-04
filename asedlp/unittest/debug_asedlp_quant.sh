ml Python/3.6.3-foss-2015b
source ~/Documents/projects/wp_ase_dlp/scripts/.env/bin/activate

chrom_id=1
fastq_id=AC43LCACXX-8-23
gff_file=~/Documents/projects/wp_ase_dlp/inputs/Ensembl_references/Homo_sapiens.GRCh37.75.gtf
sample_id_file=~/Documents/projects/wp_ase_dlp/misc/freeze2_GoNL_related_GTE_20200211_QCpassed.csv
work_dir=~/Documents/projects/wp_ase_dlp/workdir

sample_id=$(grep -w "$fastq_id" "$sample_id_file" | cut -f2)
haplotypes=$work_dir/snph5db/$chrom_id/haplotype.h5
snp_tab=$work_dir/snph5db/$chrom_id/snps_tab.h5
snp_index=$work_dir/snph5db/$chrom_id/snps_index.h5
sequence_tab=$work_dir/fastah5db/chr$chrom_id.h5
ref_read_counts=$work_dir/optdir/$fastq_id/waspOptdir/perChrom/$chrom_id/${fastq_id}_$chrom_id.refAlleleCounts.h5
alt_read_counts=$work_dir/optdir/$fastq_id/waspOptdir/perChrom/$chrom_id/${fastq_id}_$chrom_id.altAlleleCounts.h5
ase_report=ase_report.txt
train_set=matrix_and_ase.npz

python -m pdb ~/Documents/projects/wp_ase_dlp/scripts/asedlp/asedlp quant \
        --gene-ids ENSG00000188976 ENSG00000162591 ENSG00000116213 ENSG00000162408 \
        --sample-id $sample_id \
        --haplotypes $haplotypes \
        --snp-tab $snp_tab \
        --snp-index $snp_index \
        --sequence-tab $sequence_tab \
        --ref-read-counts $ref_read_counts \
        --alt-read-counts $alt_read_counts \
        --genome-annot $gff_file
