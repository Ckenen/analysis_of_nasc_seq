#!/usr/bin/env runsnakemake
import pandas as pd

runs = [
    "20220113_NASCseq_K562",
    "20220321_NASCseq_K562",
    "20220418_NASCseq_K562",
    "GSE128273_NASCseq_K562"]
dat = pd.read_excel("data/NASCseq.xlsx")
dat = dat[dat["Species"] == "Human"]
dat = dat[[r in runs for r in dat["Run"]]]
# print(len(dat))

run_cells = []
threads = 12
for run, cell in dat[["Run", "Cell"]].values:
    s = "%s/%s" % (run, cell)
    run_cells.append(s)

## 
BOWTIE2_RRNA_INDEX = "/home/chenzonggui/species/homo_sapiens/ncbi/bt2_rrna_index"
STAR_GENOME_INDEX = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.v39.STAR.index"
GENE_BED = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.genes.bed"
TRANSCRIPT_BED_GZ = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz"
ANNOTATION_TSV = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.tsv"
SNP_BED_GZ = "/home/chenzonggui/species/homo_sapiens/hg38/snp151.3.lite.bed.gz"
# config["bt2_rrna_idx"] = "/date/chenzonggui/species/homo_sapiens/ncbi/bt2_rrna_index"
# config["star_idx"] = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.v39.STAR.index"
# config["gene_bed"] = "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.genes.bed"
# config["transcript_bed"] = "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz"
# config["transcript_tsv"] = "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.tsv"
# config["snp_bed"] = "/date/chenzonggui/species/homo_sapiens/hg38/snp151.3.lite.bed.gz"