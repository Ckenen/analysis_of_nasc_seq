#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/prepare"

rule all:
    input:
        # expand(outdir + "/cutadapt/{run_cell}_1.fastq.gz", run_cell=run_cells),
        expand(outdir + "/bowtie2/{run_cell}.bam", run_cell=run_cells),


rule cutadapt:
    input:
        fq1 = "data/datasets/{cell}_1.fastq.gz",
        fq2 = "data/datasets/{cell}_2.fastq.gz"
    output:
        fq1 = outdir + "/cutadapt/{run}/{cell}_1.fastq.gz",
        fq2 = outdir + "/cutadapt/{run}/{cell}_2.fastq.gz"
    log:
        outdir + "/cutadapt/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 \
            -g GTGTATAAGAGACAG -g ATCAACGCAGAGTAC -a CTGTCTCTTATACAC -a GTACTCTGCGTTGAT \
            -G GTGTATAAGAGACAG -G ATCAACGCAGAGTAC -A CTGTCTCTTATACAC -A GTACTCTGCGTTGAT \
            -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}
        """

rule cutadapt_GSE128273:
    input:
        fq1 = "data/GSE128273_NASCseq_K562/{cell}_1.fastq.gz",
        fq2 = "data/GSE128273_NASCseq_K562/{cell}_2.fastq.gz"
    output:
        fq1 = outdir + "/cutadapt/GSE128273_NASCseq_K562/{cell}_1.fastq.gz",
        fq2 = outdir + "/cutadapt/GSE128273_NASCseq_K562/{cell}_2.fastq.gz"
    log:
        outdir + "/cutadapt/GSE128273_NASCseq_K562/{cell}.log"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -q 30 -m 20 \
            -g GTGTATAAGAGACAG -g GCAGAGTACGGG -a CTGTCTCTTATACAC -a CCCGTACTCTGC \
            -G GTGTATAAGAGACAG -G GCAGAGTACGGG -A CTGTCTCTTATACAC -A CCCGTACTCTGC \
            -o {output.fq1} -p {output.fq2} {input.fq1} {input.fq2} &> {log}
        """

rule bowtie2: # remove rRNA
    input:
        fq1 = outdir + "/cutadapt/{run}/{cell}_1.fastq.gz",
        fq2 = outdir + "/cutadapt/{run}/{cell}_2.fastq.gz",
        idx = BOWTIE2_RRNA_INDEX
    output:
        bam = outdir + "/bowtie2/{run}/{cell}.bam"
    log:
        outdir + "/bowtie2/{run}/{cell}.log"
    params:
        prefix = outdir + "/bowtie2/{run}/{cell}"
    threads:
        8
    shell:
        """(
        bowtie2 -p {threads} --local --no-unal --un-conc-gz {params.prefix}.fastq.gz -x {input.idx}/ref -1 {input.fq1} -2 {input.fq2} \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {params.prefix}_TMP - > {output.bam}
        samtools index -@ {threads} {output.bam}
        mv {params.prefix}.fastq.1.gz {params.prefix}.1.fastq.gz
        mv {params.prefix}.fastq.2.gz {params.prefix}.2.fastq.gz ) &> {log}
        """