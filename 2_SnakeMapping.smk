#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/prepare/bowtie2"
outdir = "results/mapping"

rule all:
    input:
        expand(outdir + "/star/{run_cell}", run_cell=run_cells),
        # expand(outdir + "/filtered/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/filtered/{run_cell}.flagstat", run_cell=run_cells),
        expand(outdir + "/infer_experiment/{run_cell}.txt", run_cell=run_cells),
        # expand(outdir + "/marked_duplicates/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/marked_duplicates/{run_cell}.flagstat", run_cell=run_cells),
        expand(outdir + "/marked_strand/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/marked_strand/{run_cell}.flagstat", run_cell=run_cells),
        outdir + "/all_samples_final_reads.tsv",
 
# STAR --genomeLoad LoadAndExit --genomeDir /home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.v39.STAR.index --limitBAMsortRAM 150000000000
# STAR --genomeLoad Remove --genomeDir /home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.v39.STAR.index

rule star_mapping:
    input:
        fq1 = indir + "/{run}/{cell}.1.fastq.gz",
        fq2 = indir + "/{run}/{cell}.2.fastq.gz",
        idx = STAR_GENOME_INDEX
    output:
        out = directory(outdir + "/star/{run}/{cell}")
    log:
        outdir + "/star/{run}/{cell}.log"
    params:
        prefix = outdir + "/star/{run}/{cell}/{cell}"
    threads:
        12
    shell:
        """(
        mkdir -p {output.out}
        STAR --runThreadN {threads} \
            --outFileNamePrefix {params.prefix}. \
            --genomeDir {input.idx} \
            --readFilesCommand zcat \
            --outSAMattributes All \
            --outSAMtype BAM SortedByCoordinate \
            --genomeLoad LoadAndKeep \
            --limitBAMsortRAM 150000000000 \
            --readFilesIn {input.fq1} {input.fq2}
        samtools index -@ {threads} {params.prefix}.Aligned.sortedByCoord.out.bam ) &> {log}
        """

rule filter_bam:
    input:
        bamdir = rules.star_mapping.output.out
    output:
        bam = outdir + "/filtered/{run}/{cell}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} --expr 'rname =~ "^chr([0-9]+|[XY])$"' -q 30 -d "NH:1" -f 2 -F 2308 \
            -o {output.bam} {input.bamdir}/{wildcards.cell}.Aligned.sortedByCoord.out.bam
        samtools index -@ {threads} {output.bam}
        """

rule infer_experiment:
    input:
        bam = rules.filter_bam.output.bam,
        bed = GENE_BED
    output:
        txt = outdir + "/infer_experiment/{run}/{cell}.txt"
    shell:
        """
        infer_experiment.py -s 2000000 -i {input.bam} -r {input.bed} > {output.txt} 2> /dev/null
        """

rule mark_duplicates:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/marked_duplicates/{run}/{cell}.bam",
        txt = outdir + "/marked_duplicates/{run}/{cell}_metrics.txt"
    log:
        outdir + "/marked_duplicates/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true -I {input.bam} -M {output.txt} -O {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_strand:
    input:
        bam = rules.mark_duplicates.output.bam,
        bed = TRANSCRIPT_BED_GZ
    output:
        bam = outdir + "/marked_strand/{run}/{cell}.bam",
        txt = outdir + "/marked_strand/{run}/{cell}.tsv"
    log:
        outdir + "/marked_strand/{run}/{cell}.log"
    shell:
        """
        nasctools MarkStrand --gene {input.bed} --tag ST --layout PE --strand U \
            --summary {output.txt} {input.bam} {output.bam} &> {log}
        samtools index {output.bam}
        """

rule report_reads:
    input:
        expand(outdir + "/marked_strand/{run_cell}.tsv", run_cell=run_cells),
    output:
        tsv = outdir + "/all_samples_final_reads.tsv"
    shell:
        """
        echo -e 'Name\\tForward\\tReverse\\tAmbiguous\\tUnknown' > {output.tsv}
        cat results/mapping/marked_strand/*/*.tsv | grep -v 'Name' >> {output.tsv}
        """

# common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """
