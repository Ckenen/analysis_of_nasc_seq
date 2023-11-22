#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/marked_strand"
outdir = "results/nascent"

rule all:
    input:
        #expand(outdir + "/events/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/ratio/{run_cell}.tsv", run_cell=run_cells),
        # expand(outdir + "/marked_nascent/{run_cell}.bam", run_cell=run_cells),
        # expand(outdir + "/expression/fpkm/{run_cell}.tsv", run_cell=run_cells[:1]),

# Events

rule get_events:
    input:
        bam = indir + "/{run}/{cell}.bam",
        bed = SNP_BED_GZ
    output:
        bam = outdir + "/events/{run}/{cell}.bam"
    log:
        outdir + "/events/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools GetEvent --threads {threads} --snp {input.bed} {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule report_mismatch:
    input:
        bam = rules.get_events.output.bam
    output:
        txt = outdir + "/ratio/{run}/{cell}.tsv"
    log:
        outdir + "/ratio/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools ReportMismatch --threads {threads} --strand TAG --strand-tag ST {input.bam} {output.txt} &> {log}
        """

# Nascent

rule mark_nascent:
    input:
        bam = rules.get_events.output.bam
    output:
        bam = outdir + "/marked_nascent/{run}/{cell}.bam"
    log:
        outdir + "/marked_nascent/{run}/{cell}.log"
    shell:
        """
        nasctools MarkNascent --platform NGS --layout PE {input.bam} {output.bam} &> {log}
        samtools index {output.bam}
        """

# FPKM

rule calculate_fpkm:
    input:
        bam = rules.mark_nascent.output.bam,
        bed = TRANSCRIPT_BED_GZ,
        txt = ANNOTATION_TSV
    output:
        txt = outdir + "/expression/fpkm/{run}/{cell}.tsv"
    log:
        outdir + "/expression/fpkm/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        nasctools CalculateFPKM --threads {threads} --layout PE --strand TAG --strand-tag ST --nascent \
            --annotation {input.txt} {input.bam} {input.bed} {output.txt} &> {log}
        """