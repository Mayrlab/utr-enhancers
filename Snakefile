configfile: "config.yaml"

import pandas as pd
import os

# set display options
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)


metadata = pd.read_csv(config["sample_sheet"], sep='\t',
                       index_col='sample_id').filter(regex='^(ery|lsk)', axis=0)

#print(metadata)
#exit(0)

rule all:
    input:
        expand("data/salmon/{lib_type}/{sample_id}", zip,
               sample_id=metadata.index.values,
               lib_type=metadata.library_type.values),
        expand("qc/gene_cov/{lib_type}/{sample_id}.geneBodyCoverage.curves.pdf", zip,
               sample_id=metadata.index.values,
               lib_type=metadata.library_type.values)
        #"qc/gene_cov/all.pdf"
        #expand("data/qapa/{sample_id}_pau_results.txt", sample_id=metadata.index.values)

rule download_fastq_se:        
    output:
        "data/fastq/se/{sample_id}.fastq.gz"
    params:
        tmp_dir=config['tmp_dir'],
        srr=lambda wcs: metadata.srr[wcs.sample_id]
    threads: 8
    resources:
        mem_mb=1000
    conda:
        "envs/samtools.yaml"
    shell:
        """
        TMP_FQ={params.tmp_dir}/{params.srr}.fastq
        fasterq-dump -e {threads} --split-files -t {params.tmp_dir} -O {params.tmp_dir} {params.srr}
        bgzip -c@ {threads} $TMP_FQ > {output}
        rm $TMP_FQ
        """

rule download_fastq_pe:        
    output:
        r1="data/fastq/pe/{sample_id}_1.fastq.gz",
        r2="data/fastq/pe/{sample_id}_2.fastq.gz"
    params:
        tmp_dir=config['tmp_dir'],
        srr=lambda wcs: metadata.srr[wcs.sample_id]
    threads: 8
    resources:
        mem_mb=1000
    conda:
        "envs/samtools.yaml"
    shell:
        """
        TMP_R1={params.tmp_dir}/{params.srr}_1.fastq
        TMP_R2={params.tmp_dir}/{params.srr}_2.fastq
        fasterq-dump -e {threads} --split-files -t {params.tmp_dir} -O {params.tmp_dir} {params.srr}
        bgzip -c@ {threads} $TMP_R1 > {output.r1}
        bgzip -c@ {threads} $TMP_R2 > {output.r2}
        rm $TMP_R1 $TMP_R2
        """

rule pe_to_se:
    input:
        r1="data/fastq/pe/{sample_id}_1.fastq.gz",
        script="scripts/truncate_reads.awk"
    output:
        fq="data/fastq/pe_se/{sample_id}.fastq.gz"
    params:
        len=51
    threads: 8
    resources:
        mem_mb=1000
    conda:
        "envs/samtools.yaml"
    shell:
        """
        zcat {input.r1} |\\
        awk -v len={params.len} -f {input.script} |\\
        bgzip -c@ {threads} > {output.fq}
        """

rule salmon_quant_se:
    input:
        fq="data/fastq/se/{sample_id}.fastq.gz",
        sdx=config['salmon_idx']
    output:
        directory("data/salmon/se/{sample_id}")
    params:
        extraFlags="--gcBias --validateMappings -l A"
    threads: 8
    resources:
        mem_mb=2000
    conda:
        "envs/salmon.yaml"
    shell:
        """
        salmon quant -p {threads} -i {input.sdx} {params.extraFlags} -r {input.fq} -o {output}
        """

rule salmon_quant_pe_se:
    input:
        fq="data/fastq/pe_se/{sample_id}.fastq.gz",
        sdx=config['salmon_idx']
    output:
        directory("data/salmon/pe_se/{sample_id}")
    params:
        extraFlags="--gcBias --validateMappings -l A"
    threads: 8
    resources:
        mem_mb=2000
    conda:
        "envs/salmon.yaml"
    shell:
        """
        salmon quant -p {threads} -i {input.sdx} {params.extraFlags} -r {input.fq} -o {output}
        """

rule salmon_quant_pe:
    input:
        r1="data/fastq/pe/{sample_id}_1.fastq.gz",
        r2="data/fastq/pe/{sample_id}_2.fastq.gz",
        sdx=config['salmon_idx']
    output:
        directory("data/salmon/pe/{sample_id}")
    params:
        extraFlags="--gcBias --validateMappings -l A"
    threads: 8
    resources:
        mem_mb=2000
    conda:
        "envs/salmon.yaml"
    shell:
        """
        salmon quant -p {threads} -i {input.sdx} {params.extraFlags} -1 {input.r1} -2 {input.r2} -o {output}
        """

rule ensembl_ids:
    input:
        config['ensembl_ids']
    output:
        "data/ensembl/ensembl_identifiers.txt"
    shell:
        """
        cp {input} {output}
        """
        
rule qapa_quant_utrs:
    input:
        sf=lambda wcs: "data/salmon/" + metadata.library_type[wcs.sample_id] + "/" + wcs.sample_id + "/quant.sf",
        ensembl=rules.ensembl_ids.output
    output:
        "data/qapa/{sample_id}_pau_results.txt"
    params:
        tmp_dir=config['tmp_dir']
#    conda:
#        "envs/qapa.yaml"
    shell:
        """
        source /home/fanslerm/software/miniconda3/etc/profile.d/conda.sh
        conda activate qapa
        qapa quant -t {params.tmp_dir} --db {input.ensembl} {input.sf} > {output}
        """

def strand_se (wcs):
    s = metadata.strand[wcs.sample_id]
    if s == 'f':
        return "--rna-strandness F"
    elif s == 'r':
        return "--rna-strandness R"
    else:
        return ""
    
rule hisat2_se:
    input:
        fq="data/fastq/se/{sample_id}.fastq.gz"
    output:
        bam="data/bam/se/{sample_id}.bam",
        bai="data/bam/se/{sample_id}.bam.bai"
    params:
        tmp=config['tmp_dir'],
        hisat2=config['hisat2'],
        idx=config['hisat2_idx'],
        sam=config['tmp_dir'] + "/{sample_id}.sam",
        strand=strand_se
    threads: 16
    resources:
        mem_mb=2000
    conda: "envs/samtools.yaml"
    shell:
        """
        {params.hisat2} -p {threads} -x {params.idx} {params.strand} -U {input} -S {params.sam}
        samtools sort -@ {threads} -m 1536M -T {params.tmp}/ -o {output.bam} {params.sam}
        samtools index -@ {threads} {output.bam}
        rm -f {params.sam}
        """

def strand_pe (wcs):
    s = metadata.strand[wcs.sample_id]
    if s == 'fr':
        return "--rna-strandness FR"
    elif s == 'rf':
        return "--rna-strandness RF"
    else:
        return ""

rule hisat2_pe:
    input:
        r1="data/fastq/pe/{sample_id}_1.fastq.gz",
        r2="data/fastq/pe/{sample_id}_2.fastq.gz"
    output:
        bam="data/bam/pe/{sample_id}.bam",
        bai="data/bam/pe/{sample_id}.bam.bai"
    params:
        tmp=config['tmp_dir'],
        hisat2=config['hisat2'],
        idx=config['hisat2_idx'],
        sam=config['tmp_dir'] + "/{sample_id}.sam",
        strand=strand_pe
    threads: 16
    resources:
        mem_mb=2000
    conda: "envs/samtools.yaml"
    shell:
        """
        {params.hisat2} -p {threads} -x {params.idx} {params.strand} -1 {input.r1} -2 {input.r2} -S {params.sam}
        samtools sort -@ {threads} -m 1536M -T {params.tmp}/ -o {output.bam} {params.sam}
        samtools index -@ {threads} {output.bam}
        rm -f {params.sam}
        """

rule gene_cov_all:
    input:
        bam=expand("data/bam/{lib_type}/{sample_id}.bam", zip,
               sample_id=metadata.index.values,
               lib_type=metadata.library_type.values),
        bed=config['gene_hk_bed']
    output:
        "qc/gene_cov/all.pdf"
    params:
        bams=lambda _, input: ",".join(map(str, input.bam))
    conda: "envs/rseqc.yaml"
    resources:
        mem_mb=8000,
        time_min=86400
    shell:
        """
        PDF={output}
        OUT=${{PDF%.pdf}}
        geneBody_coverage.py -r {input.bed} -i {params.bams} -o $OUT
        """

rule gene_cov_each:
    input:
        bam="data/bam/{lib_type}/{sample_id}.bam",
        bed=config['gene_hk_bed']
    output:
        "qc/gene_cov/{lib_type}/{sample_id}.geneBodyCoverage.curves.pdf"
    params:
        out="qc/gene_cov/{lib_type}/{sample_id}"
    conda: "envs/rseqc.yaml"
    resources:
        mem_mb=2000
    shell:
        """
        geneBody_coverage.py -r {input.bed} -i {input.bam} -o {params.out}
        """
