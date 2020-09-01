configfile: "config.yaml"

import pandas as pd
import os

# set display options
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)


metadata = pd.read_csv(config["sample_sheet"], sep='\t',
                       index_col='sample_id')
#print(metadata)
#exit(0)

rule all:
    input:
        expand("data/salmon/{lib_type}/{sample_id}", zip,
               sample_id=metadata.index.values,
               lib_type=metadata.library_type.values)#,
        #expand("data/qapa/{sample_id}_pau_results.txt", sample_id=metadata.index.values)

rule download_fastq_se:        
    output:
        "data/fastq/se/{sample_id}.fastq.gz"
    params:
        tmp_dir=config['tmp_dir'],
        srr=lambda wcs: metadata.srr[wcs.sample_id]
    threads: 8
    resources:
        mem=1
    shell:
        """
        FQ_DIR=$(dirname {output})
        fasterq-dump --split-files -t {params.tmp_dir} -O $FQ_DIR {params.srr}
        pigz -p {threads} $FQ_DIR/{params.srr}.fastq
        mv $FQ_DIR/{params.srr}.fastq.gz {output}
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
        mem=1
    shell:
        """
        FQ_DIR=$(dirname {output.r1})
        fasterq-dump --split-files -t {params.tmp_dir} -O $FQ_DIR {params.srr}
        pigz -p {threads} $FQ_DIR/{params.srr}_{{1,2}}.fastq
        mv $FQ_DIR/{params.srr}_1.fastq.gz {output.r1}
        mv $FQ_DIR/{params.srr}_2.fastq.gz {output.r2}
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
        mem=2
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
        mem=2
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
    conda:
        "envs/qapa.yaml"
    shell:
        """
        source /home/fanslerm/software/miniconda3/etc/profile.d/conda.sh
        conda activate qapa
        qapa quant -t {params.tmp_dir} --db {input.ensembl} {input.sf} > {output}
        """
