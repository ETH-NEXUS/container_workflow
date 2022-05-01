
rule fgbio_FastqToBam:
    input:
        unpack(get_fastq), 
    output:
        bam = "results/prepareBam/{sample}_unmapped.bam"
    log:
        "logs/fgbio_FastqToBam/{sample}.log",
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/fgbio:1.5.1--hdfd78af_0"
    params:
        variousParams=config["tools"]["fgbio"]["FastqToBam"]["variousParams"],
    threads: config['computingResources']['mediumRequirements']['threads']
    shell:
        "fgbio FastqToBam -i {input.fq1} {input.fq2} -o {output.bam} --sample {wildcards.sample} {params.variousParams} &> {log}"


rule picard_SamToFastq:
    input:
        bam = "results/prepareBam/{sample}_unmapped.bam",
    output:
        R1 = "results/prepareBam/{sample}_allocated_R1.fastq",
        R2 = "results/prepareBam/{sample}_allocated_R2.fastq"
    log:
        "logs/picard_SamToFastq/{sample}.log",
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/picard-slim:2.26.11--hdfd78af_0"
    params:
    threads: config['computingResources']['mediumRequirements']['threads']
    shell:
        "picard SamToFastq -I {input.bam} -F {output.R1} -F2 {output.R2} &> {log}"


# Create mapped bam

rule bwa_mem:
    input:
       R1 = "results/prepareBam/{sample}_allocated_R1.fastq",
       R2 = "results/prepareBam/{sample}_allocated_R2.fastq",
    output:
        bam = "results/prepareBam/{sample}_mapped.bam",
    log:
        "logs/bwa/{sample}.log",
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/bwa:0.7.17--h7132678_9"
    params:
        variousParams=config["tools"]["bwa"]["variousParams"],
        reference = config["referenceBwa"],
    threads: config['computingResources']['mediumRequirements']['threads']
    shell:
        "bwa mem {params.variousParams} -t {threads} {params.reference} {input.R1} {input.R2} > {output.bam} "

rule picard_SortSam_mapped:
    input:
        bam = "results/prepareBam/{sample}_mapped.bam", 
    output:
        bamSorted = "results/prepareBam/{sample}_mapped.sorted.bam",
    log:
        "logs/picard_SortSam_mapped/{sample}.log",
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/picard-slim:2.26.11--hdfd78af_0"
    params:
    threads: config['computingResources']['mediumRequirements']['threads']
    shell:
        "picard SortSam -I {input.bam} -O {output.bamSorted} -SORT_ORDER 'queryname' &> {log}"

rule picard_SortSam_unmapped:
    input:
        bam = "results/prepareBam/{sample}_unmapped.bam"
    output:
        bamSorted = "results/prepareBam/{sample}_unmapped.sorted.bam",
    log:
        "logs/picard_SortSam_unmapped/{sample}.log",
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/picard-slim:2.26.11--hdfd78af_0"
    params:
    threads: config['computingResources']['mediumRequirements']['threads']
    shell:
        "picard SortSam -I {input.bam} -O {output.bamSorted} -SORT_ORDER 'queryname' &> {log}"

# merge Bam
rule picard_MergeBamAlignment:
    input:
        unmappedBam = "results/prepareBam/{sample}_unmapped.sorted.bam",
        mappedBam = "results/prepareBam/{sample}_mapped.sorted.bam",
    output:
        bam = "results/prepareBam/{sample}_prepared.bam",
    log:
        "logs/picard_MergeBamAlignment/{sample}.log",
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/picard-slim:2.26.11--hdfd78af_0"
    params:
        reference = config["reference"],
    threads: config['computingResources']['mediumRequirements']['threads']
    shell:
        "picard MergeBamAlignment -UNMAPPED_BAM {input.unmappedBam} -ALIGNED_BAM {input.mappedBam} -REFERENCE_SEQUENCE {params.reference} -O {output.bam} &> {log}"

# Mark duplicates
rule picard_MarkDuplicates:
    input:
        bam = "results/prepareBam/{sample}_prepared.bam",
    output:
        bam = "results/prepareBam/{sample}_MarkedDuplicates.bam",
        metrics = "results/prepareBam/{sample}_MarkedDuplicates.Metrics.txt"
    log:
        "logs/picard_MarkDuplicates/{sample}.log",
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/picard-slim:2.26.11--hdfd78af_0"
    params:
    threads: config['computingResources']['mediumRequirements']['threads']
    shell:
        "picard MarkDuplicates BARCODE_TAG=RX I={input.bam} O={output.bam} M={output.metrics} &> {log}"

rule fgbio_SortBam:
    input:
       bam = "results/prepareBam/{sample}_MarkedDuplicates.bam",
    output:
       bam = "results/callConsensus/{sample}_MarkedDuplicates.sorted.bam"
    log:
        "logs/fgbio_SortBam/{sample}.log",
    singularity:
        "https://nexus-central.leomed.ethz.ch/repository/galaxyproject-raw/fgbio:1.5.1--hdfd78af_0"
    params:
    threads: config['computingResources']['mediumRequirements']['threads']
    shell:
        "fgbio SortBam -s TemplateCoordinate -i {input.bam} -o {output.bam} &> {log}"
