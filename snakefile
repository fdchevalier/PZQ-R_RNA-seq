import os
import fnmatch

SAMPLES = os.listdir("data/libraries")

rule all:
    input:
        expand("data/libraries/{sample}/{sample}_Aligned.toTranscriptome.out.bam", sample=SAMPLES),
        expand("data/libraries/{sample}/{sample}_rsem.genes.results", sample=SAMPLES)

rule alignment:
    input:
        read1="data/libraries/{sample}/{sample}_R1.fastq.gz",
        read2="data/libraries/{sample}/{sample}_R2.fastq.gz",
        genome="data/genome/S.mansoni_STAR_149"
    output:
        "data/libraries/{sample}/{sample}_Aligned.toTranscriptome.out.bam"
    params:
        output=r"data/libraries/{sample}/{sample}_"
    shell:
        'STAR --runMode alignReads \
              --runThreadN {resources.cores} \
              --genomeDir "{input.genome}"   \
              --readFilesIn "{input.read1}" "{input.read2}" \
              --readFilesCommand zcat \
              --outFileNamePrefix "{params.output}" \
              --outSAMtype BAM SortedByCoordinate   \
              --quantMode TranscriptomeSAM GeneCounts'

rule count:
    input:
        bam="data/libraries/{sample}/{sample}_Aligned.toTranscriptome.out.bam"
    output:
        "data/libraries/{sample}/{sample}_rsem.genes.results"
    params:
        genome="data/genome/S.mansoni_RSEM/S.mansoni",
        output=r"data/libraries/{sample}/{sample}_rsem"
    shell:
        'rsem-calculate-expression --alignments --paired-end -p {resources.cores} --no-bam-output \
                "{input.bam}"    \
                "{params.genome}" \
                "{params.output}"'

