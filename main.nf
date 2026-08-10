#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
================================================================================
PRS WGS restoration: SPRING -> germlineNGS-style CRAM

Purpose
-------
Restore ONLY the missing WGS CRAM files from archived SPRING files.

The alignment/preprocessing chain mirrors the relevant germlineNGS steps:

SPRING
  -> FASTQ.gz
  -> GATK FastqToSam
  -> GATK MarkIlluminaAdapters
  -> GATK SamToFastq
  -> BWA MEM 0.7.17 container
  -> GATK MergeBamAlignment
  -> samblaster
  -> sambamba view/sort
  -> CRAM
  -> CRAI

Only the final CRAM + CRAI are published.
FASTQ/BAM/intermediate files remain temporary in the Nextflow work directory
and are removed by the process after successful CRAM creation.

The nine SPRING files are explicitly listed below so no other samples are run.
================================================================================
*/

params.outdir      = "/lnx01_data2/shared/projects/raspau/PRS_WGS_restore_nextflow/cram"
params.spring_exec = "/lnx01_data2/shared/testdata/SPRING/build/spring"

/*
 * These variables (s_bind, simgpath, gatk_image, genome_fasta,
 * genome_version and tmpDIR) are supplied by the same shared config used by
 * germlineNGS: /data/configs/nextflowNGS.config
 */

def gatk_exec = "singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk"


process SPRING_TO_CRAM {

    tag "${meta.id}"

    /*
     * Keep whole-sample concurrency conservative. WGS alignment + sorting is
     * very I/O heavy. Increase only if the server/storage can sustain it.
     */
    maxForks 2

    cpus 12
    memory '40 GB'
    time '48h'

    /*
     * Same samtools/samblaster/sambamba environment used by germlineNGS.
     */
    conda "${params.sambamvcftools}"

    publishDir "${params.outdir}",
        mode: 'copy',
        overwrite: false,
        pattern: '*.cram*'

    input:
    tuple val(meta), path(spring_file)

    output:
    tuple val(meta),
        path("${meta.id}.${genome_version}.BWA.MD.cram"),
        path("${meta.id}.${genome_version}.BWA.MD.cram.crai"),
        emit: restored_crams

    script:
    """
    set -euo pipefail

    echo "============================================================"
    echo "Sample:      ${meta.id}"
    echo "SPRING:      ${spring_file}"
    echo "Reference:   ${genome_fasta}"
    echo "Genome ver.: ${genome_version}"
    echo "Host:        \$(hostname)"
    echo "CPUs:        ${task.cpus}"
    echo "Started:     \$(date)"
    echo "============================================================"

    if [[ ! -x "${params.spring_exec}" ]]; then
        echo "ERROR: SPRING executable not found or not executable:"
        echo "  ${params.spring_exec}"
        exit 1
    fi

    if [[ ! -f "${genome_fasta}" ]]; then
        echo "ERROR: reference FASTA not found:"
        echo "  ${genome_fasta}"
        exit 1
    fi

    mkdir -p spring_tmp
    mkdir -p local_tmp

    echo
    echo "[1/6] SPRING -> FASTQ.gz"
    echo "Started: \$(date)"

    "${params.spring_exec}" \
        -d \
        -i "${spring_file}" \
        -o "${meta.id}_R1.fastq.gz" "${meta.id}_R2.fastq.gz" \
        -g \
        -t ${task.cpus} \
        -w spring_tmp

    gzip -t "${meta.id}_R1.fastq.gz"
    gzip -t "${meta.id}_R2.fastq.gz"


    echo
    echo "[2/6] GATK FastqToSam"
    echo "Started: \$(date)"

    ${gatk_exec} FastqToSam \
        -F1 "${meta.id}_R1.fastq.gz" \
        -F2 "${meta.id}_R2.fastq.gz" \
        -SM "${meta.id}" \
        -PL illumina \
        -PU KGA_PU \
        -RG KGA_RG \
        --TMP_DIR local_tmp \
        -O "${meta.id}.unmapped.from.fq.bam"


    echo
    echo "[3/6] GATK MarkIlluminaAdapters"
    echo "Started: \$(date)"

    ${gatk_exec} MarkIlluminaAdapters \
        -I "${meta.id}.unmapped.from.fq.bam" \
        -O "${meta.id}.ubamXT.bam" \
        --TMP_DIR local_tmp \
        -M "${meta.id}.markAdapterMetrics.txt"


    echo
    echo "[4/6] SamToFastq -> BWA MEM -> MergeBamAlignment"
    echo "Started: \$(date)"

    ${gatk_exec} SamToFastq \
        -I "${meta.id}.ubamXT.bam" \
        -INTER \
        -CLIP_ATTR XT \
        -CLIP_ACT 2 \
        -NON_PF \
        -F /dev/stdout \
    | singularity run -B ${s_bind} ${simgpath}/bwa0717.sif bwa mem \
        -t ${task.cpus} \
        -p \
        "${genome_fasta}" \
        /dev/stdin \
    | ${gatk_exec} MergeBamAlignment \
        -R "${genome_fasta}" \
        -UNMAPPED "${meta.id}.ubamXT.bam" \
        -ALIGNED /dev/stdin \
        -MAX_GAPS -1 \
        -ORIENTATIONS FR \
        -SO queryname \
        --TMP_DIR local_tmp \
        -O "${meta.id}.${genome_version}.QNsort.BWA.clean.bam"


    echo
    echo "[5/6] Duplicate marking + coordinate sort -> CRAM"
    echo "Started: \$(date)"

    samtools view -h "${meta.id}.${genome_version}.QNsort.BWA.clean.bam" \
    | samblaster \
    | sambamba view \
        -t ${task.cpus} \
        -S \
        -f bam \
        /dev/stdin \
    | sambamba sort \
        -t ${task.cpus} \
        --tmpdir=local_tmp \
        -o /dev/stdout \
        /dev/stdin \
    | samtools view \
        -@ ${task.cpus} \
        -T "${genome_fasta}" \
        -C \
        -o "${meta.id}.${genome_version}.BWA.MD.cram" \
        -

    samtools index \
        -@ ${task.cpus} \
        "${meta.id}.${genome_version}.BWA.MD.cram"


    echo
    echo "[6/6] Validate CRAM"
    echo "Started: \$(date)"

    samtools quickcheck \
        -v \
        "${meta.id}.${genome_version}.BWA.MD.cram"

    echo
    echo "Final files:"
    ls -lh \
        "${meta.id}.${genome_version}.BWA.MD.cram" \
        "${meta.id}.${genome_version}.BWA.MD.cram.crai"


    echo
    echo "Cleaning temporary FASTQ/BAM files"

    rm -f \
        "${meta.id}_R1.fastq.gz" \
        "${meta.id}_R2.fastq.gz" \
        "${meta.id}.unmapped.from.fq.bam" \
        "${meta.id}.ubamXT.bam" \
        "${meta.id}.markAdapterMetrics.txt" \
        "${meta.id}.${genome_version}.QNsort.BWA.clean.bam"

    rm -rf spring_tmp local_tmp

    echo
    echo "Finished ${meta.id}: \$(date)"
    """
}


workflow {

    spring_samples = Channel.of(

        tuple(
            [id: "106378577899_WGS_Normal", npn: "106378577899"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/106378577899-WGS-Normal-48fterfgf-240515_LH00230_A227K5FLT4.spring", checkIfExists: true)
        ),

        tuple(
            [id: "107910640842_WGS_Normal", npn: "107910640842"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/107910640842-WGS-Normal-57gsphemf-240628_A01491_BHHJWFDSXC.spring", checkIfExists: true)
        ),

        tuple(
            [id: "112403936064_WGS_Normal", npn: "112403936064"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/112403936064-WGS-Normal-47byeedpf-230503_A01491_BH25K5DSX7.spring", checkIfExists: true)
        ),

        tuple(
            [id: "112428690530_WGS_Normal", npn: "112428690530"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/112428690530-WGS-Normal-85cavorcf-230519_A01491_AH327FDSX7.spring", checkIfExists: true)
        ),

        tuple(
            [id: "112479220183_WGS_Normal", npn: "112479220183"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/112479220183-WGS-Normal-59euseldf-240301_A00606_BH5K7NDSXC.spring", checkIfExists: true)
        ),

        tuple(
            [id: "112488074232_WGS_Normal", npn: "112488074232"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/112488074232-WGS-Normal-54deioldf-230911_A00606_AHK7M7DSX7.spring", checkIfExists: true)
        ),

        tuple(
            [id: "112567424356_WGS_Normal", npn: "112567424356"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/112567424356-WGS-Normal-85dwcowef-231117_LH00230_A222JH3LT4.spring", checkIfExists: true)
        ),

        tuple(
            [id: "112628532087_WGS_Normal", npn: "112628532087"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/112628532087-WGS-Normal-83eposwsf-240206_A01237_BH5FMTDSXC.spring", checkIfExists: true)
        ),

        tuple(
            [id: "D170078_WGS_Normal", npn: "D170078"],
            file("/lnx01_data2/shared/dataArchive/lnx04/springStorage/NGC_hjemtagning/D170078-WGS-Normal-60cytobrf-230821_A01237_AHHWK7DSX7.spring", checkIfExists: true)
        )
    )

    SPRING_TO_CRAM(spring_samples)
}


workflow.onComplete {
    println """
================================================================================
PRS WGS restoration finished

Success:  ${workflow.success}
Duration: ${workflow.duration}

Final CRAM + CRAI files:
${params.outdir}/
================================================================================
"""
}
