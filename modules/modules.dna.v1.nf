#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"

multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
//////////////////////////// SWITCHES ///////////////////////////////// 

switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    case 'v45':
    gatk_image="gatk4500.sif";
    default:
    gatk_image="gatk419.sif";
    break;
}


switch (params.server) {
    case 'lnx02':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/";
        simgpath="/data/shared/programmer/simg";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive/";
        refFilesDir="/fast/shared/genomes";
//        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/hg38v3/hg38v3_scatter20_BWI/*.interval_list";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/hg38v3/hg38v3_scatter10_IntervalSubdiv/*.interval_list";

        //modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
    break;
    case 'lnx01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga2/data/data.storage.archive/";
        modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
        refFilesDir="/data/shared/genomes";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/hg38v3/hg38v3_scatter10_IntervalSubdiv/*.interval_list";
        //params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";

    break;
    case 'kga01':
        simgpath="/data/shared/programmer/simg";
        s_bind="/data/:/data/";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive/";
        modules_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
        refFilesDir="/data/shared/genomes";
    break;
}

switch (params.genome) {
    case 'hg19':
        assembly="hg19"
        // Genome assembly files:
        genome_fasta = "/data/shared/genomes/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "/data/shared/genomes/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "/data/shared/genomes/hg19/human_g1k_v37.dict"
        genome_version="V1"
        break;


    case 'hg38':
        assembly="hg38"
        spliceai_assembly="grch38"
        smncaller_assembly="38"
        // Genome assembly files:
        if (params.hg38v1) {
        genome_fasta = "${refFilesDir}/hg38/GRCh38.primary.fa"
        genome_fasta_fai = "${refFilesDir}/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict = "${refFilesDir}/hg38/GRCh38.primary.dict"
        genome_version="V1"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_germline_PON/jgmr_45samples.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/"
        }
        
        if (params.hg38v2){
        genome_fasta = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.fa"
        genome_fasta_fai = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.fa.fai"
        genome_fasta_dict = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.dict"
        genome_version="V2"
        }

        // Current hg38 version (v3): NGC with masks and decoys.
        if (!params.hg38v2 && !params.hg38v1){
        genome_fasta = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.fa"
        genome_fasta_fai = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.fa.fai"
        genome_fasta_dict = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.dict"
        genome_version="V3"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/hg38v3_109samples.cnvkit.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared/genomes/hg38/inhouse_DBs/hg38v3/"
        }

        // Gene and transcript annotation files:

        gencode_gtf = "${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.gtf"
        gencode_gff3 = "${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.gff3"
     
        //Program  files:
        msisensor_list="${refFilesDir}/hg38/program_DBs/msisensor/hg38_msisensor_scan.txt"
        
      
        //Structural variants
        delly_exclude="/data/shared/genomes/hg38/program_DBs/delly/human.hg38.excl.tsv"
        
        smoove_exclude="/data/shared/genomes/hg38/interval.files/smoove/smoove.hg38.excluderegions.bed"
        smoove_gff="/data/shared/genomes/hg38/gene.annotations/GRCh38_latest_genomic.gff.gz"



        //Repeat Expansions:
        expansionhunter_catalog="/data/shared/genomes/hg38/program_DBs/expansionHunter/expansionHunter_hg38_stripy.variant_catalog.json"
        hipSTR_bed="/data/shared/genomes/hg38/interval.files/STRs/GRCh38.hipstr_reference.bed"

        // Somatic calling files (GATK Mutect2 pipeline):
        gatk_wgs_pon="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz"
        mutect_gnomad="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
        gatk_contamination_ref="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_small_exac_common_3.hg38.vcf.gz"

        // Program indexes:
        pcgr_assembly="grch38"
        sequenza_cg50_wig="/data/shared/genomes/hg38/program_DBs/sequenza/GRCh38.primary.cg50.sequenza.wig.gz"


        // Regions & variants:
        qualimap_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.6col.bed"
        gencode_exons_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"

        ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        
        //ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.bed"

        callable_regions="/data/shared/genomes/hg38/interval.files/GATK.hg38.callable.regions.bed"
        manta_callable_regions="/data/shared/genomes/hg38/interval.files/manta/GATK.hg38.callable.regions.bed.gz"

        dbsnp="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        KGindels="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
        KGindels_idx="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

        KGmills="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        KGmills_idx="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
        KG_p1_High_snps="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"

        hapmap="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
        omni="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
        AV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/av1.hg38.ROI.v2.bed"
        CV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV2_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV3_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV4_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv4.hg38.ROI.bed"
        CV5_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv5.hg38.ROI.bed"
        GV3_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/gv3.hg38.ROI.v2.bed"
        NV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/nv1.hg38.ROI.bed"
        WES_ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        MV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/muc1.hg38.coordinates.bed"
        break;
}

switch (params.panel) {
    case "AV1":
        ROI="${AV1_ROI}";
        panelID="AV1"
    break;

    case "CV5":
        ROI="${CV5_ROI}";
        panelID="CV5"
    break;

    case "GV3":
        ROI="${GV3_ROI}";
        panelID="GV3"
    break;

    case "GV_TEST":
        ROI="${GV3_ROI}";
        panelID="GV_TEST"
    break;

    case "MUC1":
        ROI="${MV1_ROI}";
        panelID="MUC1_MV1"
    break;

    case "WES_2":
        ROI="${WES_ROI}";
        panelID="WES"
    break;

    case "WES":
        ROI="${WES_ROI}";
        panelID="WES_subpanel"
    break;

    case "WGS_CNV":
        ROI="${WES_ROI}";
        panelID="WGS_CNV"
    break;
    
    default: 
        ROI="${WES_ROI}";
        panelID="WGS"
    break;
}


if (!params.archiveStorage) {
outputDir="${params.outdir}/"
}

if (params.archiveStorage) {
outputDir="${tank_storage}/alignedData/${params.genome}/${params.outdir}/"
}



channel
    .fromPath(params.intervals_list)
    .map { it -> tuple(it.baseName,it)}
    .set { haplotypecallerIntervalList }


log.info """\
===============================================
Clinical Genetics Vejle: Germline NGS v1
Panels or WGS analysis
===============================================
Genome       : $params.genome
Genome FASTA : $genome_fasta
ROI          : $ROI
AnalysisType : $params.panel
GATK ver.    : $gatk_image
Server       : $params.server
RunID        : $runID
PanelID      : $panelID
"""


/*********************************************
********** SYMLINKS (fastq or CRAM) **********
*********************************************/

process inputFiles_symlinks_fq{
    errorStrategy 'ignore'
    publishDir "${outputDir}/input_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    input:
    tuple val(sampleID), path(r1),path(r2)// from read_input2
    
    output:
    tuple path(r1),path(r2)
    script:
    """
    """
}



process inputFiles_symlinks_cram{
    errorStrategy 'ignore'
    publishDir "${outputDir}/input_symlinks/", mode: 'link', pattern: '*.{ba,cr}*'
    input:
    tuple val(sampleID), path(aln),path(index)// from symlink_input
    
    output:
    tuple path(aln),path(index)
    script:
    """
    """
}






///////////////////////////////// PREPROCESS MODULES //////////////////////// 
// input ch structure: As simple as possible: meta + actual data
process fastq_to_ubam {
    errorStrategy 'ignore'
    tag "$sampleID"
    //publishDir "${outputDir}/unmappedBAM/", mode: 'copy',pattern: '*.{bam,bai}'
    //publishDir "${outputDir}/fastq_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    cpus 20
    maxForks 10

    input:
    tuple val(sampleID), path(r1), path(r2)

    output:
    tuple val(sampleID), path("${sampleID}.unmapped.from.fq.bam")
  //  tuple path(r1),path(r2)
    
    script:
    """
    ${gatk_exec} FastqToSam \
    -F1 ${r1} \
    -F2 ${r2} \
    -SM ${sampleID} \
    -PL illumina \
    -PU KGA_PU \
    -RG KGA_RG \
    -O ${sampleID}.unmapped.from.fq.bam
    """
}

process markAdapters {

    input:
    tuple val(sampleID), path(uBAM)
    
    output:
    tuple val(sampleID), path("${sampleID}.ubamXT.bam"), path("${sampleID}.markAdapterMetrics.txt")
    
    script:

    """
    ${gatk_exec} MarkIlluminaAdapters \
    -I ${uBAM} \
    -O ${sampleID}.ubamXT.bam \
    --TMP_DIR ${tmpDIR} \
    -M ${sampleID}.markAdapterMetrics.txt
    """
}

process align {
    tag "$sampleID"

    maxForks 5
    errorStrategy 'ignore'
    cpus 20

    input:
    tuple val(sampleID), path(uBAM), path(metrics)

    output:
    tuple val(sampleID), path("${sampleID}.${params.genome}.${genome_version}.QNsort.BWA.clean.bam")
    
    script:
    """
    ${gatk_exec} SamToFastq \
    -I ${uBAM} \
    -INTER \
    -CLIP_ATTR XT \
    -CLIP_ACT 2 \
    -NON_PF \
    -F /dev/stdout \
    |  singularity run -B ${s_bind} ${simgpath}/bwa0717.sif bwa mem \
    -t ${task.cpus} \
    -p \
    ${genome_fasta} \
    /dev/stdin \
    | ${gatk_exec} MergeBamAlignment \
    -R ${genome_fasta} \
    -UNMAPPED ${uBAM} \
    -ALIGNED /dev/stdin \
    -MAX_GAPS -1 \
    -ORIENTATIONS FR \
    -SO queryname \
    -O ${sampleID}.${params.genome}.${genome_version}.QNsort.BWA.clean.bam
    """
}

process markDup_bam {
    errorStrategy 'ignore'
    maxForks 6
    tag "$sampleID"
    publishDir "${outputDir}/BAM/", mode: 'copy', pattern: "*.BWA.MD.ba*"
    publishDir "${outputDir}/CRAM/", mode: 'copy', pattern: "*.BWA.MD.cr*"
    
    input:
    tuple val(sampleID), path(aln) 
    
    output:
    tuple val(sampleID), path("${sampleID}.${params.genome}.${genome_version}.BWA.MD.bam"), path("${sampleID}.${params.genome}.${genome_version}.BWA.MD*bai")

    tuple val(sampleID), path("${sampleID}.${params.genome}.${genome_version}.BWA.MD.cram"), path("${sampleID}.${params.genome}.${genome_version}.BWA.MD*crai")
    
    script:
    """
    samtools view -h ${aln} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=/data/TMP/TMP.${user}/ -o ${sampleID}.${params.genome}.${genome_version}.BWA.MD.bam /dev/stdin
    sambamba index ${sampleID}.${params.genome}.${genome_version}.BWA.MD.bam
    
    samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${sampleID}.${params.genome}.${genome_version}.BWA.MD.cram ${sampleID}.${params.genome}.${genome_version}.BWA.MD.bam

    samtools index ${sampleID}.${params.genome}.${genome_version}.BWA.MD.cram
    """
}

process markDup_cram {
    errorStrategy 'ignore'
    maxForks 6
    tag "$sampleID"
    publishDir "${outputDir}/CRAM/", mode: 'copy', pattern: "*.BWA.MD.cr*"
    
    input:
    tuple val(sampleID), path(aln)
    
    output:
    tuple val(sampleID),  path("${sampleID}.${params.genome}.${genome_version}.BWA.MD.cram"), path("${sampleID}.${params.genome}.${genome_version}.BWA.MD*crai"), emit: markDup_output
    
    script:
    """
    samtools view -h ${aln} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=/data/TMP/TMP.${user}/ -o /dev/stdout /dev/stdin \
    |  samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${sampleID}.${params.genome}.${genome_version}.BWA.MD.cram -

    samtools index ${sampleID}.${params.genome}.${genome_version}.BWA.MD.cram
    """
}


// QC PROCESSES

process bamtools {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${outputDir}/QC/", mode: 'copy'

    input:
    val(sampleID),  path(aln), path(aln_index)
    output:
    path("${sampleID}.bamtools.sample.stats.txt"), emit: bamtools_out

    script:
    """
    bamtools stats -in ${aln} -insert > ${sampleID}.bamtools.sample.stats.txt
    """
}


process samtools {

    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${outputDir}/QC/", mode: 'copy'

    input:  
    tuple val(sampleID), path(aln), path(index)

    output:
    path("${sampleID}.samtools.sample.stats.txt"), emit: samtools

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/samtools.sif samtools \
    stats \
    ${aln} > ${sampleID}.samtools.sample.stats.txt
    """
}

// not working with CRAM:
process qualimap {
    errorStrategy 'ignore'
    tag "$sampleID"
    cpus 10
    maxForks 8
    publishDir "${outputDir}/QC/${sampleID}/qualimap/", mode: 'copy'

    input:
    tuple val(sampleID), path(aln), path(index)

    output:
    path ("${aln.baseName}/"), emit: qualimap_out
    //path ("*_results.txt") into bamQCReport

    script:
    use_bed = qualimap_ROI ? "-gff ${qualimap_ROI}" : ''
    """
    unset DISPLAY
    singularity run -B ${s_bind} ${simgpath}/qualimap.sif \
    qualimap --java-mem-size=5G bamqc \
    -nt ${task.cpus} \
    -outdir ${aln.baseName} \
    -bam ${aln} $use_bed -sd -sdmode 0 
    """
}

process fastqc_bam {
    errorStrategy 'ignore'
    tag "$sampleID"
    cpus 2
    publishDir "${outputDir}/QC/${sampleID}/", mode: 'copy'
    input:
    tuple val(sampleID), path(aln), path(index)
    
    output:
    path "*_fastqc.{zip,html}", emit: fastqc_bam

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/fastqc.sif --quiet --threads ${task.cpus} ${aln}
    """
}
// ^^^^^^^NOT WORKING WITH CRAM ^^^^^^ //

process collectWGSmetrics {

    errorStrategy 'ignore'
    tag "$sampleID"
    cpus 5
    publishDir "${outputDir}/QC/${sampleID}/picard/", mode: 'copy'

    input:
    tuple val(sampleID), path(aln), path(index)
    
    output:
    path("${sampleID}.picardWGSmetrics.txt"), emit: picard

    script:
    """
    ${gatk_exec} CollectWgsMetrics \
    -I ${aln} \
    -O ${sampleID}.picardWGSmetrics.txt \
    -R ${genome_fasta}
    """
}

process multiQC {
    
    errorStrategy 'ignore'
    publishDir "${outputDir}/QC/", mode: 'copy'

    input:
    path(inputfiles)
    //   path("_fastqc.*").collect().ifEmpty([])
    // path("${sampleID}.samtools.sample.stats.txt").collect().ifEmpty([])
    // path("bamQC/*").collect().ifEmpty([]) 
    //path("${sampleID}.picardWGSmetrics.txt").collect().ifEmpty([]) 

    output:
    path ("*multiqc_report.html")

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -c ${multiqc_config} \
    -f -q  ${launchDir}/*/QC/
    """
}

//////////////////////////// VARIANT CALLING MODULES //////////////////////////////////
process haplotypecaller{
        errorStrategy 'ignore'
        maxForks 30
        cpus 4
        tag "$sampleID"
        publishDir "${outputDir}/Variants/per_sample/", mode: 'copy', pattern: "*.HC.*"
        publishDir "${outputDir}/Variants/GVCF_files/", mode: 'copy', pattern: "*.g.*"
        publishDir "${outputDir}/HaplotypeCallerBAMout/", mode: 'copy', pattern: "*.HCbamout.*"
            
        input:
        tuple val(sampleID), path(aln), path(aln_index)
    
        output:
        path("${sampleID}.${params.genome}.${genome_version}.g.vcf"), emit: sample_gvcf

        tuple val(sampleID), path("${sampleID}.${params.genome}.${genome_version}.g.vcf"), emit: HC_sid_gvcf
    
        tuple val(sampleID), path("${sampleID}.${params.genome}.${genome_version}.HC.vcf"), path("${sampleID}.${params.genome}.${genome_version}.HC.vcf.idx"), emit: HC_sid_vcf

        path("${sampleID}.${params.genome}.${genome_version}.g.*")
        path("${sampleID}.${params.genome}.${genome_version}.HCbamout.*")
        path("${aln_index}")
        path("${aln}")
        //path("${sampleID_type}.HC.*")

        script:
        """
        ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller \
        -I ${aln} \
        -R ${genome_fasta} \
        -ERC GVCF \
        -L ${ROI} \
        --smith-waterman FASTEST_AVAILABLE \
        --native-pair-hmm-threads 4 \
        -pairHMM FASTEST_AVAILABLE \
        --dont-use-soft-clipped-bases \
        -O ${sampleID}.${params.genome}.${genome_version}.g.vcf \
        -bamout ${sampleID}.${params.genome}.${genome_version}.HCbamout.bam
    
        ${gatk_exec} GenotypeGVCFs \
        -R ${genome_fasta} \
        -V ${sampleID}.${params.genome}.${genome_version}.g.vcf \
        -O ${sampleID}.${params.genome}.${genome_version}.HC.vcf \
        -G StandardAnnotation \
        -G AS_StandardAnnotation
        """
}


process jointgenotyping {
        errorStrategy 'ignore'
        cpus 4
        publishDir "${outputDir}/Variants/", mode: 'copy', pattern: "*.VarSeq.*"
        publishDir "${outputDir}/Variants/GVCF_files/", mode: 'copy', pattern: "*.merged.g.*"
        //publishDir "tumorBoard_files", mode: 'copy', pattern: "*.VarSeq.*"

        input:
        tuple val(panelID), val(subpanel_gvcf) 
        //tuple val(sampleID),  path(vcf),path(idx) from joint_geno_dummy_ch
        output:
        //path("${params.rundir}.${params.panel}.merged.g.*") into merged_gVCF
        path("*.for.VarSeq.*")
        tuple val(panelID), path("${params.rundir}.${panelID}.${params.genome}.${genome_version}.merged.for.VarSeq.vcf"), emit: spliceAI_input

        script:
        """
        ${gatk_exec} --java-options "-Xmx64g" CombineGVCFs \
        -R ${genome_fasta} \
        ${subpanel_gvcf} \
        -O ${params.rundir}.${panelID}.${params.genome}.${genome_version}.merged.g.vcf \
        -L ${ROI} \
        -G StandardAnnotation -G AS_StandardAnnotation 

        ${gatk_exec} GenotypeGVCFs \
        -R ${genome_fasta} \
        -V ${params.rundir}.${panelID}.${params.genome}.${genome_version}.merged.g.vcf \
        -O ${params.rundir}.${panelID}.${params.genome}.${genome_version}.merged.for.VarSeq.vcf  \
        -L ${ROI} \
        -G StandardAnnotation -G AS_StandardAnnotation -A SampleList \
        -D ${dbsnp}
        """     
}

  process spliceAI {
        errorStrategy 'ignore'
        cpus 4
        publishDir "${outputDir}/Variants/", mode: 'copy', pattern: "*.spliceAI.merged.for.VarSeq.*"

        input:
        tuple val(panelID), path(vcf)// from spliceAI_input
        
        output:
        path("*.spliceAI.merged.for.VarSeq.*")// into merged_spliceAI_vcf

        when:
        !params.skipSpliceAI

        script:
        """
        singularity run -B ${s_bind} ${simgpath}/spliceai.sif spliceai \
        -R ${genome_fasta} \
        -I ${vcf} \
        -O ${params.rundir}.${panelID}.${params.genome}.${genome_version}.spliceAI.merged.for.VarSeq.vcf \
        -A ${spliceai_assembly}

        ${gatk_exec} IndexFeatureFile \
        -I ${params.rundir}.${panelID}.${params.genome}.${genome_version}.spliceAI.merged.for.VarSeq.vcf
        """
    }




//////// WGS VARIANT CALLING (HaplotypeCaller SplitIntervals)

process haplotypecallerSplitIntervals {
    errorStrategy 'ignore'
    maxForks 60

    input:
    tuple val(sampleID), path(bam), path(bai), val(sub_intID), path(sub_interval) //from HC_scatter_input_bam.combine(interval_list1)

    output:
    tuple val(sampleID), path("${sampleID}.${sub_intID}.g.vcf"), path("${sampleID}.${sub_intID}.g.vcf.idx"), emit: hc_split_output

    script:
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller \
    -I ${bam} \
    -R ${genome_fasta} \
    -ERC GVCF \
    -L ${sub_interval} \
    --smith-waterman FASTEST_AVAILABLE \
    --native-pair-hmm-threads 4 \
    -pairHMM FASTEST_AVAILABLE \
    --dont-use-soft-clipped-bases \
    -O ${sampleID}.${sub_intID}.g.vcf
    """
}

process mergeScatteredGVCF{
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${outputDir}/Variants/", mode: 'copy'
    maxForks 9

    input:

    tuple val(sampleID), path(sub_gvcf), path(sub_gvcf_idx)// from hc_split_output.groupTuple()
    
    output:
    path("${sampleID}.*")
    tuple val(sampleID), path("${sampleID}.HC.vcf"),path("${sampleID}.HC*.idx"), emit: hc_singlesamplevcf_ch1

    path("${sampleID}.g.vcf"), emit: sample_gvcf_list_scatter
    path("*.WES_ROI.*")
    script:
    """
    ${gatk_exec} SortVcf \
    ${sub_gvcf.collect { "--INPUT $it " }.join()} \
    --OUTPUT ${sampleID}.g.vcf
    
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=30" GenotypeGVCFs \
    -R ${genome_fasta} \
    -V ${sampleID}.g.vcf \
    -O ${sampleID}.HC.vcf 

    ${gatk_exec} SelectVariants \
    -R ${genome_fasta} \
    -V ${sampleID}.HC.vcf \
    -L ${ROI} \
    -O ${sampleID}.WES_ROI.vcf

    ${gatk_exec} IndexFeatureFile \
    -I ${sampleID}.WES_ROI.vcf
    """
}


process jointgenoScatter{
    errorStrategy 'ignore'
    publishDir "${outputDir}/Variants/", mode: 'copy'

    input:
    val x //from gvcfsamples_for_GATK_scatter

    output:
    path("${params.rundir}.merged.g.{vcf,idx}") //into merged_gVCF_scatter
    path("${params.rundir}.merged.RAW.{vcf,vcf.idx}")// into merged_RAW_vcf_scatter
    path("${params.rundir}.merged.WES_ROI.*")
    
    when:
    !params.single

    script:
    """
    ${gatk_exec} CombineGVCFs \
    -R ${genome_fasta} ${x} \
    -O ${params.rundir}.merged.g.vcf \
    -G StandardAnnotation -G AS_StandardAnnotation 

    ${gatk_exec} GenotypeGVCFs \
    -R ${genome_fasta} \
    -V ${params.rundir}.merged.g.vcf \
    -O ${params.rundir}.merged.RAW.vcf  \
    -G StandardAnnotation -G AS_StandardAnnotation -A SampleList -D ${dbsnp}
    
    ${gatk_exec} SelectVariants \
    -R ${genome_fasta} \
    -V ${params.rundir}.merged.RAW.vcf \
    -L ${ROI} \
    -O ${params.rundir}.merged.WES_ROI.vcf

    ${gatk_exec} IndexFeatureFile \
    -I ${params.rundir}.merged.WES_ROI.vcf
    """     
}







/////////////////////////////// SV CALLING MODULES //////////////////////

process manta {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${inhouse_SV}/manta/raw_calls/", mode: 'copy', pattern: " ${sampleID}.manta.diploidSV.*"
    publishDir "${outputDir}/structuralVariants/manta/allOutput/", mode: 'copy'

    cpus 10
    maxForks 6

    input:
    tuple val(sampleID), path(aln), path(index)

    output:

    path("${sampleID}.manta.*.{vcf,vcf.gz,gz.tbi}")
    
   //  tuple val(sampleID), path("${sampleID}.manta.INVconverted.vcf"), emit: manta
    tuple val(sampleID), path("${sampleID}.manta.diploidSV.vcf.gz"), path("${sampleID}.manta.diploidSV.vcf.gz.tbi"), emit: manta
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/manta1.6_strelka2.9.10.sif configManta.py \
    --bam ${aln} \
    --referenceFasta ${genome_fasta} \
    --callRegions ${manta_callable_regions} \
    --runDir manta

    singularity run -B ${s_bind} ${simgpath}/manta1.6_strelka2.9.10.sif ./manta/runWorkflow.py -j ${task.cpus}

    mv manta/results/variants/candidateSmallIndels.vcf.gz \
    ${sampleID}.manta.candidateSmallIndels.vcf.gz
    
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    ${sampleID}.manta.candidateSmallIndels.vcf.gz.tbi
    
    mv manta/results/variants/candidateSV.vcf.gz \
    ${sampleID}.manta.candidateSV.vcf.gz
    
    mv manta/results/variants/candidateSV.vcf.gz.tbi \
    ${sampleID}.manta.candidateSV.vcf.gz.tbi

    mv manta/results/variants/diploidSV.vcf.gz \
    ${sampleID}.manta.diploidSV.vcf.gz
    
    mv manta/results/variants/diploidSV.vcf.gz.tbi \
    ${sampleID}.manta.diploidSV.vcf.gz.tbi

    gzip -dc ${sampleID}.manta.diploidSV.vcf.gz > ${sampleID}.manta.diploidSV.vcf
    

    """
}

process filter_manta {
    tag "$sampleID"
    errorStrategy 'ignore'

    publishDir "${inhouse_SV}/manta/filtered/", mode: 'copy', pattern: "*.filtered.vcf"
    publishDir "${outputDir}/structuralVariants/manta/", mode: 'copy', pattern: "*.filtered.vcf"
    
    input:
    tuple val(sampleID), path(vcf), path(idx)

    output:
    tuple val(sampleID), path("${sampleID}.manta.filtered.vcf"), emit: mantaForSVDB
    
    script:
    """
    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${vcf} \
    --exclude-filtered \
    -O ${sampleID}.manta.filtered.vcf
    """
}

process lumpy {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${inhouse_SV}/lumpy/raw_calls/", mode: 'copy', pattern: "*.Lumpy_altmode_step1.vcf"
    publishDir "${outputDir}/structuralVariants/lumpy/", mode: 'copy'
    
    cpus 10
    maxForks 6

    input:

    tuple val(sampleID), path(aln), path(index)
    //    path(genome_fasta)
    //    path(genome_fasta_fai)
    //    path(genome_fasta_dict)

    output:
    tuple val(sampleID), path("${sampleID}.Lumpy_altmode_step1.vcf"), emit: lumpyForSVDB
    path("*.Lumpy_altmode_step1.vcf.gz") 

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/smoove.sif smoove call -d \
    --outdir ${params.rundir}.LumpyAltSingle \
    --exclude ${smoove_exclude} \
    --name ${sampleID} \
    --fasta ${genome_fasta} \
    -p ${task.cpus} \
    --genotype ${aln}
    
    mv ${params.rundir}.LumpyAltSingle/${sampleID}*.genotyped.vcf.gz \
    ${sampleID}.Lumpy_altmode_step1.vcf.gz

    gzip -dc  ${sampleID}.Lumpy_altmode_step1.vcf.gz >  ${sampleID}.Lumpy_altmode_step1.vcf

    mv ${params.rundir}.LumpyAltSingle/${sampleID}*.csi \
    ${sampleID}.Lumpy_altmode_step1.vcf.gz.csi
    """
}

process tiddit361 {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${inhouse_SV}/tiddit/PASSED_calls/", mode: 'copy', pattern: "*.PASS.vcf"
    publishDir "${inhouse_SV}/tiddit/RAW_calls/", mode: 'copy', pattern: "*.tiddit.vcf"
    publishDir "${outputDir}/structuralVariants/tiddit/", mode: 'copy'
    
    cpus 2
    maxForks 15

    input:
    tuple val(sampleID), path(aln), path(index) 

    output:
    tuple val(sampleID), path("*.{vcf,tab}"), emit: tiddit_out_ch
    tuple val(sampleID), path("*.PASS.vcf"), emit: tidditForSVDB
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/tiddit361.sif tiddit \
    --sv \
    --bam ${aln} \
    --threads ${task.cpus} \
    -q 10 \
    --ref ${genome_fasta} \
    -o ${sampleID}_tiddit

    cat ${sampleID}_tiddit.vcf | grep -E "#|PASS" > ${sampleID}_tiddit.PASS.vcf
    """
}

process cnvkit {
    errorStrategy 'ignore'
    tag "$sampleID"

    cpus 10
    maxForks 6

    publishDir "${outputDir}/structuralVariants/cnvkit/", mode: 'copy'
    publishDir "${inhouse_SV}/CNVkit/CNNfiles/", mode: 'copy', pattern: '*.cnn'

    input:
    tuple val(sampleID), path(aln), path(index)

    output:
    path("${sampleID}.cnvkit/*")
    path("*.targetcoverage.cnn"), emit: cnvkit_cnn_out
    tuple val(sampleID), path("${sampleID}.cnvkit/*.call.cns"), emit: CNVcalls
    tuple val(sampleID), path("${sampleID}.cnvkit/*.cnr"), emit: CNVcnr
    //path("${sampleID}.cnvkit/*.cnn")
    
    // touch ${index}
    script:
    """
    mv ${index} intermediate_crai
    cp intermediate_crai ${index}
    rm intermediate_crai
    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py batch \
    ${aln} \
    -m wgs \
    -p ${task.cpus} \
    -r ${cnvkit_germline_reference_PON} \
    --scatter --diagram \
    -d ${sampleID}.cnvkit/
    cp ${sampleID}.cnvkit/*.targetcoverage.cnn .
    """
}

process cnvkitExportFiles {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${inhouse_SV}/CNVkit/raw_calls/", mode: 'copy', pattern: '*.cnvkit.vcf'
    publishDir "${outputDir}/structuralVariants/cnvkit/", mode: 'copy'

    input:
    tuple val(sampleID), path(cnvkit_calls)// from cnvkit_calls_out
    tuple val(sampleID), path(cnvkit_cnr)// from cnvkit_cnr_out

    output:
    path("*.vcf")
    path("*.seg")
    tuple val(sampleID), path("${sampleID}.cnvkit.vcf"), emit: cnvkitForSVDB

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py export vcf \
    ${cnvkit_calls} \
    -i ${sampleID} \
    -o ${sampleID}.cnvkit.vcf

    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py export seg \
    ${cnvkit_cnr} \
    -o ${sampleID}.cnvkit.cnr.seg
    """
}

process merge4callerSVDB {
    tag "$sampleID"
    errorStrategy 'ignore'

    //publishDir "${outputDir}/all_callers_merged/", mode: 'copy'
   // publishDir "${outputDir}/structuralVariants/SVDB_merged/", mode: 'copy', pattern: "*.4caller.SVDB.merged.*"
    publishDir "${outputDir}/structuralVariants/SVDB_merged/60pctOverlap/", mode: 'copy', pattern: "*.60pctOverlap.*"
    publishDir "${outputDir}/structuralVariants/SVDB_merged/80pctOverlap/", mode: 'copy', pattern: "*.80pctOverlap.*"
    publishDir "${outputDir}/structuralVariants/SVDB_merged/100pctOverlap/", mode: 'copy', pattern: "*.100pctOverlap.*"

    //publishDir "${outputDir}/", mode: 'copy', pattern: '*.vcf'
    //container 'kfdrc/manta:1.6.0'
    maxForks 12
    input:
    tuple val(sampleID), path(manta_vcf), path(lumpy_vcf),path(cnvkit_vcf),path(tiddit_vcf) // from single_4caller_for_svdb
    output:
    path("${sampleID}.4caller.SVDB.merged.*")

    script:
    """
    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --merge \
    --overlap 0.6 \
    --vcf ${manta_vcf}:MANTA ${lumpy_vcf}:LUMPY ${cnvkit_vcf}:CNVKIT ${tiddit_vcf}:TIDDIT \
    --priority LUMPY,MANTA,CNVKIT,TIDDIT > ${sampleID}.4caller.SVDB.merged.60pctOverlap.vcf

    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --merge \
    --overlap 0.8 \
    --vcf ${manta_vcf}:MANTA ${lumpy_vcf}:LUMPY ${cnvkit_vcf}:CNVKIT ${tiddit_vcf}:TIDDIT \
    --priority LUMPY,MANTA,CNVKIT,TIDDIT > ${sampleID}.4caller.SVDB.merged.80pctOverlap.vcf


    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --merge \
    --overlap 1.0 \
    --vcf ${manta_vcf}:MANTA ${lumpy_vcf}:LUMPY ${cnvkit_vcf}:CNVKIT ${tiddit_vcf}:TIDDIT \
    --priority LUMPY,MANTA,CNVKIT,TIDDIT > ${sampleID}.4caller.SVDB.merged.100pctOverlap.vcf
    """
}

process expansionHunter {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${outputDir}/repeatExpansions/expansionHunter/", mode: 'copy'
    cpus 10
    input:
    tuple val(sampleID), path(aln), path(index)

    output:
    path("*.{vcf,json,bam}")
    script:
    """
    /data/shared/programmer/BIN/ExpansionHunter500 \
    --reads ${aln} \
    --reference ${genome_fasta} \
    --variant-catalog ${expansionhunter_catalog} \
    -n ${task.cpus} \
    --output-prefix ${sampleID}.expansionHunter
    """
}

process stripy {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${outputDir}/repeatExpansions/STRipy/", mode: 'copy'

    input:
    tuple val(sampleID), path(aln), path(index)

    output:
    path("${sampleID}.stripy/*.html")
    script:
    """
    mkdir ${sampleID}.stripy/ 
    sleep 5
    python3 /data/shared/programmer/stripy-pipeline-main/stri.py \
    --genome ${params.genome} \
    --reference ${genome_fasta} \
    --locus AFF2,AR,ARX_1,ARX_2,ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,ATXN8OS,BEAN1,C9ORF72,CACNA1A,CBL,CNBP,COMP,DAB1,DIP2B,DMD,DMPK,FGF14,FMR1,FOXL2,FXN,GIPC1,GLS,HOXA13_1,HOXA13_2,HOXA13_3,HOXD13,HTT,JPH3,LRP12,MARCHF6,NIPA1,NOP56,NOTCH2NLC,NUTM2B-AS1,PABPN1,PHOX2B,PPP2R2B,PRDM12,RAPGEF2,RFC1,RILPL1,RUNX2,SAMD12,SOX3,STARD7,TBP,TBX1,TCF4,TNRC6A,XYLT1,YEATS2,ZIC2,ZIC3 \
    --output ${sampleID}.stripy/ \
    --input ${aln}
    """
}

process prepareManifestSMN {
    
    publishDir "${outputDir}/SMNcaller/", mode: 'copy'
    
    input:
    path(samplesheet) // from smn_input_ch
    
    output:
    path("SMNmanifest.txt")//, emit: SMN_manifest_ch

    shell:
    '''
    cat !{samplesheet} | cut -f2 > SMNmanifest.txt
    '''
}

process smnCopyNumberCaller {
    publishDir "${outputDir}/SMNcaller/", mode: 'copy'
    errorStrategy "ignore"
    cpus 12

    input:
    path(manifest)// from SMN_manifest_ch

    output:
    path("*.{tsv,pdf,json}")
    
    script:
    """
    python3 /data/shared/programmer/SMNCopyNumberCaller-1.1.2/smn_caller.py \
    --manifest ${manifest} \
    --genome ${smncaller_assembly} \
    --prefix ${params.rundir} \
    --threads ${task.cpus} \
    --outDir .

    python /data/shared/programmer/SMNCopyNumberCaller-1.1.2/smn_charts.py \
    -s ${params.rundir}.json \
    -o .
    """
}



/////////////////////////////////////////////////////////////
/// SUBWORKFLOWS meta r1 r2 input channel///////
/////////////////////////////////////////////////////////////

workflow SUB_PREPROCESS {

    take:
    fq_read_input
    
    main:
    inputFiles_symlinks_fq(fq_read_input)
    fastq_to_ubam(fq_read_input)
    markAdapters(fastq_to_ubam.out[0])
    align(markAdapters.out)
    markDup_cram(align.out)
    //markDup_v3_cram.out.markDup_output.view()
    emit:
    finalAln=markDup_cram.out.markDup_output
}


/////////////////////////////////////////////////////////////
/// SUBWORKFLOWS meta-aln-index input channel///////
/////////////////////////////////////////////////////////////

workflow SUB_VARIANTCALL {
    take:
    meta_aln_index  // sampleID, aln, index
    main:
    haplotypecaller(meta_aln_index)
    haplotypecaller.out.sample_gvcf
    .map{ tuple(it.simpleName, it) }
    .set { gvcf_list }
    if (panelID=="AV1"){
        gvcf_list
            .filter {it =~/_CV6/}
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileTEST_CV6.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple("CV6",it) }
            .set { cv6_gatk }

        gvcf_list
            .filter {it =~/_FSGS1/}
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileTEST_FSGS1.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple("FSGS1",it) }
            .set {fsgs1_gatk}

        gvcf_list
            .filter {it =~/_NV2/}
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileTEST_NV1.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple("NV2",it) }
            .set {nv2_gatk}

        gvcf_list
            .filter {it =~/_FH1/}
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileTEST_FH1.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple("FH1",it) }
            .set {fh1_gatk}

        gvcf_list
            .filter {it =~/_GV4/}
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileTEST_GV3.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple("GV4",it) }
            .set {gv4_gatk}

        gvcf_list
            .filter {it =~/_OBS/}
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileTEST_OBS.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple("OBS",it) }
            .set {obs_gatk}

        gvcf_list
            .filter {it =~/AV1/}
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileTEST_AV1.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple("AV1_ALL",it) }
            .set {av1_gatk}

        cv6_gatk.concat(fsgs1_gatk).concat(nv2_gatk).concat(fh1_gatk).concat(gv4_gatk).concat(av1_gatk).concat(obs_gatk)
                .set { gvcfsamples_for_GATK }

    }

    if (panelID=="WES_subpanel"){
        gvcf_list
                .filter {it =~/_ONK/}
                .map{" -V "+ it[1] }
                .collectFile(name: "collectfileTEST_EV8_ONK.txt", newLine: false)
                .map {it.text.trim()}
                .map { tuple("EV8_ONK",it) }
                .set {onk_ev8_gatk}

        gvcf_list
                .filter {it =~/_ALM/}
                .map{" -V "+ it[1] }
                .collectFile(name: "collectfileTEST_EV8_ALM.txt", newLine: false)
                .map {it.text.trim()}
                .map { tuple("EV8_ALM",it) }
                .set {alm_ev8_gatk}

        onk_ev8_gatk.concat(alm_ev8_gatk)
                .set {gvcfsamples_for_GATK}
    }

    if (panelID != "AV1" && panelID!= "WES_subpanel") {
        gvcf_list
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileNOTAV1.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple(panelID, it) }
            .set {gvcfsamples_for_GATK}
    }

    jointgenotyping(gvcfsamples_for_GATK)
    
    if (panelID=="AV1" && params.server!="kga01"){
        spliceAI(jointgenotyping.out.spliceAI_input)
    }
}

workflow SUB_VARIANTCALL_WGS {
    take:
    meta_aln_index
    main:
    haplotypecallerSplitIntervals(meta_aln_index.combine(haplotypecallerIntervalList))
    mergeScatteredGVCF(haplotypecallerSplitIntervals.out)
    
    mergeScatteredGVCF.out.sample_gvcf_list_scatter
    .map{" -V "+ it }
    .set{gvcflist_scatter_done}

    gvcflist_scatter_done
    .collectFile(name: "collectfileTEST_scatter.txt", newLine: false)
    .map {it.text.trim()}.set {gvcfsamples_for_GATK_scatter}

    //jointgenoScatter(gvcfsamples_for_GATK_scatter)
}

workflow SUB_CNV_SV {
    take:
    meta_aln_index
    main:
    manta(meta_aln_index)
    filter_manta(manta.out.manta)   // mantafiltered for SVDB
    lumpy(meta_aln_index)
    cnvkit(meta_aln_index)
    cnvkitExportFiles(cnvkit.out.CNVcalls, cnvkit.out.CNVcnr)
    tiddit361(meta_aln_index)
    merge4callerSVDB(filter_manta.out.mantaForSVDB.join(lumpy.out.lumpyForSVDB).join(cnvkitExportFiles.out.cnvkitForSVDB).join(tiddit361.out.tidditForSVDB))
}

workflow SUB_STR {
    take:
    meta_aln_index
    main:
    expansionHunter(meta_aln_index) 
    stripy(meta_aln_index)
}

workflow SUB_SMN {
    take:
    meta_aln_index
    main:
    
    meta_aln_index
    .map {"TEST"+'\t'+it[1]}
    .collectFile(name: "smncaller_manifest.txt", newLine: true, storeDir: "${launchDir}/")
    .set{smn_input_ch}
    
    prepareManifestSMN(smn_input_ch)
    smnCopyNumberCaller(prepareManifestSMN.out)
}
