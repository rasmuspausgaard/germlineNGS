#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
date2=new Date().format( 'yyMMdd HH:mm:ss' )
user="$USER"
runID="${date}.${user}"


// Server and reference independent variables
multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
vntyperREF="/data/shared/genomes/hg19/program_DBs/vntyper"
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
        if (params.panel=="AV1" || params.panel=="GV3" || params.panel=="CV5"){
            gatk_image="gatk419.sif";
        }
        else {
            gatk_image="gatk4400.sif";
        }
    break;
}


switch (params.server) {

    case 'lnx01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        dataStorage="/lnx01_data3/storage/";
        refFilesDir="/data/shared/genomes";
    break;

    default:
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/,/lnx01_data4/:/lnx01_data4/";
        simgpath="/data/shared/programmer/simg";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        dataStorage="/lnx01_data3/storage/";
        dataArchive="/lnx01_data2/shared/dataArchive";
        refFilesDir="/fast/shared/genomes";
    break;
}

switch (params.genome) {
    case 'hg19':
        assembly                   ="hg19"
        // Genome assembly files:
      //  genome_fasta            ="/data/shared/genomes/hg19/human_g1k_v37.fasta"
      //  genome_fasta_fai        ="/data/shared/genomes/hg19/human_g1k_v37.fasta.fai"
      //  genome_fasta_dict       ="/data/shared/genomes/hg19/human_g1k_v37.dict"
            genome_fasta           = "/data/shared/genomes/hg19/pacbio_hg19/human_hs37d5.fasta"
            genome_fasta_fai       = "/data/shared/genomes/hg19/pacbio_hg19/human_hs37d5.fasta.fai"
            genome_fasta_dict      = "/data/shared/genomes/hg19/pacbio_hg19/human_hs37d5.fasta.dict"
        genome_version             = "hg19v2"
            ROI                    ="/data/shared/genomes/hg19/interval.files/200108.NCBIrefseq.codingexons.nocontig.20bp.merged.sorted.bed"
            WES_ROI                ="/data/shared/genomes/hg19/interval.files/200108.NCBIrefseq.codingexons.nocontig.20bp.merged.sorted.bed"
        smoove_exclude          ="/data/shared/genomes/hg19/interval.files/smoove.hg19.excluderegions.bed"

        //sv databases wont work, but output not used. Required to get single sample SV calls:
        svdb_databases          ="/data/shared/genomes/hg38/inhouse_DBs/hg38v3/svdb_AF"    
        mantaSVDB               ="${svdb_databases}/mantaSVDB315.db"
        lumpySVDB               ="${svdb_databases}/lumpySVDB218.db"
        cnvkitSVDB              ="${svdb_databases}/cnvkitSVDB313.db"
        dellySVDB               ="${svdb_databases}/dellySVDB112.db"
    
    break;


    case 'hg38':
        assembly                ="hg38"
        spliceai_assembly       ="grch38"
        smncaller_assembly      ="38"
        svdb_databases          ="/data/shared/genomes/hg38/inhouse_DBs/hg38v3/svdb_AF"
        // Genome assembly files:
        if (params.hg38v1) {
        genome_fasta            ="${refFilesDir}/hg38/GRCh38.primary.fa"
        genome_fasta_fai        ="${refFilesDir}/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict       ="${refFilesDir}/hg38/GRCh38.primary.dict"
        genome_version          ="hg38v1"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_germline_PON/jgmr_45samples.reference.cnn"
        cnvkit_inhouse_cnn_dir  ="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV              ="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/"
        }
        
        if (params.hg38v2){
        genome_fasta            ="${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.fa"
        genome_fasta_fai        ="${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.fa.fai"
        genome_fasta_dict       ="${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.dict"
        genome_version          ="hg38v2"
        }

        // Current hg38 version (v3): NGC with masks and decoys.
        if (!params.hg38v2 && !params.hg38v1){
        genome_fasta            ="${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.fa"
        genome_fasta_fai        ="${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.fa.fai"
        genome_fasta_dict       ="${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.dict"
        genome_version          ="hg38v3"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/hg38v3_109samples.cnvkit.reference.cnn"
        cnvkit_inhouse_cnn_dir  ="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/wgs_persample_cnn/"
    //  inhouse_SV              ="/data/shared/genomes/hg38/inhouse_DBs/hg38v3/"
        inhouse_SV              ="/lnx01_data4/archive/variantStorage/hg38/KGVejle/structuralVariants/"
        }

        // Gene and transcript annotation files:

        gencode_gtf             = "${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.gtf"
        gencode_gff3            = "${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.gff3"
     
        //Program  files:
        msisensor_list          ="${refFilesDir}/hg38/program_DBs/msisensor/hg38_msisensor_scan.txt"
        intervalList_GATK       ="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/hg38v3/hg38v3_scatter20_BWI/*.interval_list";             
      
        //Structural variants
        delly_exclude           ="/data/shared/genomes/hg38/program_DBs/delly/human.hg38.excl.tsv"
        smoove_exclude          ="/data/shared/genomes/hg38/interval.files/smoove/smoove.hg38.excluderegions.bed"
        smoove_gff              ="/data/shared/genomes/hg38/gene.annotations/GRCh38_latest_genomic.gff.gz"


        //inhouse SV AF databases: 
        mantaSVDB               ="${svdb_databases}/mantaSVDB315.db"
        lumpySVDB               ="${svdb_databases}/lumpySVDB218.db"
        cnvkitSVDB              ="${svdb_databases}/cnvkitSVDB313.db"
        //tidditSVDB="${svdb_databases}/tidditSVDB.db"
        dellySVDB               ="${svdb_databases}/dellySVDB112.db"


        //Repeat Expansions:
        expansionhunter_catalog ="/data/shared/genomes/hg38/program_DBs/expansionHunter/expansionHunter_hg38_stripy.variant_catalog.json"
        hipSTR_bed              ="/data/shared/genomes/hg38/interval.files/STRs/GRCh38.hipstr_reference.bed"

        // Somatic calling files (GATK Mutect2 pipeline):
        gatk_wgs_pon            ="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz"
        mutect_gnomad           ="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
        gatk_contamination_ref  ="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_small_exac_common_3.hg38.vcf.gz"

        // Program indexes:
        pcgr_assembly           ="grch38"
        sequenza_cg50_wig       ="/data/shared/genomes/hg38/program_DBs/sequenza/GRCh38.primary.cg50.sequenza.wig.gz"


        // Regions & variants:
        qualimap_ROI            ="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.6col.bed"
        gencode_exons_ROI       ="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"

        ROI                     ="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        
        //ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.bed"

        callable_regions        ="/data/shared/genomes/hg38/interval.files/GATK.hg38.callable.regions.bed"
        manta_callable_regions  ="/data/shared/genomes/hg38/interval.files/manta/GATK.hg38.callable.regions.bed.gz"

        dbsnp                   ="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        KGindels                ="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
        KGindels_idx            ="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

        KGmills                 ="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        KGmills_idx             ="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
        KG_p1_High_snps         ="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"

        hapmap                  ="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
        omni                    ="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
        AV1_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/av1.hg38.ROI.v2.bed"
        CV1_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV2_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV3_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV4_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/cv4.hg38.ROI.bed"
        CV5_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/cv5.hg38.ROI.bed"
        GV3_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/gv3.hg38.ROI.v2.bed"
        NV1_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/nv1.hg38.ROI.bed"
        WES_ROI                 ="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        MV1_ROI                 ="/data/shared/genomes/${params.genome}/interval.files/panels/muc1.hg38.coordinates.bed"
    break;

    case 't2t':
        assembly                   = "t2t"
        genome_mmi                 = "/data/shared/genomes/t2t/pacbio_t2t/chm13v2p0_maskedY_rCRS.mmi"
        genome_fasta               = "/data/shared/genomes/t2t/pacbio_t2t/chm13v2p0_maskedY_rCRS.fasta"
        genome_fasta_fai           = "/data/shared/genomes/t2t/pacbio_t2t/chm13v2p0_maskedY_rCRS.fasta.fai"
        genome_fasta_dict          = "/data/shared/genomes/t2t/pacbio_t2t/chm13v2p0_maskedY_rCRS.dict"
        genome_version             = "T2Tv1"
        intervalList_GATK          ="/data/shared/genomes/t2t/interval.files/WGS_splitIntervals/t2t_15intervals_BWI_OVERFLOW/*.interval_list"

        pbSV_trf                   = "/data/shared/genomes/t2t/pacbio_t2t/chm13v2p0_maskedY_rCRS.trf.bed"
        ROI                        = "/data/shared/genomes/t2t/interval.files/250314.T2T.RefSeq.fullExons+50bp.bed"
        WES_ROI                    = "/data/shared/genomes/t2t/interval.files/250314.T2T.RefSeq.fullExons+50bp.bed"

    break;

}

switch (params.panel) {
    case "AV1":
        ROI="${AV1_ROI}";
        panelID="AV1";
        panelID_storage="AV1"
    break;

    case "CV5":
        ROI="${CV5_ROI}";
        panelID="CV5";
        panelID_storage="deprecated_panels"
    break;

    case "GV3":
        ROI="${GV3_ROI}";
        panelID="GV3";
        panelID_storage="deprecated_panels"
    break;

    case "GV_TEST":
        ROI="${GV3_ROI}";
        panelID="GV_TEST";
        panelID_storage="deprecated_panels"
    break;

    case "MV1":
        ROI="${MV1_ROI}";
        panelID="MV1"
        panelID_storage="MV1"
    break;

    case "WES_2":
        ROI="${WES_ROI}";
        panelID="WES";
        panelID_storage="WES"
    break;

    case "WES":
        ROI="${WES_ROI}";
        panelID="WES_subpanel";
        panelID_storage="WES"
    break;

    case "WGS_CNV":
        ROI="${WES_ROI}";
        panelID="WGS_CNV";
        panelID_storage="WGS"
    break;

    case "WGS_NGC":
        ROI="${WES_ROI}";
        panelID="WGS_NGC";
        panelID_storage="WGS"
    break;

    case "WGS_ALL":
        ROI="${WES_ROI}";
        panelID="WGS_ALL";
        panelID_storage="WGS"
    break;

    default: 
        ROI="${WES_ROI}";
        panelID="WGS";
        panelID_storage="undetermined"
    break;
}



outputDir="${params.outdir}/"
variantStorage="${dataStorage}/variantStorage/${params.genome}/"
cramStorage="${dataStorage}/alignedData/${params.genome}/"



channel
    .fromPath(intervalList_GATK)
    .map { it -> tuple(it.baseName,it)}
    .set { haplotypecallerIntervalList }
/*
*/

log.info """\
======================================================
Clinical Genetics Vejle: GermlineNGS v3 + SLURM
======================================================
Genome       : $params.genome
Genome FASTA : $genome_fasta
ROI          : $ROI
AnalysisType : $params.panel
GATK ver.    : $gatk_image
Server       : $params.server
RunID        : $runID
PanelID      : $panelID

Script start : $date2
"""
//IntervalList : $intervalList_GATK

/*********************************************
********** SYMLINKS (fastq or CRAM) **********
*********************************************/

process inputFiles_symlinks_fq{
    label 'low'
    
    publishDir "${outputDir}/input_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    
    input:
    tuple val(meta), path(reads)// from read_input2
    
    output:
    tuple val(meta), path(reads)
    
    script:
    """
    """
}



process inputFiles_symlinks_cram{
    label 'low'
   
    publishDir "${outputDir}/input_symlinks/", mode: 'link', pattern: '*.{ba,cr}*'
    publishDir "${outputDir}/Variants/CRAM_symlinks/", mode: 'link', pattern: '*.{ba,cr}*'
   
    input:
    tuple val(meta), path(aln)// from symlink_input
    
    output:
    tuple val(meta), path(aln)
    
    script:
    """
    """
}


process inputFiles_symlinks_spring{
    label 'low'

    publishDir "${outputDir}/input_symlinks/", mode: 'link', pattern: '*.spring'

    input:
    tuple val(meta), path(spring)    
    
    output:
    path(spring)
    
    script:
    """
    """
}

process spring_decompress {
    tag "$meta.id"
    label 'medium'
    conda '/data/shared/programmer/miniconda3/envs/spring'

    publishDir "${outputDir}/fastqFromSpring/", mode: 'copy', pattern:"*.fastq.gz"

    input:
    tuple val(meta), path(springfile)

    output:
    tuple val(meta), path("*_R1.fastq.gz"), path("*_R2.fastq.gz"),emit: spring_fastq

    script:
    """
    spring -d \
    -i ${springfile} \
    -o ${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz \
    -g

    """
}





///////////////////////////////// PREPROCESS MODULES //////////////////////// 

process fastq_to_ubam {
    label 'medium'
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.unmapped.from.fq.bam")
    
    script:
    """
    ${gatk_exec} FastqToSam \
    -F1 ${reads[0]} \
    -F2 ${reads[1]} \
    -SM ${meta.id} \
    -PL illumina \
    -PU KGA_PU \
    -RG KGA_RG \
    --TMP_DIR ${tmpDIR} \
    -O ${meta.id}.unmapped.from.fq.bam
    """
}

process markAdapters {
    label 'medium'
    tag "$meta.id"

    input:
    tuple val(meta), path(uBAM)
    
    output:
    tuple val(meta), path("${meta.id}.ubamXT.bam"), path("${meta.id}.markAdapterMetrics.txt")
    
    script:

    """
    ${gatk_exec} MarkIlluminaAdapters \
    -I ${uBAM} \
    -O ${meta.id}.ubamXT.bam \
    --TMP_DIR ${tmpDIR} \
    -M ${meta.id}.markAdapterMetrics.txt
    """
}

process align {
    tag "$meta.id"
    label 'veryHigh'

    input:
    tuple val(meta), path(uBAM), path(metrics)

    output:
    tuple val(meta), path("${meta.id}.${genome_version}.QNsort.BWA.clean.bam")
    
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
    --TMP_DIR ${tmpDIR} \
    -O ${meta.id}.${genome_version}.QNsort.BWA.clean.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GATK: ${gatk_image}
    END_VERSIONS
    """
}

process markDup_bam {
    label 'high'
    tag "$meta.id"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    publishDir "${outputDir}/BAM/", mode: 'copy', pattern: "*.BWA.MD.ba*"
    publishDir "${outputDir}/CRAM/", mode: 'copy', pattern: "*.BWA.MD.cr*"
    publishDir "${outputDir}/Variants/Alignment_symlinks/", mode: 'link', pattern: "*.BWA.MD.cr*"


    input:
    tuple val(meta), path(aln) 
    
    output:
    tuple val(meta), path("${meta.id}${genome_version}.BWA.MD.bam"), path("${meta.id}.${genome_version}.BWA.MD*bai")
    tuple val(meta), path("${meta.id}.${genome_version}.BWA.MD.cram"), path("${meta.id}.${genome_version}.BWA.MD*crai")
    
    script:
    """
    samtools view -h ${aln[0]} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=${tmpDIR} -o ${meta.id}.${genome_version}.BWA.MD.bam /dev/stdin
    sambamba index ${meta.id}.${genome_version}.BWA.MD.bam
    
    samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${meta.id}.${genome_version}.BWA.MD.cram ${meta.id}.${genome_version}.BWA.MD.bam

    samtools index ${meta.id}.${genome_version}.BWA.MD.cram
    """
}

process markDup_cram {
    label 'high'
    tag "$meta.id"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    publishDir "${outputDir}/CRAM/", mode: 'copy', pattern: "*.BWA.MD.cr*"

    input:
    tuple val(meta), path(aln)
    
    output:
    tuple val(meta),  path("${meta.id}.${genome_version}.BWA.MD.cram"), path("${meta.id}.${genome_version}.BWA.MD*crai"), emit: markDup_output
    path  "versions.yml"                       , emit: versions

    script:
    """
    samtools view -h ${aln[0]} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=${tmpDIR} -o /dev/stdout /dev/stdin \
    |  samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${meta.id}.${genome_version}.BWA.MD.cram -

    samtools index ${meta.id}.${genome_version}.BWA.MD.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GATK: ${gatk_image}
    END_VERSIONS
    """
}

// QC PROCESSES

process bamtools {
    label 'low'
    tag "$meta.id"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    publishDir "${outputDir}/QC/", mode: 'copy'
    
    input:
    tuple val(meta),  path(aln)

    output:
    path("${meta.id}.bamtools.sample.stats.txt"), emit: multiqc

    script:
    """
    bamtools stats -in ${aln[0]} -insert > ${meta.id}.bamtools.sample.stats.txt
  
    """
  /*cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bamtools: bamtools --version| grep "bamtools" | sed 's/bamtools //g'
    END_VERSIONS
    */
}


process samtools {
    label 'low'
    tag "$meta.id"

    publishDir "${outputDir}/QC/samtools/", mode: 'copy'

    input:  
    tuple val(meta), path(aln)

    output:
    path("${meta.id}.samtools.sample.stats.txt"), emit: multiqc

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/samtools.sif samtools \
    stats \
    ${aln[0]} > ${meta.id}.samtools.sample.stats.txt
    """
}


process qualimap {
    label 'high'
    tag "$meta.id"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/qualimapSamtools/' 

    publishDir "${outputDir}/QC/qualimap/", mode: 'copy'

    input:
    tuple val(meta), path(aln)

    output:
    path ("${meta.id}/"), emit: multiqc
    path ("versions.yml"), emit: versions

    script:
    use_bed = qualimap_ROI ? "-gff ${qualimap_ROI}" : ''
    """
    qualimap --java-mem-size=40G bamqc\
    -nt ${task.cpus} \
    -outdir ${meta.id} \
    -bam ${aln[0]} $use_bed -sd -sdmode 0

    """
}

process collectWGSmetrics {
    label 'medium'
    tag "$meta.id"
    
    publishDir "${outputDir}/QC/picard/", mode: 'copy'

    input:
    tuple val(meta), path(aln)
    
    output:
    path("${meta.id}.picardWGSmetrics.txt"), emit: multiqc

    script:
    """
    ${gatk_exec} CollectWgsMetrics \
    -I ${aln[0]} \
    -O ${meta.id}.picardWGSmetrics.txt \
    -R ${genome_fasta}
    """
}

process multiQC {
    label 'low'

    publishDir "${outputDir}/QC/", mode: 'copy'

    input:
    path(inputfiles)

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
        tag "$meta.id"

        publishDir "${outputDir}/Variants/per_sample/", mode: 'copy', pattern: "*.HC.*"
        publishDir "${outputDir}/Variants/GVCF_files/", mode: 'copy', pattern: "*.g.*"
        publishDir "${outputDir}/HaplotypeCallerBAMout/", mode: 'copy', pattern: "*.HCbamout.*"
        
        if (!params.panel=="CV5") {
            publishDir "${variantStorage}/gVCF/${panelID_storage}/", mode: 'copy', pattern:'*.g.vc*' //
        }


        input:
        tuple val(meta), path(aln)
    
        output:
        path("${meta.id}.${genome_version}.g.vcf.gz"), emit: sample_gvcf

        tuple val(meta), path("${meta.id}.${genome_version}.g.vcf.gz"), emit: HC_sid_gvcf
    
        tuple val(meta), path("${meta.id}.${genome_version}.HC.*")

        path("${meta.id}.${genome_version}.g.*")
        path("${meta.id}.${genome_version}.HCbamout.*")
        tuple path("${aln[0]}"),path("${aln[1]}")
        path("versions.yml"), emit: versions

        script:
        """
        ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller \
        -I ${aln[0]} \
        -R ${genome_fasta} \
        -ERC GVCF \
        -L ${meta.roi} \
        --smith-waterman FASTEST_AVAILABLE \
        --native-pair-hmm-threads 4 \
        -pairHMM FASTEST_AVAILABLE \
        --dont-use-soft-clipped-bases \
        -O ${meta.id}.${genome_version}.g.vcf.gz \
        -bamout ${meta.id}.${genome_version}.HCbamout.bam
    
        ${gatk_exec} GenotypeGVCFs \
        -R ${genome_fasta} \
        -V ${meta.id}.${genome_version}.g.vcf.gz \
        -O ${meta.id}.${genome_version}.HC.vcf.gz \
        -G StandardAnnotation \
        -G AS_StandardAnnotation

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        GATK: ${gatk_image}
        END_VERSIONS
        """
}


process jointgenotyping {
        errorStrategy 'ignore'
        cpus 4
        //publishDir "${outputDir}/Variants/", mode: 'copy', pattern: "*.VarSeq.*"
        //publishDir "${outputDir}/Variants/GVCF_files/", mode: 'copy', pattern: "*.merged.g.*"
       
        publishDir "${outputDir}/Variants/", mode: 'copy', pattern: "*.VarSeq.*"
        publishDir "${outputDir}/Variants/GVCF_files/", mode: 'copy', pattern: "*.merged.g.*"
        //publishDir "tumorBoard_files", mode: 'copy', pattern: "*.VarSeq.*"

        input:
        tuple val(panelID), val(subpanel_gvcf) 
        //tuple val(meta),  path(vcf),path(idx) from joint_geno_dummy_ch
        output:

        path("*.for.VarSeq.*")

        script:
        """
        ${gatk_exec} --java-options "-Xmx64g" CombineGVCFs \
        -R ${genome_fasta} \
        ${subpanel_gvcf} \
        -O ${params.rundir}.${panelID}.${genome_version}.merged.g.vcf \
        -L ${ROI} \
        -G StandardAnnotation -G AS_StandardAnnotation 

        ${gatk_exec} GenotypeGVCFs \
        -R ${genome_fasta} \
        -V ${params.rundir}.${panelID}.${genome_version}.merged.g.vcf \
        -O ${params.rundir}.${panelID}.${genome_version}.merged.for.VarSeq.vcf.gz  \
        -L ${ROI} \
        -G StandardAnnotation -G AS_StandardAnnotation -A SampleList \
        -D ${dbsnp}
        """     
}

//////// WGS VARIANT CALLING (HaplotypeCaller SplitIntervals)

process haplotypecallerSplitIntervals {
    label 'low'

    input:
    tuple val(meta), path(aln), val(sub_intID), path(sub_interval) //from HC_scatter_input_bam.combine(interval_list1)

    output:
    tuple val(meta), path("${meta.id}.${sub_intID}.g.vcf"), path("${meta.id}.${sub_intID}.g.vcf.idx"), emit: hc_split_output

    script:
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller \
    -I ${aln[0]} \
    -R ${genome_fasta} \
    -ERC GVCF \
    -L ${sub_interval} \
    --smith-waterman FASTEST_AVAILABLE \
    --native-pair-hmm-threads 4 \
    -pairHMM FASTEST_AVAILABLE \
    --dont-use-soft-clipped-bases \
    -O ${meta.id}.${sub_intID}.g.vcf
    """
}


process combineGVCF {
    label 'medium'
    tag "$meta.id"

    publishDir "${variantStorage}/gVCF/${panelID_storage}/", mode: 'copy', pattern:'*.g.*' // storageDir= /lnx01_data3/storage/alignedData/hg38/


    input:

    tuple val(meta), path(sub_gvcf), path(sub_gvcf_idx)// from hc_split_output.groupTuple()
    
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.g.vcf.gz"), path("${meta.id}.*.tbi"), emit: singleGVCF
    path("${meta.id}.${genome_version}.g.vcf.gz"), emit: sample_gvcf_list_scatter
    script:
    """
    ${gatk_exec} CombineGVCFs \
    -R ${genome_fasta} \
    ${sub_gvcf.collect { "-V $it " }.join()} \
    -O ${meta.id}.${genome_version}.g.vcf.gz

        """
}

process genotypeSingle {
    label 'medium'
    tag "$meta.id"

    publishDir "${outputDir}/Variants/", mode: 'copy'
    
    input:
    tuple val(meta), path(gvcf),path(index)
    output:
    path("${meta.id}.*")
    script:
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=30" GenotypeGVCFs \
    -R ${genome_fasta} \
    -V ${gvcf} \
    -O ${meta.id}.${genome_version}.fullWGS.vcf.gz 

    ${gatk_exec} SelectVariants \
    -R ${genome_fasta} \
    -V ${meta.id}.${genome_version}.fullWGS.vcf.gz \
    -L ${ROI} \
    -O ${meta.id}.${genome_version}.WES_ROI.vcf.gz

    ${gatk_exec} IndexFeatureFile \
    -I ${meta.id}.${genome_version}.WES_ROI.vcf.gz
    """
}

process jointgenoScatter{
    label 'medium'
    publishDir "${outputDir}/Variants/", mode: 'copy'

    input:
    val x //from gvcfsamples_for_GATK_scatter

    output:
    path("*.merged.fullWGS.*")// into merged_RAW_vcf_scatter
    path("*.merged.WES_ROI.*")
    
    script:

    """
    ${gatk_exec} CombineGVCFs \
    -R ${genome_fasta} ${x} \
    -O ${params.rundir}.${genome_version}.merged.g.vcf.gz \
    -G StandardAnnotation -G AS_StandardAnnotation 

    ${gatk_exec} GenotypeGVCFs \
    -R ${genome_fasta} \
    -V ${params.rundir}.${genome_version}.merged.g.vcf.gz \
    -O ${params.rundir}.${genome_version}.merged.fullWGS.vcf.gz  \
    -G StandardAnnotation -G AS_StandardAnnotation -A SampleList 
    
    ${gatk_exec} SelectVariants \
    -R ${genome_fasta} \
    -V ${params.rundir}.${genome_version}.merged.fullWGS.vcf.gz \
    -L ${ROI} \
    -O ${params.rundir}.${genome_version}.merged.WES_ROI.vcf.gz

    """     
}
//-D ${dbsnp}


/////////////////////////////// SV CALLING MODULES //////////////////////

process manta {
    label 'high'
    tag "$meta.id"

    publishDir "${inhouse_SV}/manta/raw_calls/", mode: 'copy', pattern: "${meta.id}.${genome_version}.manta.diploidSV.vcf.gz"
    publishDir "${outputDir}/structuralVariants/manta/allOutput/", mode: 'copy'
    publishDir "${outputDir}/structuralVariants/manta/", mode: 'copy', pattern: "*.{AFanno,filtered}.*"

    input:
    tuple val(meta), path(aln)

    output:
    path("${meta.id}.${genome_version}.manta.*.{vcf,vcf.gz,gz.tbi}")
  //  tuple val("${meta.id}"), path("${meta.id}.${genome_version}.manta.AFanno.frq_below5pct.vcf"), emit: mantaForSVDB
    tuple val(meta), path("${meta.id}.${genome_version}.manta.AFanno.frq_below5pct.vcf"), emit: mantaForSVDB
    
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/manta1.6_strelka2.9.10.sif configManta.py \
    --bam ${aln[0]} \
    --referenceFasta ${genome_fasta} \
    --callRegions ${manta_callable_regions} \
    --runDir manta

    singularity run -B ${s_bind} ${simgpath}/manta1.6_strelka2.9.10.sif ./manta/runWorkflow.py -j ${task.cpus}

    mv manta/results/variants/candidateSmallIndels.vcf.gz \
    ${meta.id}.${genome_version}.manta.candidateSmallIndels.vcf.gz
    
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    ${meta.id}.${genome_version}.manta.candidateSmallIndels.vcf.gz.tbi
    
    mv manta/results/variants/candidateSV.vcf.gz \
    ${meta.id}.${genome_version}.manta.candidateSV.vcf.gz
    
    mv manta/results/variants/candidateSV.vcf.gz.tbi \
    ${meta.id}.${genome_version}.manta.candidateSV.vcf.gz.tbi

    mv manta/results/variants/diploidSV.vcf.gz \
    ${meta.id}.${genome_version}.manta.diploidSV.vcf.gz
    
    mv manta/results/variants/diploidSV.vcf.gz.tbi \
    ${meta.id}.${genome_version}.manta.diploidSV.vcf.gz.tbi

    gzip -dc ${meta.id}.${genome_version}.manta.diploidSV.vcf.gz > ${meta.id}.${genome_version}.manta.diploidSV.vcf

    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --query \
    --query_vcf ${meta.id}.${genome_version}.manta.diploidSV.vcf \
    --sqdb ${mantaSVDB} > ${meta.id}.${genome_version}.manta.AFanno.vcf 

    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${meta.id}.${genome_version}.manta.AFanno.vcf \
    --exclude-filtered \
    -select "FRQ>0.05" \
    -invert-select \
    -O ${meta.id}.${genome_version}.manta.AFanno.frq_below5pct.vcf

    """
}

process lumpy {
    label 'high'
    tag "$meta.id"

    publishDir "${inhouse_SV}/lumpy/raw_calls/", mode: 'copy', pattern: "${meta.id}.${genome_version}.Lumpy.all.vcf.gz"
    publishDir "${outputDir}/structuralVariants/lumpy/", mode: 'copy'

    input:
    tuple val(meta), path(aln)

    output:
    path("${meta.id}.${genome_version}.Lumpy.all.vcf.*") 
    tuple val(meta), path("${meta.id}.${genome_version}.lumpy.AFanno.frq_below5pct.vcf"), emit: lumpyForSVDB

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/smoove.sif smoove call -d \
    --outdir ${params.rundir}.LumpyAltSingle \
    --exclude ${smoove_exclude} \
    --name ${meta.id}.${genome_version} \
    --fasta ${genome_fasta} \
    -p ${task.cpus} \
    --genotype ${aln[0]}
    
    mv ${params.rundir}.LumpyAltSingle/${meta.id}.${genome_version}*.genotyped.vcf.gz \
    ${meta.id}.${genome_version}.Lumpy.raw.vcf.gz

    bcftools annotate -x \
    INFO/PREND,INFO/PRPOS \
    ${meta.id}.${genome_version}.Lumpy.raw.vcf.gz \
    | bgzip > ${meta.id}.${genome_version}.Lumpy.all.vcf.gz

    tabix -p vcf ${meta.id}.${genome_version}.Lumpy.all.vcf.gz

    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --query \
    --query_vcf ${meta.id}.${genome_version}.Lumpy.all.vcf.gz \
    --sqdb ${lumpySVDB} > ${meta.id}.${genome_version}.lumpy.AFanno.vcf 

    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${meta.id}.${genome_version}.lumpy.AFanno.vcf  \
    -select "FRQ>0.05" \
    -invert-select \
    -O ${meta.id}.${genome_version}.lumpy.AFanno.frq_below5pct.vcf

    """
}



process delly126 {
    tag "$meta.id"
    label 'medium'

    publishDir "${inhouse_SV}/delly/raw_calls/", mode: 'copy', pattern: "*.raw.*"
    publishDir "${outputDir}/structuralVariants/delly/", mode: 'copy'

    input:
    tuple val(meta), path(aln)

    output:
    tuple val(meta), path("${meta.id}.${genome_version}.delly.raw.vcf")
    tuple val(meta), path("${meta.id}.${genome_version}.delly.AFanno.frq_below5pct.vcf"), emit: dellyForSVDB

    script:
    """
    /data/shared/programmer/BIN/delly126 call \
    -g ${genome_fasta} \
    ${aln[0]} > ${meta.id}.${genome_version}.delly.raw.vcf

    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --query \
    --query_vcf ${meta.id}.${genome_version}.delly.raw.vcf \
    --sqdb ${dellySVDB} > ${meta.id}.${genome_version}.delly.AFanno.vcf 

    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${meta.id}.${genome_version}.delly.AFanno.vcf  \
    --exclude-filtered \
    -select "FRQ>0.05" \
    -invert-select \
    -O ${meta.id}.${genome_version}.delly.AFanno.frq_below5pct.vcf
    """
}

process cnvkit {
    label 'medium'
    tag "$meta.id"

    publishDir "${outputDir}/structuralVariants/cnvkit/", mode: 'copy'
    publishDir "${inhouse_SV}/CNVkit/CNNfiles/", mode: 'copy', pattern: '*.cnn'

    input:
    tuple val(meta), path(aln)

    output:
    path("${meta.id}.${genome_version}.cnvkit/*")
    path("*.targetcoverage.cnn"), emit: cnvkit_cnn_out
    tuple val(meta), path("${meta.id}.${genome_version}.cnvkit/*.call.cns"),path("${meta.id}.${genome_version}.cnvkit/*.cnr"), emit: CNVcalls

    script:
    """
    mv ${aln[1]} intermediate_crai
    cp intermediate_crai ${aln[1]}
    rm intermediate_crai
    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py batch \
    ${aln[0]} \
    -m wgs \
    -p ${task.cpus} \
    -r ${cnvkit_germline_reference_PON} \
    --scatter --diagram \
    -d ${meta.id}.${genome_version}.cnvkit/
    cp ${meta.id}.${genome_version}.cnvkit/*.targetcoverage.cnn .
    """
}

process cnvkitExportFiles {
    label 'medium'
    tag "$meta.id"

    publishDir "${inhouse_SV}/CNVkit/raw_calls/", mode: 'copy', pattern: '*.cnvkit.vcf'
    publishDir "${outputDir}/structuralVariants/cnvkit/", mode: 'copy'

    input:
    tuple val(meta), path(cnvkit_calls),path(cnvkit_cnr)

    output:
    path("*.vcf")
    path("*.seg")
    tuple val(meta), path("${meta.id}.${genome_version}.cnvkit.AFanno.frq_below5pct.vcf"), emit: cnvkitForSVDB

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py export vcf \
    ${cnvkit_calls} \
    -i ${meta.id} \
    -o ${meta.id}.${genome_version}.cnvkit.vcf

    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py export seg \
    ${cnvkit_cnr} \
    -o ${meta.id}.${genome_version}.cnvkit.cnr.seg

    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --query \
    --query_vcf ${meta.id}.${genome_version}.cnvkit.vcf \
    --sqdb ${cnvkitSVDB} > ${meta.id}.${genome_version}.cnvkit.AFanno.vcf 

    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${meta.id}.${genome_version}.cnvkit.AFanno.vcf  \
    -select "FRQ>0.05" \
    -invert-select \
    -O ${meta.id}.${genome_version}.cnvkit.AFanno.frq_below5pct.vcf

    """
}

process merge4callerSVDB {
    tag "$meta.id"
    label 'low'

    publishDir "${outputDir}/structuralVariants/SVDB_merged/60pctOverlap/", mode: 'copy', pattern: "*.60pctOverlap.*"
    publishDir "${outputDir}/structuralVariants/SVDB_merged/80pctOverlap/", mode: 'copy', pattern: "*.80pctOverlap.*"
    publishDir "${outputDir}/structuralVariants/SVDB_merged/100pctOverlap/", mode: 'copy', pattern: "*.100pctOverlap.*"

    input:
    tuple val(meta), path(manta_vcf), path(lumpy_vcf),path(cnvkit_vcf),path(delly_vcf)

    output:
    path("${meta.id}.4callerNEW.SVDB.*")
    path("${meta.id}.*.SVDB.*")

    script:
    """
    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --merge \
    --overlap 0.6 \
    --vcf ${manta_vcf}:MANTA ${lumpy_vcf}:LUMPY ${cnvkit_vcf}:CNVKIT ${delly_vcf}:DELLY \
    --priority LUMPY,MANTA,CNVKIT,DELLY > ${meta.id}.4callerNEW.SVDB.5pctAF.60pctOverlap.vcf

    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --merge \
    --overlap 0.8 \
    --vcf ${manta_vcf}:MANTA ${lumpy_vcf}:LUMPY ${cnvkit_vcf}:CNVKIT ${delly_vcf}:DELLY \
    --priority LUMPY,MANTA,CNVKIT,DELLY > ${meta.id}.4callerNEW.SVDB.5pctAF.80pctOverlap.vcf


    singularity exec  \
    --bind ${s_bind} /data/shared/programmer/FindSV/FindSV.simg svdb \
    --merge \
    --overlap 1.0 \
    --vcf ${manta_vcf}:MANTA ${lumpy_vcf}:LUMPY ${cnvkit_vcf}:CNVKIT ${delly_vcf}:DELLY \
    --priority LUMPY,MANTA,CNVKIT,DELLY > ${meta.id}.4callerNEW.SVDB.5pctAF.100pctOverlap.vcf
    """
}

process expansionHunter {
    label 'low'
    tag "$meta.id"

    publishDir "${outputDir}/repeatExpansions/expansionHunter/", mode: 'copy'
    
    input:
    tuple val(meta), path(aln)

    output:
    path("*.{vcf,json,bam}")
    
    script:
    """
    /data/shared/programmer/BIN/ExpansionHunter500 \
    --reads ${aln[0]} \
    --reference ${genome_fasta} \
    --variant-catalog ${expansionhunter_catalog} \
    -n ${task.cpus} \
    --output-prefix ${meta.id}.expansionHunter
    """
}

process stripy {
    label 'low'
    tag "$meta.id"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/py38/' // contains python modules required by stripy

    publishDir "${outputDir}/repeatExpansions/STRipy_ALL/", mode: 'copy',pattern:"*.ALL.html"
    publishDir "${outputDir}/repeatExpansions/STRipy_ataksi/", mode: 'copy',pattern:"*.ataksi.html"
    publishDir "${outputDir}/repeatExpansions/STRipy_myotoni/", mode: 'copy',pattern:"*.myotoni.html"
    publishDir "${outputDir}/repeatExpansions/STRipy_neuropati/", mode: 'copy',pattern:"*.neuropati.html"
    publishDir "${outputDir}/repeatExpansions/STRipy_ALS_FTD/", mode: 'copy',pattern:"*.ALS_FTD.html"
    publishDir "${outputDir}/repeatExpansions/STRipy_myopati/", mode: 'copy',pattern:"*.myopati.html"
    publishDir "${outputDir}/repeatExpansions/STRipy_epilepsi/", mode: 'copy',pattern:"*.epilepsi.html"

    input:
    tuple val(meta), path(aln)

    output:
    path("*.html")

    script:
    """
    mkdir ${meta.id}.stripy/
    sleep 5

    python3 /data/shared/programmer/stripy-pipeline-main/stri.py \
    --genome ${params.genome} \
    --reference ${genome_fasta} \
    --locus ABCD3,AFF2,AR,ARX_1,ARX_2,ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,ATXN8OS,BEAN1,C9ORF72,CACNA1A,CBL,CNBP,COMP,CSTB,DAB1,DIP2B,DMD,DMPK,EIF4A3,FGF14,FMR1,FOXL2,FXN,GIPC1,GLS,HOXA13_1,HOXA13_2,HOXA13_3,HOXD13,HTT,JPH3,LRP12,MARCHF6,NIPA1,NOP56,NOTCH2NLC,NUTM2B-AS1,PABPN1,PHOX2B,PPP2R2B,PRDM12,PPNP,RAPGEF2,RFC1,RILPL1,RUNX2,SAMD12,SOX3,STARD7,TBP,TBX1,TCF4,THAP11,TNRC6A,VWA1,XYLT1,YEATS2,ZFHX3,ZIC2,ZIC3 \
    --output ${meta.id}.stripy/ \
    --input ${aln[0]}

    mv ${meta.id}.stripy/${aln[0]}.html ${meta.id}.stripy.ALL.html

    python3 /data/shared/programmer/stripy-pipeline-main/stri.py \
    --genome ${params.genome} \
    --reference ${genome_fasta} \
    --locus ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,ATXN8OS,BEAN1,CACNA1A,CSTB,DAB1,FGF14,FMR1,FXN,NOP56,NOTCH2NLC,PPP2R2B,RFC1,TBP,YEATS2 \
    --output ${meta.id}.stripy/ \
    --input ${aln[0]}

    mv ${meta.id}.stripy/${aln[0]}.html ${meta.id}.stripy.ataksi.html

    python3 /data/shared/programmer/stripy-pipeline-main/stri.py \
    --genome ${params.genome} \
    --reference ${genome_fasta} \
    --locus ATN1,ATXN1,ATXN2,ATXN3,ATXN10,ATXN80S,C9ORF72,CACNA1A,FXN,JPH3,NOTCH2NLC,PPP2R2B,TBP \
    --output ${meta.id}.stripy/ \
    --input ${aln[0]}

    mv ${meta.id}.stripy/${aln[0]}.html ${meta.id}.stripy.myotoni.html

    python3 /data/shared/programmer/stripy-pipeline-main/stri.py \
    --genome ${params.genome} \
    --reference ${genome_fasta} \
    --locus RFC1 \
    --output ${meta.id}.stripy/ \
    --input ${aln[0]}

    mv ${meta.id}.stripy/${aln[0]}.html ${meta.id}.stripy.neuropati.html

    python3 /data/shared/programmer/stripy-pipeline-main/stri.py \
    --genome ${params.genome} \
    --reference ${genome_fasta} \
    --locus AR,ATXN2,C9ORF72,NOP56,NOTCH2NLC \
    --output ${meta.id}.stripy/ \
    --input ${aln[0]}
    mv ${meta.id}.stripy/${aln[0]}.html ${meta.id}.stripy.ALS_FTD.html

    python3 /data/shared/programmer/stripy-pipeline-main/stri.py \
    --genome ${params.genome} \
    --reference ${genome_fasta} \
    --locus CNBP,DMD,DMPK,GIPC1,LRP12,NOTCH2NLC,NUTM2B-AS1,PABPN1,RILPL1 \
    --output ${meta.id}.stripy/ \
    --input ${aln[0]}
    mv ${meta.id}.stripy/${aln[0]}.html ${meta.id}.stripy.myopati.html

    python3 /data/shared/programmer/stripy-pipeline-main/stri.py \
    --genome ${params.genome} \
    --reference ${genome_fasta} \
    --locus CSTB,MARCHF6,RAPGEF2,SAMD12,STARD7,TNRC6A,YEATS2 \
    --output ${meta.id}.stripy/ \
    --input ${aln[0]}
    mv ${meta.id}.stripy/${aln[0]}.html ${meta.id}.stripy.epilepsi.html

    """
}

process prepareManifestSMN {
    label 'low'  

    input:
    path(samplesheet) 
    
    output:
    path("SMNmanifest.txt")//, emit: SMN_manifest_ch

    shell:
    '''
    cat !{samplesheet} | cut -f2 > SMNmanifest.txt
    '''
}

process smnCopyNumberCaller {
    label 'medium'
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/py38/' // contains python modules required by smncopynumbercaller

    publishDir "${outputDir}/SMNcallerWGS/", mode: 'copy'

    input:
    path(manifest)

    output:
    path("*.{tsv,pdf,json}")
    
    script:
    """    
    python /data/shared/programmer/SMNCopyNumberCaller-1.1.2/smn_caller.py \
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

process vntyper_newRef {
    label 'high'

    publishDir "${outputDir}/MUC1-VNTR_kestrel/", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*/*.{tsv,vcf}")

    script:
    def reads_command = "--fastq1 ${reads[0]} --fastq2 ${reads[1]}"
  
    """
    singularity run -B ${s_bind} ${simgpath}/vntyper120.sif \
    -ref ${vntyperREF}/chr1.fa \
    --fastq1 ${r1} --fastq2 ${r2} \
    -t ${task.cpus} \
    -w vntyper \
    -m ${vntyperREF}/hg19_genic_VNTRs.db \
    -o ${meta.id} \
    -ref_VNTR ${vntyperREF}/MUC1-VNTR_NEW.fa \
    --fastq \
    --ignore_advntr \
    -p /data/shared/programmer/vntyper/VNtyper/
    """
}



/////////////////////////////////////////////////////////////
/// SUBWORKFLOWS meta r1 r2 input channel///////
/////////////////////////////////////////////////////////////

workflow SUB_SPRING_DECOMPRESS {

    take:
    spring_input_ch

    main:
    inputFiles_symlinks_spring(spring_input_ch)
    spring_decompress(spring_input_ch)
    emit:
    fq_read_input_spring=spring_decompress.out.spring_fastq

}


workflow SUB_PREPROCESS {

    take:
    readsInputFinal

    main:
    inputFiles_symlinks_fq(readsInputFinal)
    fastq_to_ubam(readsInputFinal)
    markAdapters(fastq_to_ubam.out[0])
    align(markAdapters.out)
    markDup_cram(align.out)
    emit:
    finalAln=markDup_cram.out.markDup_output
}
workflow SUB_QC {
    take:
    meta_aln_index
    main:
    collectWGSmetrics(meta_aln_index)
    samtools(meta_aln_index)
    //multiQC(samtools.out.multiqc.ifEmpty([]).mix(qualimap.out.multiqc.ifEmpty([])).mix(collectWGSmetrics.out.multiqc.ifEmpty([])))

}


/////////////////////////////////////////////////////////////
/// SUBWORKFLOWS meta-aln-index input channel///////
/////////////////////////////////////////////////////////////

workflow SUB_VARIANTCALL {
    take:
    meta_aln_index  // alnInputFinal
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

    if (panelID=="MV1"){
        gvcf_list
                .filter {it =~/MV1/}
                .map{" -V "+ it[1] }
                .collectFile(name: "collectfileTEST_MV1.txt", newLine: false)
                .map {it.text.trim()}
                .map { tuple("MV1",it) }
                .set {gvcfsamples_for_GATK}    
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

    if (panelID != "AV1" && panelID!= "WES_subpanel" && panelID!= "MV1") {
        gvcf_list
            .map{" -V "+ it[1] }
            .collectFile(name: "collectfileNOTAV1.txt", newLine: false)
            .map {it.text.trim()}
            .map { tuple(panelID, it) }
            .set {gvcfsamples_for_GATK}
    }

    if (!params.skipJointGenotyping) {
        jointgenotyping(gvcfsamples_for_GATK)
    }

}


workflow SUB_VARIANTCALL_WGS {
    take:
    meta_aln_index
    main:
    haplotypecallerSplitIntervals(meta_aln_index.combine(haplotypecallerIntervalList))

    if (params.genome=="hg38") {
    combineGVCF(haplotypecallerSplitIntervals.out.groupTuple(size:17))
    }
    if (params.genome=="t2t") {
    combineGVCF(haplotypecallerSplitIntervals.out.groupTuple(size:15))
    }

    genotypeSingle(combineGVCF.out.singleGVCF)

    //    if (!params.single) {
    combineGVCF.out.sample_gvcf_list_scatter
    .map{" -V "+ it }
    .set{gvcflist_scatter_done}
    
    gvcflist_scatter_done
    .collectFile(name: "collectfileTEST_scatter.txt", newLine: false)
    .map {it.text.trim()}.set {gvcfsamples_for_GATK_scatter}

    if (!params.skipJointGenotyping) {
        jointgenoScatter(gvcfsamples_for_GATK_scatter)
    }
}
workflow SUB_CNV_SV {
    take:
    meta_aln_index
    main:
    manta(meta_aln_index)
   // filter_manta(manta.out.manta)   // mantafiltered for SVDB
    lumpy(meta_aln_index)
    delly126(meta_aln_index)
    if (params.genome=="hg38") { 
    cnvkit(meta_aln_index)
    cnvkitExportFiles(cnvkit.out.CNVcalls)
    merge4callerSVDB(manta.out.mantaForSVDB.join(lumpy.out.lumpyForSVDB).join(cnvkitExportFiles.out.cnvkitForSVDB).join(delly126.out.dellyForSVDB))
    }



    //merge4callerSVDB(filter_manta.out.mantaForSVDB.join(lumpy.out.lumpyForSVDB).join(cnvkitExportFiles.out.cnvkitForSVDB).join(tiddit361.out.tidditForSVDB))
    /*

    manta.out.mantaForSVDB.join(lumpy.out.lumpyForSVDB).join(cnvkitExportFiles.out.cnvkitForSVDB).join(delly126.out.dellyForSVDB)
    |set {4callerMerge_ch}
*/

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
    smn_input_ch 
    main:

    prepareManifestSMN(smn_input_ch)
    smnCopyNumberCaller(prepareManifestSMN.out)
}
