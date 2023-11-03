#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"


//Unset parameters
params.help                     =false
params.samplesheet              =null
params.preprocessOnly           =null

params.hg38v1                   =null
params.hg38v2                   =null
params.outputBam                =null
params.cram                     =null
params.fastq                    =null
params.archiveStorage           =null
params.lnx01_storage            =null
params.skipSpliceAI             =null
params.fastqInput               =null
//Preset parameters:
params.gatk                     ="new"
params.panel                    = "WGS" // preset for WGS analysis
params.server                   = "lnx01"
params.genome                   = "hg38"
params.outdir                   = "${launchDir.baseName}.Results"
params.rundir                   = "${launchDir.baseName}"
params.intervals_list           ="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";




switch (params.server) {
    case 'lnx01':
        modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        subworkflow_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows";
        dataArchive="/lnx01_data2/shared/dataArchive";
    break;
    case 'kga01':
        modules_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        subworkflow_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows";
        dataArchive="/data/shared/dataArchive";

    break;
}



switch (params.panel) {
    case "WGS":
        reads_pattern_cram="*{-,.,_}{WG3,WG4,LIB,WG4_CNV}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{WG3,WG4,LIB,WG4_CNV}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{WG3,WG4,LIB,WG4_CNV}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        panelID="WGS"
    break;

    case "WGS_CNV":
        reads_pattern_cram="*{-,.,_}{WG4_CNV}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{WG4_CNV}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{WG4_CNV}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        panelID="WGS"
    break;

    default: 
        reads_pattern_cram="*{-,.,_}{WG3,WG4,LIB,WG4_CNV}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{WG3,WG4,LIB,WG4_CNV}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{WG3,WG4,LIB,WG4_CNV}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        panelID="WGS"
    break;
}



////////////////////////////////////////////////////
////// INPUT DATA (fastq or CRAM) channels //////////
////////////////////////////////////////////////////
if (params.fastq) {

    params.reads="${params.fastq}/${reads_pattern_fastq}"
}

if (!params.fastq && params.fastqInput) {

    params.reads="${dataArchive}/{lnx01,kga01_novaRuns,tank_kga_external_archive}/**/${reads_pattern_fastq}"
}

// if fastq input, set reads input channels
if (!params.cram && params.fastqInput) {
    channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { it -> [it[0], file(it[1][0]),file(it[1][1])] }
    .set { read_pairs_ch }
}


if (params.cram) {

    cramfiles="${params.cram}/${reads_pattern_cram}"
    craifiles="${params.cram}/${reads_pattern_crai}"

    Channel
    .fromPath(cramfiles)
    .map { tuple(it.baseName.tokenize('_').get(0),it) }
    .set { sampleID_cram }

    Channel
    .fromPath(craifiles)
    .map { tuple(it.baseName.tokenize('_').get(0),it) }
    .set {sampleID_crai }
}

// If only samplesheet is provided, use CRAM from archive as input (default setup)!

if (params.samplesheet && !params.cram && !params.fastqInput) {
    cramfiles="${dataArchive}/{lnx01,kga01_novaRuns,tank_kga_external_archive}/${reads_pattern_cram}"
    craifiles="${dataArchive}/{lnx01,kga01_novaRuns,tank_kga_external_archive}/${reads_pattern_crai}"

    Channel
    .fromPath(cramfiles)
    .map { tuple(it.baseName.tokenize('_').get(0),it) }
    .set { sampleID_cram }

    Channel
    .fromPath(craifiles)
    .map { tuple(it.baseName.tokenize('_').get(0),it) }
    .set {sampleID_crai }
}


////////////////////////////////////////////////////
///////////// SAMPLESHEET channels /////////////////
////////////////////////////////////////////////////
if (params.samplesheet) {
    channel.fromPath(params.samplesheet)
        .splitCsv(sep:'\t')
        .map { row -> tuple(row[1], row[0],row[2],row[3])}
        .set { full_samplesheet }
    //above: NPN, caseID, relation, samplestatus

    channel.fromPath(params.samplesheet)
        .splitCsv(sep:'\t')
        .map { row -> row[0]}
        .unique()
        .collect()
        .set { caseID_ch }

    channel.fromPath(params.samplesheet)
        .splitCsv(sep:'\t')
        .map { row -> tuple(row[0],row[1])}
        .set {caseID_sampleID}
}


////////////////////////////////////////////////////
///////////// set final input channels   ///////////
////////////////////////////////////////////////////



if (!params.samplesheet && params.fastq) {
    read_pairs_ch
    .set { fq_read_input }
}

if (!params.samplesheet && params.cram) {
    sampleID_cram.join(sampleID_crai)
    .set { meta_aln_index }
}



if (params.samplesheet && !params.cram && params.fastqInput) {
    full_samplesheet.join(read_pairs_ch)
    .map {tuple (it[0]+"_"+it[1]+"_"+it[2],it[4],it[5])}
    .set { fq_read_input }
}

if (params.samplesheet && !params.fastqInput) {

    full_samplesheet.join(sampleID_cram).join(sampleID_crai)
    .map {tuple (it[0]+"_"+it[1]+"_"+it[2],it[4],it[5])}
    .set {meta_aln_index}
}
//////// END: Combine input and samplesheet //////////

///// Haplotypecaller splitintervals channel: /////
channel
    .fromPath(params.intervals_list)
    .map { it -> tuple(it.baseName,it)}
    .set { haplotypecallerIntervalList }
////////////////////////////////////////////////////


include { 
         // Symlinks:
         inputFiles_symlinks_fq;
         inputFiles_symlinks_cram;
         // Preprocess tools:

         fastq_to_ubam;
         align;
         markAdapters; 
         markDup_bam;
         markDup_cram;
         //QC tools
         samtools;
         qualimap;
         fastqc_bam;
         collectWGSmetrics;
         multiQC;
         //Variant call tools:
         haplotypecaller;
         jointgenotyping;
         haplotypecallerSplitIntervals;         // WGS only
         mergeScatteredGVCF;                    // WGS only
         jointgenoScatter;                      // WGS only
         // SV + CNV tools:
         tiddit361;
         manta;
         filter_manta;
         lumpy;
         cnvkit;
         cnvkitExportFiles;
         merge4callerSVDB;         
         // STR tools
         expansionHunter;                   
         stripy;
         // SMN caller tools
         prepareManifestSMN;
         smnCopyNumberCaller } from modules/KGVejle.DNA.modules.nf" 


workflow PREPROCESS {

    take:
    fq_read_input
    
    main:

    fastq_to_ubam(fq_read_input)
    markAdapters(fastq_to_ubam.out[0])
    align(markAdapters.out)
    markDup_cram(align.out)
    //markDup_v3_cram.out.markDup_output.view()
    emit:
    finalAln=markDup_cram.out.markDup_output
}


workflow QC {
    take: 
    meta_aln_index
    main:
    samtools(meta_aln_index)
//    qualimap(meta_aln_index)
//    fastqc_bam(meta_aln_index)
    multiQC(samtools.out.ifEmpty([]).mix(qualimap.out.ifEmpty([])).mix(fastqc_bam.out.ifEmpty([])).collect())

}

workflow VARIANTCALL {
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

    jointgenoScatter(gvcfsamples_for_GATK_scatter)
}

workflow CNV_SV {
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


workflow STR {
    take:
    meta_aln_index
    main:
    expansionHunter(meta_aln_index) 
    stripy(meta_aln_index)
  

}

workflow SMN {
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



workflow {

    if (params.preprocessOnly) {
        inputFiles_symlinks_fq(fq_read_input)
        PREPROCESS(fq_read_input)
    }

    if (params.fastqInput) {
        inputFiles_symlinks_fq(fq_read_input)
        PREPROCESS(fq_read_input)
        VARIANTCALL(PREPROCESS.out.finalAln)
        CNV_SV(PREPROCESS.out.finalAln)
        STR(PREPROCESS.out.finalAln)
        SMN(PREPROCESS.out.finalAln)
    }

    if (!params.fastqInput) {
        inputFiles_symlinks_cram(meta_aln_index)
        VARIANTCALL(meta_aln_index)
        CNV_SV(meta_aln_index)
        STR(meta_aln_index)
        SMN(meta_aln_index)
    }

}


