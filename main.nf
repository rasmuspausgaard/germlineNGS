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


def helpMessage() {
    log.info"""

    Usage:

    The only requirement is a samplesheet containing 4 columns without headerline in this specific order:
    famID/projektNavn, NPN, Relation, SampleStatus

    Example samplesheet for standard trio:
    johnDoe 123456789012    index   affected
    johnDoe 234567890123    mater   normal
    johnDoe 345678901234    pater   normal

    The above information can usually be extracted directly from the sample overview excel file

    If the inputdata (FastQ or CRAM) have been transferred to the data archive (which it is by default), the script will automatically find the relevant inputdata  and create symlinks for them in the output (results) directory.

    The script will automatically look for FastQ or CRAM files in subfolders at /lnx01_data2/shared/dataArchive/. This location contains read-only access to the data archive, containing all FastQ and CRAM files. There's no need to copy or move any data.

    The user can point to a specific folder containing raw data (FastQ) using the --fastq option  or alignment data (CRAM) using the --cram option
    This is only needed if input data (FastQ or CRAM) exists outside the data archive (e.g. if data are in personal folders), or if the script is run without samplesheet.

    If the script is run without samplesheet, the user MUST point to a folder containing inputdata with either the --fastq or --cram option.

    COMMON USAGE:
    If CRAM files exist:
    nextflow run KGVejle/wgsPipeline -r main --samplesheet /path/to/samplesheet/

    If CRAM files do not exist, or the user wants to start from FastQ:
    nextflow run KGVejle/wgsPipeline -r main --samplesheet /path/to/samplesheet/ --fastqInput

    Main options:
      --help            Print this help message
      
      --genome          hg19 or hg38
                            Default: hg38 v3 (NGC version incl. masked + decoys)

      --hg38v1          Use primary (full) hg38 assembly (UCSC primary).

      --hg38v2          Use hg38 v2 (ucsc.hg38.NGS.analysisSet.fa).
      
      --samplesheet     Path to samplesheet for samples to be analyzed  
      
      --fastq            Path to folder with wgs fastq files
                            Default: /lnx01_data2/shared/dataArchive/{all subfolders}

      --fastqInput      Use fastq as input (i.e. perform trimming, and alignment)
                            Default: Not set - use CRAM as input.
      
      --cram             Path to folder with wgs CRAM files
                            Default: /lnx01_data2/shared/dataArchive/{all subfolders}

      --single           Analyze all samples in samplesheet as single WGS (e.g no jointgenotyping, merging of SV calls etc.)
                            Default: Not set - analyze all samples together
      
      --outdir          Manually set output directory
                            Default: {current_dir}/WGS_results.{DATE}
    
    Select or modify analysis steps (optional - do not use them for regular standard analyses):
      --skipVariants    Do not call SNPs and INDELs at all
                            Default: Call SNPs and INDELs using GATK HaplotypeCaller

      --skipSV          Do not call Structural Variants (SV incl. CNVs) at all
                            Default: Call SVs using Manta, Lumpy, CNVNator and CNVKit

      --skipSTR         Do not call repeat expansions.
                            Default: Calls repeat expansions using Stripy and ExpansionHunter

      --skipQC          Do not run QC module (e.g. Picard Metrics, samtools, multiQC etc.)
                            Default: Run QC module

      --skipSMN         Do not call SMN1 and SMN2 variants
                            Default: Call SMN variants with SMNCopyNumberCaller

    """.stripIndent()
}
if (params.help) exit 0, helpMessage()

def errorMessage1() {

    log.info"""

    USER INPUT ERROR: If no samplesheet is selected, the user needs to point to a folder containing relevant fastq or CRAM files... 
    Run the script with the --help parameter to see available options
    
    """.stripIndent()
}

if (!params.samplesheet && !params.fastq && !params.cram) exit 0, errorMessage1()

def FastqCRAM_error() {
    log.info"""
    USER INPUT ERROR: The user should point to either FastQ (--fastq parameter) or CRAM (--cram parameter) as input - not both! 
    """.stripIndent()
}
if (params.cram && params.fastq) exit 0, FastqCRAM_error()



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
         smnCopyNumberCaller } from "./modules/modules.dna.v1.nf" 


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


