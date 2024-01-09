#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"


//Unset parameters
params.help                     =false
params.panel                    =null
params.samplesheet              =null
params.preprocessOnly           =null

params.hg38v1                   =null
params.hg38v2                   =null
params.cram                     =null
params.fastq                    =null
params.archiveStorage           =null
params.lnx01_storage            =null
params.skipSpliceAI             =null
params.fastqInput               =null
params.skipSV                   =null
params.skipVariants             =null
params.skipQC                   =null
params.skipSTR                  =null
params.skipSMN                  =null
//Preset parameters:
params.gatk                     =null

params.server                   = "lnx01"
params.genome                   = "hg38"
params.outdir                   = "${launchDir.baseName}.Results"
params.rundir                   = "${launchDir.baseName}"
params.intervals_list           ="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";



def helpMessage() {
    log.info"""

    Usage:

    KG Vejle Germline script (WGS, WES or panels)

    PANEL ANALYSIS:


    WGS ANALYSIS:

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

    Main options:
      --help            Print this help message
      
      --genome          hg19 or hg38
                            Default: hg38 v3 (masked + decoys)

      --hg38v1          Use primary (full) hg38 assembly (UCSC primary).

      --hg38v2          Use hg38 v2 (ucsc.hg38.NGS.analysisSet.fa).

      --gatk            "danak" (v.4.1.9) or "new" (v.4.4.0.0)
                            Default: danak  
      
      --samplesheet     Path to samplesheet for samples to be analyzed (Only required for WGS analysis)
      
      --fastq            Path to folder with wgs fastq files
                            Default: /lnx01_data2/shared/dataArchive/{all subfolders}

      --fastqInput      Use fastq as input (i.e. perform trimming, and alignment)
                            Default: Not set - use CRAM as input.
      
      --cram             Path to folder with wgs CRAM files
                            Default: /lnx01_data2/shared/dataArchive/{all subfolders}
      
      --outdir          Manually set output directory
                            Default: {current_dir}/Results

    Panel analysis:

      --skipSpliceAI    Do not run SpliceAI (NB: spliceAI onlu for AV1 panel data currently)
                            Default: not set - run SpliceAI



    WGS Analysis: Select or modify analysis steps:

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

    case "AV1":
        reads_pattern_cram="*{.,-,_}{AV1}{.,-,_}*.cram";
        reads_pattern_crai="*{.,-,_}{AV1}{.,-,_}*.crai";
        reads_pattern_fastq="*{.,-,_}{AV1}{.,-,_}*R{1,2}*{fq,fastq}.gz";
        panelID="AV1"
    break;

    case "CV5":
        reads_pattern_cram="*{.,-,_}{CV5}{.,-,_}*.cram";
        reads_pattern_crai="*{.,-,_}{CV5}{.,-,_}*.crai";
        reads_pattern_fastq="*{.,-,_}{CV5}{.,-,_}*R{1,2}*{fq,fastq}.gz";
        panelID="CV5"
    break;

    case "GV3":
        reads_pattern_cram="*{GV1,GV2,GV3}*.cram";
        reads_pattern_crai="*{GV1,GV2,GV3}*.crai";
        reads_pattern_fastq="*{GV1,GV2,GV3}*R{1,2}*{fq,fastq}.gz";
        panelID="GV3"
    break;

    case "GV_TEST":
        reads_pattern_cram="*.cram";
        reads_pattern_crai="*.crai";
        reads_pattern_fastq="*R{1,2}*{fq,fastq}.gz";
        panelID="GV_TEST"
    break;


    case "WES_2":
        reads_pattern_cram="*{-,.,_}{EV8,EV7,EV6}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{EV8,EV7,EV6}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{EV8,EV7,EV6}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        panelID="WES"
    break;

    case "WES":
        reads_pattern_cram="*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        panelID="WES_subpanel"
    break;

    case "WGS_CNV":
        reads_pattern_cram="*{-,.,_}{WG4_CNV}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{WG4_CNV}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{WG4_CNV}{-,.,_}*R{1,2}*{fq,fastq}.gz";
        panelID="WGS"
    break;

    default: 
        reads_pattern_cram="*{-,.,_}{WG3,WG4,LIB,WG4_CNV,WGSmerged}{-,.,_}*.cram";
        reads_pattern_crai="*{-,.,_}{WG3,WG4,LIB,WG4_CNV,WGSmerged}{-,.,_}*.crai";
        reads_pattern_fastq="*{-,.,_}{WG3,WG4,LIB,WG4_CNV,WGSmerged}{-,.,_}*R{1,2}*{fq,fastq}.gz";
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

// Standard use: point to fastq folder for paneldata

if (!params.samplesheet && params.fastq) {
// If NOT samplesheet (std panel run), set sampleID == NPN_PANEL_SUBPANEL
    Channel
    .fromPath(params.reads, checkIfExists: true)
    .filter {it =~/_R1_/}
    .map { tuple(it.baseName.tokenize('-').get(0)+"_"+it.baseName.tokenize('-').get(1),it) }
    .set { sampleid_R1}

    Channel
    .fromPath(params.reads, checkIfExists: true)
    .filter {it =~/_R2_/}
    .map { tuple(it.baseName.tokenize('-').get(0)+"_"+it.baseName.tokenize('-').get(1),it) }
    .set { sampleid_R2 }

    sampleid_R1.join(sampleid_R2)
    .set { read_pairs_ch }

}

if (params.samplesheet && params.fastq || params.fastqInput) {
// If samplesheet, reduce sampleID to NPN only (no panel/subpanel info!)
    Channel
    .fromPath(params.reads, checkIfExists: true)
    .filter {it =~/_R1_/}
    .map { tuple(it.baseName.tokenize('-').get(0),it) }
    .set { sampleid_R1}

    Channel
    .fromPath(params.reads, checkIfExists: true)
    .filter {it =~/_R2_/}
    .map { tuple(it.baseName.tokenize('-').get(0),it) }
    .set { sampleid_R2 }

    sampleid_R1.join(sampleid_R2)
    .set { read_pairs_ch }

}


// Standard use: Point to fastq for WGS ana




if (params.cram && !params.panel) {

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

//if (params.cram && (params.panel || params.samplesheet)) {
    if (params.cram && params.panel) {
    cramfiles="${params.cram}/${reads_pattern_cram}"
    craifiles="${params.cram}/${reads_pattern_crai}"

    Channel
    .fromPath(cramfiles)
    .map { tuple(it.baseName.tokenize('.').get(0),it) }
    .set { sampleID_cram }

    Channel
    .fromPath(craifiles)
    .map { tuple(it.baseName.tokenize('.').get(0),it) }
    .set {sampleID_crai }
}


// If only samplesheet is provided, use CRAM from archive as input (default setup)!

if (params.samplesheet && !params.cram && !params.fastqInput && !params.fastq) {
    cramfiles="${dataArchive}/{lnx01,tank_kga_external_archive}/**/${reads_pattern_cram}"
    craifiles="${dataArchive}/{lnx01,tank_kga_external_archive}/**/${reads_pattern_crai}"

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

if (params.samplesheet && !params.cram && (params.fastqInput||params.fastq)) {
    full_samplesheet.join(read_pairs_ch)
    .map {tuple (it[0]+"_"+it[1]+"_"+it[2],it[4],it[5])}
    .set { fq_read_input }
}

if (params.samplesheet && !params.fastqInput && !params.fastq) {

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
         inputFiles_symlinks_cram;
         // Preprocess tools:
         //QC tools
         samtools;
         qualimap;
         fastqc_bam;
         collectWGSmetrics;
         multiQC;
         //subworkflows:
         SUB_PREPROCESS;
         SUB_VARIANTCALL;
         SUB_VARIANTCALL_WGS;
         SUB_CNV_SV;
         SUB_STR;
         SUB_SMN } from "./modules/modules.dna.v1.nf" 




workflow QC {
    take: 
    meta_aln_index
    main:
    samtools(meta_aln_index)
//    qualimap(meta_aln_index)
//    fastqc_bam(meta_aln_index)
    multiQC(samtools.out.ifEmpty([]).mix(qualimap.out.ifEmpty([])).mix(fastqc_bam.out.ifEmpty([])).collect())

}



workflow {

    if (!params.panel) { 

        if (params.fastqInput||params.fastq) {
            SUB_PREPROCESS(fq_read_input)
            
            if (!params.skipVariants) {
                SUB_VARIANTCALL_WGS(SUB_PREPROCESS.out.finalAln)
            }
            if (!params.skipSV) {
                SUB_CNV_SV(SUB_PREPROCESS.out.finalAln)
            }
            if (!params.skipSTR) {
                SUB_STR(SUB_PREPROCESS.out.finalAln)
            }
            
            if (!params.skipSMN) {
            SUB_SMN(SUB_PREPROCESS.out.finalAln)
            }
        }

        if (!params.fastqInput && !params.fastq) {
            inputFiles_symlinks_cram(meta_aln_index)

            if (!params.skipVariants) {
                SUB_VARIANTCALL_WGS(meta_aln_index)
            }
            if (!params.skipSV) {
                SUB_CNV_SV(meta_aln_index)
            }
            if (!params.skipSTR) {
                SUB_STR(meta_aln_index)
            }
            if (!params.skipSMN) {
            SUB_SMN(meta_aln_index)
            }
        }
    }

    if (params.panel) {

        if (params.fastqInput||params.fastq) {
            SUB_PREPROCESS(fq_read_input)
            SUB_VARIANTCALL(SUB_PREPROCESS.out.finalAln)
        }
        if (!params.fastqInput && !params.fastq) {
            inputFiles_symlinks_cram(meta_aln_index)
            SUB_VARIANTCALL(meta_aln_index)
        }
    }



}
