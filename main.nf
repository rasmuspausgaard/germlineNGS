#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Full Nextflow DSL2 script, focusing on CRAM inputs (optionally FASTQ).
 * Removed samplesheet logic. Collects sample names from CRAM, sends them in email on completion.
 */

/* -----------------------------------------------------------------
   Basic definitions and parameter handling
   ----------------------------------------------------------------- */
date = new Date().format('yyMMdd')
user = System.getenv('USER') ?: "unknownUser"
runID = "${date}.${user}"

// Unset parameters
params.help                = params.help                ?: false
params.panel               = params.panel               ?: null
params.preprocessOnly      = params.preprocessOnly      ?: null
params.keepwork            = params.keepwork            ?: null
params.nomail              = params.nomail              ?: null
params.hg38v1              = params.hg38v1              ?: null
params.hg38v2              = params.hg38v2              ?: null
params.cram                = params.cram                ?: null
params.fastq               = params.fastq               ?: null
params.archiveStorage       = params.archiveStorage      ?: null
params.lnx01_storage        = params.lnx01_storage       ?: null
params.skipSpliceAI         = params.skipSpliceAI        ?: null
params.skipJointGenotyping  = params.skipJointGenotyping ?: null
params.fastqInput           = params.fastqInput          ?: null
params.skipSV               = params.skipSV              ?: null
params.skipVariants         = params.skipVariants        ?: null
params.skipQC               = params.skipQC              ?: null
params.skipSTR              = params.skipSTR             ?: null
params.skipSMN              = params.skipSMN             ?: null
// Preset parameters:
params.gatk                 = params.gatk                ?: null
params.copyCram             = params.copyCram            ?: null
params.single               = params.single              ?: null
params.server               = params.server              ?: "lnx01"
params.genome               = params.genome              ?: "hg38"
params.outdir               = params.outdir              ?: "${launchDir.baseName}.Results"
params.rundir               = params.rundir              ?: "${launchDir.baseName}"
params.cram_date   = params.cram_date   ?: null
// Example intervals (if needed):
// params.intervals_list    = "/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/..."
params.reference   = params.reference   ?: '/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa'

if( !params.cram_date ) {
    if( params.cram ) {
        // strip trailing slash if present
        def cramDir = params.cram.endsWith('/') ? params.cram[0..-2] : params.cram
        def m = (new File(cramDir).name =~ /(\d{6})/)
        assert m, "Cannot find YYMMDD in cram folder name: ${params.cram}"
        params.cram_date = m[0][1]        // e.g. "250617"
        log.info "cram_date derived as ${params.cram_date}"
    } else {
        throw new IllegalArgumentException(
            "Need --cram_date or --cram from which to derive it")
    }
}
/* folder where VarSeq CNV expects its VCFs & sample-fields */
def variants_dir = "/lnx01_data2/shared/patients/hg38/WGS.CNV/2025/${params.cram_date}/${params.cram_date}.Results/Variants"
/* -----------------------------------------------------------------
   Usage / Help messages
   ----------------------------------------------------------------- */
def helpMessage() {
    log.info """
    Usage: nextflow run this_script.nf [options]

    This pipeline processes WGS or panel data from CRAM or FASTQ input.

    REQUIRED:
      --cram  <folder>   Path to folder containing CRAM (and CRAI) files
         OR
      --fastq <folder>   Path to folder containing FASTQ files
         (At least one of --cram or --fastq must be provided)

    OPTIONAL:
      --fastqInput       Use FASTQ as input (perform trimming + alignment)
                         If not set, pipeline assumes CRAM as input

      --panel            Type of panel data to analyze (AV1, CV5, MV1, etc.)
                         Default: not set => assumes WGS

      --skipVariants     Skip SNP/INDEL calling
      --skipSV           Skip structural variant calling
      --skipSTR          Skip repeat expansions
      --skipQC           Skip QC
      --skipSMN          Skip SMN calling

      --keepwork         Keep the Nextflow work folder
      --nomail           Do not send email on completion
      --server           lnx01 or lnx02 (affects how email is sent)
      --copyCram         If set, CRAM files will be physically copied
                         instead of symlinked

    EXAMPLE:
      nextflow run this_script.nf --cram /path/to/cram --panel AV1
    """.stripIndent()
}
if (params.help) {
    helpMessage()
    exit 0
}

// If user provides neither CRAM nor FASTQ, error
if (!params.cram && !params.fastq) {
    log.error "ERROR: Must provide --cram or --fastq. Use --help for usage."
    exit 1
}
// If user provides both CRAM and FASTQ, error
if (params.cram && params.fastq) {
    log.error "ERROR: Cannot provide both --cram and --fastq. Use --help for usage."
    exit 1
}


/* -----------------------------------------------------------------
   Determine server-based paths, e.g. dataArchive
   ----------------------------------------------------------------- */
switch (params.server) {
    case 'lnx02':
        dataArchive = "/lnx01_data2/shared/dataArchive"
        break

    case 'lnx01':
        dataArchive = "/lnx01_data2/shared/dataArchive"
        // modules_dir     = "/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules"
        // subworkflow_dir = "/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/subworkflows"
        break

    case 'kga01':
        dataArchive = "/data/shared/dataArchive"
        break

    default:
        dataArchive = "/lnx01_data2/shared/dataArchive"
        break
}


/* -----------------------------------------------------------------
   Panel logic: define patterns for CRAM / FASTQ
   ----------------------------------------------------------------- */
switch (params.panel) {

    case "AV1":
        reads_pattern_cram  = "*{.,-,_}{AV1}{.,-,_}*.cram"
        reads_pattern_crai  = "*{.,-,_}{AV1}{.,-,_}*.crai"
        reads_pattern_fastq = "*{.,-,_}{AV1}{.,-,_}*R{1,2}*{fq,fastq}.gz"
        panelID = "AV1"
        break

    case "CV5":
        reads_pattern_cram  = "*{.,-,_}{CV5}{.,-,_}*.cram"
        reads_pattern_crai  = "*{.,-,_}{CV5}{.,-,_}*.crai"
        reads_pattern_fastq = "*{.,-,_}{CV5}{.,_,-}*R{1,2}*{fq,fastq}.gz"
        panelID = "CV5"
        break

    case "MV1":
        reads_pattern_cram  = "*{MV1}*.cram"
        reads_pattern_crai  = "*{MV1}*.crai"
        reads_pattern_fastq = "*{MV1}*R{1,2}*{fq,fastq}.gz"
        panelID = "MV1"
        break

    case "WES":
        reads_pattern_cram  = "*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*.cram"
        reads_pattern_crai  = "*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*.crai"
        reads_pattern_fastq = "*{-,.,_}{EV8_ALM,EV8_ONK}{-,.,_}*R{1,2}*{fq,fastq}.gz"
        panelID = "WES_subpanel"
        break

    case "WGS_CNV":
        reads_pattern_cram  = "*{-,.,_}{WG4_CNV}{-,.,_}*.cram"
        reads_pattern_crai  = "*{-,.,_}{WG4_CNV}{-,.,_}*.crai"
        reads_pattern_fastq = "*{-,.,_}{WG4_CNV}{-,.,_}*R{1,2}*{fq,fastq}.gz"
        panelID = "WGS"
        break

    case "NGC":
        reads_pattern_cram  = "*{-,.,_}{WG4_NGC}{-,.,_}*.cram"
        reads_pattern_crai  = "*{-,.,_}{WG4_NGC}{-,.,_}*.crai"
        reads_pattern_fastq = "*{-,.,_}{WG4_NGC}{-,.,_}*R{1,2}*{fq,fastq}.gz"
        panelID = "NGC"
        break

    default:
        // Default: WGS
        reads_pattern_cram  = "*{-,.,_}{WG3,WG4,A_WG4,LIB,WG4_CNV,WGSmerged}{-,.,_}*.cram"
        reads_pattern_crai  = "*{-,.,_}{WG3,WG4,A_WG4,LIB,WG4_CNV,WGSmerged}{-,.,_}*.crai"
        reads_pattern_fastq = "*{-,.,_}{WG3,WG4,A_WG4,LIB,WG4_CNV,WGSmerged,WGS}{-,.,_}*R{1,2}*{fq,fastq}.gz"
        panelID = "WGS"
        break
}


/* -----------------------------------------------------------------
   INPUT DATA CHANNELS (CRAM or FASTQ)
   ----------------------------------------------------------------- */

/** 1) CRAM handling **/
if (params.cram) {
    // CRAM + CRAI
    cramfiles = "${params.cram}/${reads_pattern_cram}"
    craifiles = "${params.cram}/${reads_pattern_crai}"

    Channel
        .fromPath(cramfiles, checkIfExists: true)
        .map { file ->
            // sampleID from baseName
            def sampleID = file.baseName.tokenize('.').get(0)
            tuple(sampleID, file)
        }
        .set { sampleID_cram }

    Channel
        .fromPath(craifiles, checkIfExists: true)
        .map { file ->
            def sampleID = file.baseName.tokenize('.').get(0)
            tuple(sampleID, file)
        }
        .set { sampleID_crai }

    // Join CRAM + CRAI => meta_aln_index
    sampleID_cram.join(sampleID_crai)
        .set { meta_aln_index }
   def joinedFiles = sampleID_cram.join(sampleID_crai)
}


/** 2) FASTQ handling **/
if (params.fastq) {
    if (!params.fastqInput) {
        // If user gave --fastq but didn't set --fastqInput, warn or handle?
        log.warn "FASTQ provided but --fastqInput not set. The pipeline expects CRAM by default."
    }

    // If user sets fastq, define the path
    params.reads = "${params.fastq}/${reads_pattern_fastq}"

    Channel
        .fromPath(params.reads, checkIfExists: true)
        .filter { it.name =~ /R1/ }
        .map { file ->
            def tokens = file.baseName.tokenize('-')
            def sampleID = tokens[0] + "_" + tokens[1]
            tuple(sampleID, file)
        }
        .set { sampleid_R1 }

    Channel
        .fromPath(params.reads, checkIfExists: true)
        .filter { it.name =~ /R2/ }
        .map { file ->
            def tokens = file.baseName.tokenize('-')
            def sampleID = tokens[0] + "_" + tokens[1]
            tuple(sampleID, file)
        }
        .set { sampleid_R2 }

    // Combine R1 + R2
    sampleid_R1.join(sampleid_R2)
        .set { read_pairs_ch }
}

/* -----------------------------------------------------------------
   FINAL INPUT CHANNELS (tie CRAM or FASTQ to pipeline)
   ----------------------------------------------------------------- */

// If purely CRAM-based (no FASTQ):
if (params.cram && !params.fastq) {
    // meta_aln_index is the final input channel for alignment-based steps
    // e.g. .set { meta_aln_index }
}

// If purely FASTQ-based (no CRAM):
if (!params.cram && params.fastq) {
    // read_pairs_ch is the final input channel for alignment steps
}

// If you do a hybrid scenario, adapt as needed. (But your pipeline typically does one or the other.)


/* -----------------------------------------------------------------
   SUBWORKFLOWS / MODULES
   ----------------------------------------------------------------- */
// This part references your modules file, if you still use them.
include { 
    // Tools:
    inputFiles_symlinks_cram
    inputFiles_cramCopy
    samtools
    qualimap
    fastqc_bam
    collectWGSmetrics
    multiQC
    vntyper_newRef
    // Subworkflows:
    CalculateCoverage
    VarSeqCNV
    SUB_PREPROCESS
    SUB_VARIANTCALL
    SUB_VARIANTCALL_WGS
    SUB_CNV_SV
    SUB_STR
    SUB_SMN
} from "./modules/modules.dna.v1.nf"


/* -----------------------------------------------------------------
   QC Workflow example
   ----------------------------------------------------------------- */
workflow QC {
    take:
    meta_aln_index

    main:
    samtools(meta_aln_index)
    // qualimap(meta_aln_index)
    // fastqc_bam(meta_aln_index)

    multiQC(
        samtools.out.ifEmpty([])
        .mix(qualimap.out.ifEmpty([]))
        .mix(fastqc_bam.out.ifEmpty([]))
        .collect()
    )
}
/* -----------------------------------------------------------------
   storage list for coverage calculation
   ----------------------------------------------------------------- */
def coverageList = []
/* -----------------------------------------------------------------
   MAIN WORKFLOW
   ----------------------------------------------------------------- */

workflow {

    /* ──────────────────────────────────────────────────────────────
       0 ▪ Coverage calculation
       ────────────────────────────────────────────────────────────── */
    coverageResults = CalculateCoverage(meta_aln_index)
    coverageResults.subscribe { res ->
        coverageList << res
        println "Coverage for sample '${res[0]}': ${res[1]}"
    }

    /* ──────────────────────────────────────────────────────────────
       1 ▪ MAIN logic for WGS / WGS_CNV / NGC
       ────────────────────────────────────────────────────────────── */
    if ( !params.panel || params.panel in ['WGS_CNV','NGC'] ) {

        /* ───── FASTQ branch ───── */
        if ( params.fastqInput || params.fastq ) {

            SUB_PREPROCESS(read_pairs_ch)

            if ( !params.preprocessOnly ) {
                if ( !params.skipVariants )
                    SUB_VARIANTCALL_WGS( SUB_PREPROCESS.out.finalAln )
                if ( !params.skipSV )
                    SUB_CNV_SV( SUB_PREPROCESS.out.finalAln )
                if ( !params.skipSTR )
                    SUB_STR( SUB_PREPROCESS.out.finalAln )
                if ( !params.skipSMN )
                    SUB_SMN( SUB_PREPROCESS.out.finalAln )
            }
        }

        /* ───── CRAM branch ───── */
        else if ( params.cram ) {

            /* 1 ▪ Build (sid , sampleFile) – kun CRAM-filer der matcher *WG4_CNV*.cram */
            Channel
                .fromPath("${params.cram}/*WG4_CNV*.cram", checkIfExists:true)
                .map { cramFile ->
                    def sid = cramFile.baseName.split('\\.')[0]
                    def sampleFile = file("${variants_dir}/${sid}_sample_file.txt")
            
                    if ( !sampleFile.exists() ) {
                        sampleFile.text = "Sample,BAM Path\n${sid},${cramFile}\n"
                        log.info "sample_fields file written for $sid"
                    }
                    tuple(sid, sampleFile)               // (sid , sampleFile)
                }
                .set { sid_sampleFile }


            /* 2 ▪ CRAM→CRAI join and downstream calls */
            sampleID_cram.join(sampleID_crai)
                         .set { meta_aln_index_with_files }

            def vcalls = !params.skipVariants ?
                         SUB_VARIANTCALL_WGS(meta_aln_index_with_files) :
                         null

            if ( !params.copyCram ) {
                inputFiles_symlinks_cram(meta_aln_index_with_files)
                if ( !params.skipSV  ) SUB_CNV_SV(meta_aln_index_with_files)
                if ( !params.skipSTR ) SUB_STR(meta_aln_index_with_files)
                if ( !params.skipSMN ) SUB_SMN(meta_aln_index_with_files)
            } else {
                inputFiles_symlinks_cram(meta_aln_index_with_files)
                inputFiles_cramCopy(meta_aln_index_with_files)
                if ( !params.skipSV  ) SUB_CNV_SV(inputFiles_cramCopy.out)
                if ( !params.skipSTR ) SUB_STR(inputFiles_cramCopy.out)
                if ( !params.skipSMN ) SUB_SMN(inputFiles_cramCopy.out)
            }

            /* 3 ▪ Run VarSeq once per sample (WGS_CNV only) */
            if ( params.panel == 'WGS_CNV' && vcalls ) {
                vcalls.hc_vcfs
                      .join(sid_sampleFile)               // (sid , vcf , sampleFile)
                      .map { sid, vcf, sampleFile ->
                          tuple(vcf, sid, sampleFile)      // → VarSeqCNV
                      } | VarSeqCNV
            }
        }
    }

    /* ──────────────────────────────────────────────────────────────
       2 ▪ OTHER PANELS
       ────────────────────────────────────────────────────────────── */
    if ( params.panel && !(params.panel in ['WGS_CNV','NGC']) ) {

        if ( params.fastqInput || params.fastq ) {
            SUB_PREPROCESS(read_pairs_ch)
            SUB_VARIANTCALL( SUB_PREPROCESS.out.finalAln )

            if ( params.panel == 'MV1' )
                vntyper_newRef(read_pairs_ch)

        } else if ( params.cram ) {
            SUB_VARIANTCALL( meta_aln_index )
        }
    }
}



/* -----------------------------------------------------------------
   ON COMPLETE: send email with sample names, etc.
   ----------------------------------------------------------------- */
workflow.onComplete {
    def currentYear = new Date().format('yyyy')

    println "Coverage summary:\n${coverageSummary}"

    // (Optional) If you have an IP file for lnx02 emailing:
    def ipFilePath = '/lnx01_data2/shared/testdata/test_scripts/ip_file'
    def ip = ""
    if (params.server == 'lnx02') {
        if (new File(ipFilePath).exists()) {
            ip = new File(ipFilePath).text.trim()
        } else {
            println "Warning: IP file not found at ${ipFilePath}, might fail sending email from lnx02."
        }
    }

    // Build the sample names string
    def coverageSummary = coverageList
        .collect { tuple -> "${tuple[0]}: ${tuple[1].trim()}" }
        .join('\n')

    println "Coverage summary:\n${coverageSummary}"

    // Email conditions: pipeline success, duration > 5 minutes(300000), user is "mmaj" or "raspau", etc.
    if (!params.nomail && workflow.success && workflow.duration > 3) {
        if (user in ["mmaj", "raspau"]) {

            // Example: derive "sequencingRun" from the CRAM folder name
            def sequencingRun = params.cram
                ? new File(params.cram).getName().take(6)
                : params.fastq
                    ? new File(params.fastq).getName().take(6)
                    : 'Not provided'

            // Check for OBS sample if panel == AV1
            def obsSampleMessage = ""
            if (params.panel == "AV1" && params.cram) {
                def cramDir = new File(params.cram)
                def obsSamples = cramDir.listFiles().findAll { it.name.contains("OBS") }
                if (obsSamples.size() > 0) {
                    obsSampleMessage = "\nTHERE IS AN OBS SAMPLE IN THIS RUN"
                }
            }

            def workDirMessage = params.keepwork ? "WorkDir: ${workflow.workDir}" : "WorkDir: Deleted"
            def outputDir = "${launchDir}/${launchDir.baseName}.Results"

            def body = """|Pipeline execution summary
            |---------------------------
            |Pipeline completed: ${params.panel}
            |Sequencing run: ${sequencingRun}${obsSampleMessage}
            |Duration: ${workflow.duration}
            |Success: ${workflow.success}
            |${workDirMessage}
            |OutputDir: ${outputDir}
            |Exit status: ${workflow.exitStatus}
            |${obsSampleMessage}
            |
            |Samples included in the pipeline:
            |${coverageSummary}
            """.stripMargin('|')


            // Example recipients
            def recipients = 'Andreas.Braae.Holmgaard@rsyd.dk,Annabeth.Hogh.Petersen@rsyd.dk,Isabella.Almskou@rsyd.dk,Jesper.Graakjaer@rsyd.dk,Lene.Bjornkjaer@rsyd.dk,Martin.Sokol@rsyd.dk,Mads.Jorgensen@rsyd.dk,Rasmus.Hojrup.Pausgaard@rsyd.dk,Signe.Skou.Tofteng@rsyd.dk,Amalie.Schirmer.Ahlgreen.Larsen@rsyd.dk,Sara.Kaczor.Elbaek@rsyd.dk'

            // Send mail depending on server
            if (params.server == 'lnx01') {
                // Nextflow's built-in mail
                sendMail(to: recipients, subject: 'CRAM-based pipeline Update', body: body)
            }
            else if (params.server == 'lnx02') {
                // Use external command
                def emailCommand = "ssh ${ip} 'echo \"${body}\" | mail -s \"CRAM-based pipeline Update\" ${recipients}'"
                def proc = ['bash', '-c', emailCommand].execute()
                proc.waitFor()
                if (proc.exitValue() != 0) {
                    println("Error sending email from lnx02: ${proc.err.text}")
                } else {
                    println("Email successfully sent from lnx02.")
                }
            }

            // Move WGS_CNV from lnx02 to lnx01 if success
            if (params.server == 'lnx02' && params.panel == 'WGS_CNV' && workflow.success) {
                def moveWGSCNVCommand = "mv ${launchDir} /lnx01_data2/shared/patients/hg38/WGS.CNV/${currentYear}/"
                def moveWGSCNVProcess = ['bash', '-c', moveWGSCNVCommand].execute()
                moveWGSCNVProcess.waitFor()
                if (moveWGSCNVProcess.exitValue() != 0) {
                    println("Error moving WGS_CNV files: ${moveWGSCNVProcess.err.text}")
                }
            }

            // Move WGS_NGC from lnx02 to lnx01 if success
            if (params.server == 'lnx02' && params.panel == 'NGC' && workflow.success) {
                def moveWGSNGCCommand = "mv ${launchDir} /lnx01_data2/shared/patients/hg38/WGS_NGC/${currentYear}/"
                def moveWGSNGCProcess = ['bash', '-c', moveWGSNGCCommand].execute()
                moveWGSNGCProcess.waitFor()
                if (moveWGSNGCProcess.exitValue() != 0) {
                    println("Error moving WGS_NGC files: ${moveWGSNGCProcess.err.text}")
                }
            }

            // Move WES from lnx02 to lnx01 if success
            if (params.server == 'lnx02' && params.panel == 'WES' && workflow.success) {
                def moveWESCommand = "mv ${launchDir} /lnx01_data2/shared/patients/hg38/WES_ALM_ONK/${currentYear}/"
                def moveWESProcess = ['bash', '-c', moveWESCommand].execute()
                moveWESProcess.waitFor()
                if (moveWESProcess.exitValue() != 0) {
                    println("Error moving WES files: ${moveWESProcess.err.text}")
                }
            }
        }
    }

    // Delete work directory if not keeping it
    if (!params.keepwork) {
        println("Deleting work directory: ${workflow.workDir}")
        def deleteWorkDirCommand = "rm -rf ${workflow.workDir}".execute()
        deleteWorkDirCommand.waitFor()
        if (deleteWorkDirCommand.exitValue() != 0) {
            println("Error deleting work directory: ${deleteWorkDirCommand.err.text}")
        }
    }
}
