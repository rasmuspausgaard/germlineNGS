#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Full Nextflow DSL2 script, focusing on CRAM inputs (optionally FASTQ).
 * Collects sample names from CRAM, sends them in email on completion.
 */

/* -----------------------------------------------------------------
   Basic definitions and parameter handling
   ----------------------------------------------------------------- */
date  = new Date().format('yyMMdd')
user  = System.getenv('USER') ?: "unknownUser"
runID = "${date}.${user}"

params.help                = params.help                ?: false
params.panel               = params.panel               ?: null
params.preprocessOnly      = params.preprocessOnly      ?: null
params.keepwork            = params.keepwork            ?: null
params.nomail              = params.nomail              ?: null
params.cram                = params.cram                ?: null
params.fastq               = params.fastq               ?: null
params.fastqInput          = params.fastqInput          ?: null
params.copyCram            = params.copyCram            ?: null
params.server              = params.server              ?: "lnx01"
params.genome              = params.genome              ?: "hg38"
params.outdir              = params.outdir              ?: "${launchDir.baseName}.Results"
params.rundir              = params.rundir              ?: "${launchDir.baseName}"
// ... plus any other params you need

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

    OPTIONAL:
      --fastqInput
      --panel
      --skipVariants
      --skipSV
      --skipSTR
      --skipQC
      --skipSMN
      --keepwork
      --nomail
      --server (lnx01 or lnx02)
      --copyCram

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
if (params.cram) {
    def cramfiles = "${params.cram}/${reads_pattern_cram}"
    def craifiles = "${params.cram}/${reads_pattern_crai}"

    Channel
        .fromPath(cramfiles, checkIfExists: true)
        .map { file ->
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

    // Join CRAM + CRAI => meta_aln_index (4 elements)
    sampleID_cram.join(sampleID_crai)
        .set { meta_aln_index }

    // Create a 3-element channel for calcMeanDepth
    Channel
        .from(meta_aln_index)
        .map { sampleID1, cramFile, sampleID2, craiFile ->
            tuple(sampleID1, cramFile, craiFile)
        }
        .set { meta_aln_index_for_calcMeanDepth }
}

if (params.fastq) {
    if (!params.fastqInput) {
        log.warn "FASTQ provided but --fastqInput not set. The pipeline expects CRAM by default."
    }
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

    sampleid_R1.join(sampleid_R2)
        .set { read_pairs_ch }
}

/* -----------------------------------------------------------------
   SUBWORKFLOWS / MODULES
   ----------------------------------------------------------------- */
include { 
    inputFiles_symlinks_cram
    inputFiles_cramCopy
    samtools
    qualimap
    fastqc_bam
    collectWGSmetrics
    multiQC
    vntyper_newRef
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
   CALCULATING DEPTH
   ----------------------------------------------------------------- */
/**
 * We'll do the older "automatic binding" style for the process,
 * so we do NOT call it from the workflow block.
 * Instead, we define from meta_aln_index_for_calcMeanDepth in `input:`,
 * and emit into meanDepthChannel in `output:`.
 */
process calcMeanDepth {
    tag "$sampleID"

    input:
        tuple val(sampleID), file(cram), file(crai) \
            from meta_aln_index_for_calcMeanDepth

    output:
        tuple val(sampleID), file("mean_depth.txt") \
            into meanDepthChannel

    script:
    """
    samtools depth -a ${cram} | awk '{sum += \$3} END {if (NR>0) print sum/NR; else print 0}' > mean_depth.txt
    """
}

/* -----------------------------------------------------------------
   DEFINE A GLOBAL VARIABLE FOR MEAN DEPTH
   ----------------------------------------------------------------- */
def meanDepthSummary = ''

/* -----------------------------------------------------------------
   MAIN WORKFLOW
   ----------------------------------------------------------------- */
workflow {
    /**
     * Because we used "input: from meta_aln_index_for_calcMeanDepth" in calcMeanDepth,
     * Nextflow automatically triggers that process when data arrives on that channel.
     * We do NOT call `calcMeanDepth(...)` in the workflow block.
     */

    // Subscribe to the meanDepthChannel to build the summary
    if (params.cram) {
        meanDepthChannel.collect().subscribe { results ->
            meanDepthSummary = results.collect { sample, depthFile ->
                def depth = depthFile.text.trim()
                return "${sample}: ${depth}X"
            }.join('\n')
        }
    }

    // Example panel logic
    if (!params.panel || params.panel == 'WGS_CNV' || params.panel == 'NGC') {
        // If we have FASTQ input
        if (params.fastqInput || params.fastq) {
            SUB_PREPROCESS(read_pairs_ch)
            if (!params.preprocessOnly) {
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
        }
        // If we have CRAM input
        else if (params.cram) {
            if (!params.copyCram) {
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
            else {
                inputFiles_symlinks_cram(meta_aln_index)
                inputFiles_cramCopy(meta_aln_index)
                if (!params.skipVariants) {
                    SUB_VARIANTCALL_WGS(inputFiles_cramCopy.out)
                }
                if (!params.skipSV) {
                    SUB_CNV_SV(inputFiles_cramCopy.out)
                }
                if (!params.skipSTR) {
                    SUB_STR(inputFiles_cramCopy.out)
                }
                if (!params.skipSMN) {
                    SUB_SMN(inputFiles_cramCopy.out)
                }
            }
        }
    }

    if (params.panel && params.panel != "WGS_CNV" || params.panel != 'NGC') {
        if (params.fastqInput || params.fastq) {
            SUB_PREPROCESS(read_pairs_ch)
            SUB_VARIANTCALL(SUB_PREPROCESS.out.finalAln)
            if (params.panel == "MV1") {
                vntyper_newRef(read_pairs_ch)
            }
        }
        else if (params.cram) {
            inputFiles_symlinks_cram(meta_aln_index)
            SUB_VARIANTCALL(meta_aln_index)
        }
    }
}

/* -----------------------------------------------------------------
   COLLECT SAMPLE NAMES FROM CRAM
   ----------------------------------------------------------------- */
def sampleNamesList = []

if (params.cram) {
    sampleID_cram
        .map { it[0] }
        .collect()
        .subscribe { allSampleIDs ->
            sampleNamesList = allSampleIDs.unique()
        }
}

/* -----------------------------------------------------------------
   ON COMPLETE: send email with sample names, etc.
   ----------------------------------------------------------------- */
workflow.onComplete {
    def currentYear = new Date().format('yyyy')

    // Optional: lnx02 emailing
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
    def sampleNamesString = sampleNamesList.join('\n')

    if (!params.nomail && workflow.success && workflow.duration > 3) {
        if (user in ["mmaj", "raspau"]) {
            def sequencingRun = params.cram
                ? new File(params.cram).getName().take(6)
                : params.fastq
                    ? new File(params.fastq).getName().take(6)
                    : 'Not provided'

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
                          |
                          |Samples included in the pipeline:
                          |${sampleNamesString}
                          |
                          |Mean Depth per Sample:
                          |${meanDepthSummary}
                          """.stripMargin('|')

            def recipients = 'Rasmus.Hojrup.Pausgaard@rsyd.dk'

            if (params.server == 'lnx01') {
                sendMail(to: recipients, subject: 'CRAM-based pipeline Update', body: body)
            }
            else if (params.server == 'lnx02') {
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
