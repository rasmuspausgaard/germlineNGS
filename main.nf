#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Full Nextflow DSL2 script, focusing on CRAM inputs (optionally FASTQ).
 * Includes a mosdepth coverage process called `doMosdepthCoverage`
 * that is invoked in the main workflow on `meta_aln_index`.
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

// Preset parameters
params.gatk                 = params.gatk                ?: null
params.copyCram             = params.copyCram            ?: null
params.single               = params.single              ?: null
params.server               = params.server              ?: "lnx01"
params.genome               = params.genome              ?: "hg38"
params.outdir               = params.outdir              ?: "${launchDir.baseName}.Results"
params.rundir               = params.rundir              ?: "${launchDir.baseName}"

// For coverage
params.reference = params.reference ?: '/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa'

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
    // CRAM + CRAI
    cramfiles = "${params.cram}/${reads_pattern_cram}"
    craifiles = "${params.cram}/${reads_pattern_crai}"

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

    // Join CRAM + CRAI => meta_aln_index
    sampleID_cram.join(sampleID_crai).set { meta_aln_index }
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

    sampleid_R1.join(sampleid_R2).set { read_pairs_ch }
}

// If purely CRAM-based (no FASTQ):
if (params.cram && !params.fastq) {
    // meta_aln_index is the final input channel for alignment-based steps
}

// If purely FASTQ-based (no CRAM):
if (!params.cram && params.fastq) {
    // read_pairs_ch is the final input channel for alignment steps
}

/* -----------------------------------------------------------------
   SUBWORKFLOWS / MODULES
   ----------------------------------------------------------------- */
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
    SUB_PREPROCESS
    SUB_VARIANTCALL
    SUB_VARIANTCALL_WGS
    SUB_CNV_SV
    SUB_STR
    SUB_SMN
} from "./modules/modules.dna.v1.nf"

/* -----------------------------------------------------------------
   PROCESS: doMosdepthCoverage (mosdepth)
   Takes [sampleID, cramFile, craiFile] => outputs [sampleID, coverageValue]
   ----------------------------------------------------------------- */
process doMosdepthCoverage {
    input:
        tuple val(sampleID), path(cramFile), path(craiFile)

    tag { sampleID }

    output:
        tuple val(sampleID), stdout

    script:
        """
        sampleID='${sampleID}'
        cram='${cramFile}'
        crai='${craiFile}'
        ref='${params.reference}'

        if [ ! -f "\${crai}" ]; then
            echo "No CRAI index found for \${cram}. Indexing..."
            samtools index "\${cram}"
        fi

        prefix="\${sampleID}_mosdepth"
        singularity run -B /data/:/data/,/lnx01_data2/:/lnx01_data2/ \\
          /lnx01_data2/shared/testdata/mosdepth.sif \\
          mosdepth -n --fast-mode -t 4 --fasta "\${ref}" \\
          "\${prefix}" "\${cram}"

        summaryFile="\${prefix}.mosdepth.summary.txt"
        if [ -f "\${summaryFile}" ]; then
          grep '^total' "\${summaryFile}" | awk '{print \$4}'
        else
          echo "0"
        fi
        """
}

/* -----------------------------------------------------------------
   Collect sample names from CRAM
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
   MAIN WORKFLOW
   ----------------------------------------------------------------- */
// A normal Groovy list to hold coverage data from coverage process
def coverageList = []

workflow {
    // First block: run if panel is null or panel == 'WGS_CNV' or panel == 'NGC'
    if (!params.panel || params.panel == 'WGS_CNV' || params.panel == 'NGC') {
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
            // COVERAGE CALCULATION
            def coverageResults = doMosdepthCoverage(meta_aln_index)
            coverageResults.subscribe { result ->
                coverageList << result
                println "Coverage for sample '${result[0]}': ${result[1]}"
            }
        }
    } // End of first if-block

    // Second block: run if panel is set, not 'WGS_CNV', and not 'NGC'
    if (params.panel && params.panel != "WGS_CNV" && params.panel != "NGC") {
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
} // End of workflow

/* -----------------------------------------------------------------
   ON COMPLETE: send email with sample names, coverage, etc.
   ----------------------------------------------------------------- */
workflow.onComplete {
    def currentYear = new Date().format('yyyy')
    def sampleNamesString = sampleNamesList.join('\n')
    def coverageSummary = coverageList
        .collect { tuple -> "${tuple[0]}: ${tuple[1].trim()}" }
        .join('\n')

    if (!params.nomail && workflow.success && workflow.duration > 3 && user in ["mmaj", "raspau"]) {
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
                      |${obsSampleMessage}
                      |
                      |Samples included in the pipeline:
                      |${sampleNamesString}
                      |
                      |Coverage results:
                      |${coverageSummary}
                      """.stripMargin('|')

        def recipients = 'Rasmus.Hojrup.Pausgaard@rsyd.dk'
        if (params.server == 'lnx01') {
            sendMail(to: recipients, subject: 'CRAM-based pipeline Update', body: body)
        }
        else if (params.server == 'lnx02') {
            def ipFilePath = '/lnx01_data2/shared/testdata/test_scripts/ip_file'
            def ip = new File(ipFilePath).exists() ? new File(ipFilePath).text.trim() : ''
            def emailCommand = "ssh ${ip} 'echo \"${body}\" | mail -s \"CRAM-based pipeline Update\" ${recipients}'"
            def proc = ['bash', '-c', emailCommand].execute()
            proc.waitFor()
            if (proc.exitValue() != 0) {
                println("Error sending email from lnx02: ${proc.err.text}")
            } else {
                println("Email successfully sent from lnx02.")
            }
        }
    }

    // Delete work directory if not keeping
    if (!params.keepwork) {
        println("Deleting work directory: ${workflow.workDir}")
        def deleteWorkDirCommand = "rm -rf ${workflow.workDir}".execute()
        deleteWorkDirCommand.waitFor()
        if (deleteWorkDirCommand.exitValue() != 0) {
            println("Error deleting work directory: ${deleteWorkDirCommand.err.text}")
        }
    }
}
