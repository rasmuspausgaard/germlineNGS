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
// Example intervals (if needed):
// params.intervals_list    = "/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/..."


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
   MAIN WORKFLOW
   ----------------------------------------------------------------- */
workflow {
    /*
     * Panel logic for WGS_CNV, NGC, or if panel is null => WGS,
     * or if panel is set => do subworkflow for that panel, etc.
     */

    if (!params.panel || params.panel == 'WGS_CNV' || params.panel == 'NGC') {
        // If we have FASTQ input
        if (params.fastqInput || params.fastq) {
            // Align FASTQ => CRAM
            SUB_PREPROCESS(read_pairs_ch)

            if (!params.preprocessOnly) {
                // e.g. run variant calling, SV, STR, SMN
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
                // Symlink CRAM
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
                // Physically copy CRAM
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
        // Panel logic
        if (params.fastqInput || params.fastq) {
            SUB_PREPROCESS(read_pairs_ch)
            SUB_VARIANTCALL(SUB_PREPROCESS.out.finalAln)

            if (params.panel == "MV1") {
                // For example, run vntyper for MV1
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

// Only do this if we have CRAM input
if (params.cram) {
    // sampleID_cram emits (sampleID, cramFile)
    // We'll map to just sampleID, collect them all, and store in sampleNamesList
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
    // Determine the current year dynamically
    def currentYear = new Date().format('yyyy')
    
    // Read IP address from file
    def ipFilePath = '/lnx01_data2/shared/testdata/test_scripts/ip_file'
    def ip = ""

    if (new File(ipFilePath).exists()) {
        println("IP file exists. Reading IP address.")
        ip = new File(ipFilePath).text.trim()
        println("IP address read from file: ${ip}")
    } else {
        println("Error: IP address file not found at ${ipFilePath}")
        return
    }

    // Kør AV1-kontrol (chr2:47414420) på VarSeq-VCF
    def av1PositionMsg = ""
    // Correctly set the outputDir
    def outputDir = "${launchDir}/${launchDir.baseName}.Results"

    if (params.panel == 'AV1' && workflow.success) {
        try {
            def variantsDir = new File(outputDir, "Variants")
            if (variantsDir.exists()) {
                // Find alle AV1 VarSeq-VCF’er for denne kørsel
                def vcfFiles = variantsDir.listFiles()?.findAll { f ->
                    f.name.endsWith("AV1_ALL.hg38.V3.merged.for.VarSeq.vcf.gz") ||
                    f.name.endsWith("AV1_ALL.hg38.V3.merged.for.VarSeq.vcf")
                } ?: []

                if (!vcfFiles.isEmpty()) {
                    def allHits = []
                    vcfFiles.each { f ->
                        // Hvis der en dag bliver brug for 2:… i stedet, kan vi udvide, men dine filer bruger chr2
                        def cmd = "bcftools query -r chr2:47414420 -f '[%SAMPLE\\t%GT\\n]' '${f.absolutePath}'"
                        def proc = ['bash', '-c', cmd].execute()
                        proc.waitFor()

                        if (proc.exitValue() == 0) {
                            def lines = proc.in.text.readLines().findAll { line ->
                                def toks = line.split("\\t")
                                toks.size() >= 2 && (toks[1] == "0/0" || toks[1] == "./.")
                            }
                            lines.each { hit ->
                                allHits << "${f.name}: ${hit}"
                            }
                        } else {
                            println "Error running bcftools on ${f.absolutePath}: ${proc.err.text}"
                        }
                    }

                    if (!allHits.isEmpty()) {
                        av1PositionMsg = """
                        
                        AV1 check at chr2:47414420
                        --------------------------
                        Samples with GT 0/0 or ./.:
                        ${allHits.collect { " - ${it}" }.join("\n")}
                        """.stripIndent()
                    } else {
                        println "AV1 check: no 0/0 or ./. at chr2:47414420 for this run"
                    }
                } else {
                    println "AV1 check: no VarSeq VCF found in ${variantsDir.absolutePath}"
                }
            } else {
                println "AV1 check: Variants dir not found: ${variantsDir.absolutePath}"
            }
        }
        catch (Exception e) {
            println "Exception during AV1 chr2:47414420 check: ${e.message}"
        }
    }

    // Only send email if --nomail is not specified, the user is mmaj or raspau, and duration is longer than 5 minutes / 300000 milliseconds
    if (!params.nomail && workflow.duration > 300000 && workflow.success) {
        if (System.getenv("USER") in ["raspau", "mmaj"]) {
            def sequencingRun = params.cram ? new File(params.cram).getName().take(6) :
                               params.fastq ? new File(params.fastq).getName().take(6) : 'Not provided'

            // Checks if there are OBS samples in the cram folder
            def obsSampleMessage = ""
            if (params.panel == "AV1" && params.cram) {
                def cramDir = new File(params.cram)
                def obsSamples = cramDir.listFiles()?.findAll { it.name.contains("OBS") } ?: []
                if (obsSamples.size() > 0) {
                    obsSampleMessage = "\nTHERE IS AN OBS SAMPLE IN THIS RUN"
                }
            }

            def workDirMessage = params.keepwork ? "WorkDir: ${workflow.workDir}" : "WorkDir: Deleted"

            // Her bruger vi outputDir fra ovenfor

            def body = """\
            Pipeline execution summary
            ---------------------------
            Pipeline completed: ${params.panel}
            Sequencing run: ${sequencingRun}${obsSampleMessage}
            Duration: ${workflow.duration}
            Success: ${workflow.success}
            ${workDirMessage}
            OutputDir: ${outputDir}
            Exit status: ${workflow.exitStatus}
            ${obsSampleMessage}
            ${av1PositionMsg}
            """.stripIndent()

            // Construct the email sending command
            def subject = 'GermlineNGS pipeline Update'
            def recipients = 'Andreas.Braae.Holmgaard@rsyd.dk,Annabeth.Hogh.Petersen@rsyd.dk,Isabella.Almskou@rsyd.dk,Jesper.Graakjaer@rsyd.dk,Lene.Bjornkjaer@rsyd.dk,Martin.Sokol@rsyd.dk,Mads.Jorgensen@rsyd.dk,Rasmus.Hojrup.Pausgaard@rsyd.dk,Signe.Skou.Tofteng@rsyd.dk,Amalie.Schirmer.Ahlgreen.Larsen@rsyd.dk,Sara.Kaczor.Elbaek@rsyd.dk'

            if (params.server == 'lnx02') {
                // Use Nextflow's built-in sendMail function when on lnx01
                sendMail(to: recipients, subject: subject, body: body)
            }

            // Check if --keepwork was specified
            if (!params.keepwork) {
                // If --keepwork was not specified, delete the work directory
                println("Deleting work directory: ${workflow.workDir}")
                def deleteWorkDirCommand = "rm -rf ${workflow.workDir}".execute()
                deleteWorkDirCommand.waitFor()
                if (deleteWorkDirCommand.exitValue() != 0) {
                    println("Error deleting work directory: ${deleteWorkDirCommand.err.text}")
                }
            }

            // Move WGS.CNV from lnx02 to lnx01
            if (params.server == 'lnx02' && params.panel == 'WGS_CNV' && workflow.success) {
                def moveWGSCNVCommand = "mv ${launchDir} /lnx01_data2/shared/patients/hg38/WGS.CNV/${currentYear}/"
                def moveWGSCNVProcess = ['bash', '-c', moveWGSCNVCommand].execute()
                moveWGSCNVProcess.waitFor()

                if (moveWGSCNVProcess.exitValue() != 0) {
                    println("Error moving WGS_CNV files: ${moveWGSCNVProcess.err.text}")
                } 
            }

            // Move WES from lnx02 to lnx01
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
}




