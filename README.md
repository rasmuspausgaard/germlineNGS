# KGA Vejle DSL2 WGS pipeline, v1

## Vejledning

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
