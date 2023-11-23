# KG Vejle Germline pipeline

## General info:
This script is used for KG Vejle germline analysis of the following NGS designs:
1. Panels (AV1 and CV5 currently supported)
2. WES
3. WGS CNV (simple pipeline for WGS CNV)
4. WGS (full analysis)

The script takes mapped data (CRAM files) OR raw data (fastq) as input.
Hg38 (v3) assembly is used by default. Hg19 is no longer supported. 

# Usage

The tools used and output generated depends on how the script is used, i.e. if it is used for panels or WGS. See below for instructions.


## Panels (AV1, WGS_CNV, CV5, WES)
if fastq is used as input (using --fastq /path/to/fastq/) the script will perform read preprocessing, alignment, variantcalling and joint genotyping.

If CRAM is used as input (using --cram /path/to/cram/) the script will perform variantcalling and joint genotyping.

Creates a merged (joint Genotyped) VCF for all analyzed samples.
For AV1 panels and WES (ALM and ONK), joint genotyping is also run per subpanel.
for AV1 panels, SpliceAI is run by default (spliceAI can be skipped with "--skipSpliceAI", see usage below

NOTE: CRAM should be used as input, if possible.

NOTE: This script can be run from both servers (kga01 and lnx01) when analyzing paneldata. The script assumes the lnx01 server is used by default. If this script is run from kga01, make sure to set the "--server kga01" parameter.

## WGS

WGS analysis requires a tab-delimited samplesheet containing 4 columns without headerline in this specific order:

caseID/projectID, NPN, Relation, SampleStatus

Example samplesheet for standard trio:

johnDoe 123456789012    index   affected 

johnDoe 234567890123    mater   normal

johnDoe 345678901234    pater   normal

The above information can usually be extracted directly from the sample overview excel file

If the inputdata (FastQ or CRAM) have been transferred to the data archive (which it is by default), the script will automatically find the relevant inputdata  and create symlinks for them in the output (results) directory, if --samplesheet /path/to/samplesheet/ is used.

The user can point to a specific folder containing raw data (FastQ) using the --fastq option or alignment data (CRAM) using the --cram option
This is only needed if input data (FastQ or CRAM) exists outside the data archive (e.g. if data are in personal folders), or if the script is run without samplesheet.

## Default settings:

By default, the pipeline uses the hg38 assembly and consists of the following sub workflows:

1. Preprocessing if FastQ is used as input (fastq --> CRAM or BAM).
2. Variantcalling (GATK Haplotypecaller + jointGenotyping)
3. QC module 
4. CNV and SV calling (CNVkit, Manta, Lumpy and Tiddit + merging of SV calls using SVDB)
5. STR analysis (ExpansionHunter and Stripy)
6. SMNcaller (SMNCopyNumberCaller)


## Most common use cases:
This pipeline uses hg38 (KG Vejle V3) by default.

Analyze all WGS data in folder (no samplesheet), starting with CRAM files:

    nextflow run KGVejle/wgsPipeline -r main --cram /path/to/cram/ 

Analyze only "WGS CNV" samples, starting with CRAM files: 

    nextflow run KGVejle/wgsPipeline -r main --cram /path/to/cram/ --panel WGS_CNV 


Analyze samples in samplesheet, starting with fastq:
   
    nextflow run KGVejle/wgsPipeline -r main --samplesheet /path/to/samplesheet/ --fastqInput 


Analyze samples in samplesheet, starting with CRAM:
   
    nextflow run KGVejle/wgsPipeline -r main --samplesheet /path/to/samplesheet/

Run the script with --help to see available options and default parameters:

    nextflow run KGVejle/wgsPipeline -r main --help




