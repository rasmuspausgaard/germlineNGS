# wgsPipeline
Standard pipeline, WGS KG Vejle

## General info:
The only requirement is a samplesheet containing 4 columns without headerline in this specific order:
famID/projektNavn, NPN, Relation, SampleStatus

Example samplesheet for standard trio:
johnDoe 123456789012    index   affected
johnDoe 234567890123    mater   normal
johnDoe 345678901234    pater   normal

The above information can usually be extracted directly from the sample overview excel file

If the inputdata (FastQ or CRAM) have been transferred to the data archive (which it is by default), the script will automatically find the relevant inputdata  and create symlinks for them in the output (results) directory.

The script will automatically look for FastQ or CRAM files in subfolders at the dataArchive. 

The user can point to a specific folder containing raw data (FastQ) using the --fastq option or alignment data (CRAM) using the --cram option
This is only needed if input data (FastQ or CRAM) exists outside the data archive (e.g. if data are in personal folders), or if the script is run without samplesheet.

## Default settings:

By default, the pipeline uses the hg38 assembly and consists of the following sub workflows:

1. Preprocessing if FastQ is used as input (fastq --> CRAM or BAM).
2. Variantcalling (GATK Haplotypecaller + jointGenotyping)
3. QC module 
4. CNV and SV calling (CNVkit, Manta, Lumpy and Tiddit)
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


