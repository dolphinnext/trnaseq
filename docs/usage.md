# dolphinnext/trnaseq: Usage

## Introduction
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run dolphinnext/trnaseq -profile docker --DOWNDIR /path/to/save/trnaseq --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build human_hg38_gencode_v30
```

If you're running for the first time, you need to enable `--run_checkAndBuild` paramater as follows:

```bash
nextflow run dolphinnext/trnaseq -profile docker --DOWNDIR /path/to/save/trnaseq --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build human_hg38_gencode_v30 --run_checkAndBuild 'yes'
```


This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results 
.nextflow_log   # Log file from Nextflow
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version. In order to download latest version of the pipeline you need to run following command:

```bash
nextflow pull dolpinnext/trnaseq
```

## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from Dockerhub: [`dolphinnext/trnaseq`](http://hub.docker.com/r/dolphinnext/trnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*.fastq' --mate 'single'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character


### `--mate`
You can specify the single-end data by entering mate parameter as 'single'. For example:

```bash
--reads 'path/to/data/sample_*.fastq' --mate 'single'
```


## Reference genomes

### `--genome_build` 
To run the pipeline, you must specify which to use with the `--genome_build` flag.

List of genomes that are supported are:

* Human
  * `--genome_build human_hg19`
  * `--genome_build human_hg19_eColitK_suppressor_oligostd`
  * `--genome_build mouse_mm10`
  * `--genome_build mouse_mm10_suppressor_oligostd`
  * `--genome_build custom`
  

Note: For new genome requests, please send e-mail to UMMS-Biocore(biocore@umassmed.edu).

```

## Adapter Removal
If specific Adapter Removal is required, you can enable trimmomatic and enter the adapter sequence. 

```bash
To enable adapter_removal: 
--run_Adapter_Removal "yes"

--Adapter_Trimmer_Quality_Module_Adapter_Removal.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

--Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence [string]
# You can enter a single sequence or multiple sequences in different lines. Reverse sequences will not be removed.

--Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length [int @default:10]
# Specifies the minimum length of reads to be kept

--Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches [int @default:1]
# Specifies the maximum mismatch count which will still allow a full match to be performed

--Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold [int @default:30]
# Specifies how accurate the match between the two -adapter ligated- reads must be for PE palindrome read alignment

--Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold [int @default:5]
# Specifies how accurate the match between any adapter etc. sequence must be against a read.

--Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped [@options:"yes","no" @default:"yes"]
# Discard_non_clipped sequences (keep only sequences which contained the adapter)
```

## Trimmer
Optianally, you can trim your reads by defining trimming lenghts as shown at below: 

```bash

--run_Trimmer [@options:"yes","no" @default:"no"]
# Enables Trimmer by setting this parameter as "yes"

--Adapter_Trimmer_Quality_Module_Trimmer.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

For Single End Reads  : 
--Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads "single"
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime [int]

For Paired End Reads  : 
--Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads "pair"
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2 [int]
```

## Quality Filtering
Optianally, you can trim your reads based on their quality. Trimmomatic works on both paired-end and single ended data. Alternatively fastx option (fastx_toolkit) could be used for single reads. 

```bash
To use Trimmomatic  : 
--run_Quality_Filtering "yes"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.tool "trimmomatic"

--Adapter_Trimmer_Quality_Module_Quality_Filtering.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

--Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size [int @default:10]
# Performs a sliding window trimming approach. It starts scanning at the 5' end and clips the read once the average quality within the window falls below a threshold (=required_quality).

--Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming [int @default:15]
# Specifies the average quality required for window trimming approach

--Adapter_Trimmer_Quality_Module_Quality_Filtering.leading [int @default:5]
# Cut bases off the start of a read, if below a threshold quality

--Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing [int @default:5]
# Cut bases off the end of a read, if below a threshold quality

--Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen [int @default:36]
# Specifies the minimum length of reads to be kept
```

```bash
To use fastx_toolkit  : 
--run_Quality_Filtering "yes"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.tool "fastx"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality [int @default:20]
# Minimum quality score to keep reads

--Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent [int @default:100]
# Minimum percent of bases that must have entered minQuality
```

```bash
Process Parameters for UMIextract :
params.UMIextract.single_or_paired_end_reads = "single" //* @dropdown @options:"single","pair" 
params.UMIextract.barcode_pattern1 = "(?P<umi_1>(CAAGATCGGAAGAGCACACGTCTGAA)){s<=1}.+" //* @input
params.UMIextract.barcode_pattern2 = "" //* @input
params.UMIextract.UMIqualityFilterThreshold = "13" //* @input @description:"Quality (phred quality score) cutoff for UMI. Default is 13, that is UMI with qualities >= 13 will be kept."
params.UMIextract.phred = 33 //*  @dropdown @options:"33","64" @description:"Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data"
params.UMIextract.remove_duplicates_based_on_UMI = "true" //* @checkbox @description:"Removes duplicate reads by checking UMI."
```

```bash
Process Parameters for mimseq :
params.mimseq_groups.control_group_name = "ctrl" //* @input @description:"Control/Wild-type group name"
params.mimseq_groups.control_group = "Ctrl" //* @input @description:"Comma-separated list of sample names (please don't include the extension of the file e.g. for control.rep1.fq -> use control.rep1)" 
params.mimseq_groups.treatment_group_name = ["exp"] //* @title:"Treatment Groups" @input @description:"Treatment group name for figures"
params.mimseq_groups.treatment_group = ["Expr"] //*  @input @description:"Comma-separated list of samples (please don't include the extension of the file) (e.g. for first group: treat1.rep1, treat1.rep2 and for second group click add button and enter: treat2.rep1, treat2.rep2)" 
params.mimseq.experiment_name = "experiment1" //* @input @description:"Name of experiment. Note, output files and indices will have this as a prefix"
params.mimseq.mimseq_parameters = "--cluster --cluster-id 0.95  --threads 10 --min-cov 2000 --max-mismatches 0.1 --max-multi 1 --remap --remap-mismatches 0.075" //* @input @description:"Mimseq parameters"
```

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). 

Note - you can use this to override pipeline defaults.