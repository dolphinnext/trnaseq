tRNA-Seq Pipeline is an automated analysis pipeline for the quantitation and analysis of tRNA expression and modifications. It accepts an adapter sequence for adapter removal and UMI pattern for UMI extraction. 

Main Features:

1. Cluster tRNAs, index modifications, and perform SNP-tolerant read alignment with GSNAP
2. Deconvolve cluster aligned reads back into unique tRNA transcript-level reads
3. Calculate coverage information and plots (useful for QC)
4. Quantify expression
5. Calculate tRNA differential expression with DESeq2.
6. Analyze functional tRNA pools and tRNA completeness via 3'-CCA analysis
7. Comprehensive modification quantification and misincorporation signature analysis


This pipeline adapted from following study: <a class="link-underline" href="https://doi.org/10.1016/j.molcel.2021.01.028" target="_blank">Paper</a>  and <a class="link-underline" href="https://github.com/nedialkova-lab/mim-tRNAseq" target="_blank">Code</a> 


##### Citation:

* If you use DolphinNext in your research, please cite: 
Yukselen, O., Turkyilmaz, O., Ozturk, A.R. et al. DolphinNext: a distributed data processing platform for high throughput genomics. BMC Genomics 21, 310 (2020). https://doi.org/10.1186/s12864-020-6714-x
* Behrens et al., 2021, High-resolution quantitative profiling of tRNA abundance and modification status in eukaryotes by mim-tRNaseq. Molecular Cell 81, 1–14 https://doi.org/10.1016/j.molcel.2021.01.028

##### Steps:

1. For Quality Control, we use FastQC to create QC outputs. 
2. Adapter Removal: Adapter removal is performed by trimmomatic.
3. There are optional read quality filtering (trimmomatic) and read quality trimming (trimmomatic) processes available after adapter removal. 
4. Optionally UMI-extraction and duplicate read are removed.
5. Cluster tRNAs, index modifications, and perform SNP-tolerant read alignment with GSNAP.
6. Calculate coverage information and plots (useful for QC).
7. Quantify expression.
8. Calculate tRNA differential expression with DESeq2.
9. Analyze functional tRNA pools and tRNA completeness via 3'-CCA analysis.
10. Comprehensive modification quantification and misincorporation signature analysis.

##### Inputs:

* **Reads**: Specify the location of your input FastQ file. <a class="link-underline" href="https://dolphinnext.readthedocs.io/en/latest/dolphinNext/quick.html#adding-files" target="_blank">Need Help?</a>
* **Adapter Sequence**: Please enter the adapter sequence(s) in the settings of `run_Adaper_Removal`.
* **Umi extract**: Please use following pattern for umi extract. Barcode pattern 1:`.+(?P<umi_1>(NNNCGANNNTACNNN)|(NNNATCNNNAGTNNN)){s<=1}.+` and enable `remove_duplicates_based_on_UMI` option
 
##### Program Versions:
  - python=3.7
  - mimseq=0.3.4
  - umi_tools=1.1.1
  - trimmomatic=0.39
  - fastqc=0.11.8
  - fastx_toolkit=0.0.14

##### Run through DolphinNext User Interface:

To start using the dolphinnext/trnaseq Pipeline please go to <a class="link-underline" href="https://dolphinnext.umassmed.edu/index.php?np=1&id=744" target="_blank">DolphinNext Web page</a> and click run button.

##### Run through Command Line:

To install and start using the dolphinnext/trnaseq pipeline by using command line, please follow these steps: <a class="link-underline" href="https://github.com/dolphinnext/trnaseq/blob/1.0/docs/local.md" target="_blank">Installation</a> .