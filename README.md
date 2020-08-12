
# functional and ecological determinants of evolutionary dynamics in crAss-like phages
contains Python code.  
It could be merged with Linux bash scripts, which call the Python programs and 
therefore also contain clear examples on parameter usage

to do
- installation instructions
- all Conda dependencies (refer to generated file by Conda?)

# Summary 


# Examples 

Accompanying shell scripts are in the mgx folder of:
https://github.com/resharp/phages_scripts

# Installation
this repository can be best put in local directory:

source/phages

the shell scripts reference the Python scripts by this path 

# Annotation pipeline

- AnnotateCrassGenomes.py
    - All results of annotation pipeline are integrated, input:
        - HMM search results for Yutin HMM profiles
        - searched all genes against HMM profiles of pVOG database (Grazziotin 2017, date accessed 14/1/2020).
        - to do: (all-to-all pairwise blastp search between the genes)
   
*help scripts*
- create_proteins_from_gbk.py
    - We extracted the proteins from NC_024711.1 (genus 1)
    - output: crassphage_refseq.proteins.faa
    - used for Yutin hmm-searching, pairwise blastp and pVOG searches 
- extract_positions_from_prodigal_predictions.py
    - generates gene position files (as input for DiversiTools) from the prodigal gene predictions 
- create_categories_from_gbk.py
    - For genus 1 we created the gene positions from the GenBank file
     
# measures calculated per sample

- CalcDiversiMeasures.py
    - aggregate the output of DiversiTools 
    - calculated the following measures per gene per sample 
    - (also per sliding window/bin per sample):
        -	Mean and standard deviation (dev) of coverage
        -	Mean and standard dev of top amino acid (AA) count, second AA count and third AA count
        -	Sum of synonymous and non-synonymous counts
        -	Total number of SCPs 
        -	dN/dS
        -	entropy 
        -	pN/pS
    
- MakeGenePlots.py
    - aggregates the output of CalcDiversiMeasures.py. It contains
        - analysis for one reference genome 
        - analysis for gene families based on multiple ref genomes
    - statistical testing and plotting
    - (uses matplotlib.pyplot and seaborn)

creates figures 2 and 3 of main text      
and supplementary tables [to do]        
        
# measures calculated across samples
- CalcCodonMeasures.py
    - For cross-sample measures all reads of the samples were first stacked on top of each other, ... 
- MakeDiversityPlots.py 
    - ... after which we calculated the following measures for every codon position
        -   the overall coverage
        -	the total number of possible codon types in all samples (equal to one when no variation is present; we accepted a codon polymorphism (SCP) if it occurred at least in three reads  in all samples and at a minimum of 1% of the total depth at that codon position.)
        -	the entropy 
        -	pN/pS   

creates figure 4 of main text  
and supplementary tables [to do]

*help files*
- MakeSynProbabilitiesTable.py
    - creates table for codon bias
    - output: codon_syn_non_syn_probabilities.txt

# macro diversity
- MakeSamplePlots.py
    - calculated  macro diversity per sample based on the sample/genome mapping statistics produced by samtools idxstats (Li 2009) on the bam files

creates figure 1 and 5 of main text

# Summary 
- MakeGeneSummary.py
    - integrates different statistics for genus 1 genes
    and gene families in one table
    - also contains a pangenome analysis based on the hmm hits from the protein family profiles
 
# extra stuff
Scripts used for selecting samples
- ExtractSraMetadata.py
    - parses metadata files in MGXDB database
    - extract samples containg a certain taxonid from metadata files
User for proof of concept:
- DownloadMgnifySamples.py
    - download metadata for samples and runs based on some properties
        - human gut fecal samples
        - only metagenomics data sets (and being able to filter on this in runs)
        - also add information from analysis:
        - number of reads
- download_fastq.py
    - download ...    
 
    
# clustering proteins and prophages from prophage predictions by VirSorter
This is obsolete, and should be moved to different GIT repository:
- MclClusterEvaluation.py
- PhageSharedContent.py
- PhageClusterEvaluation.py
- VirsorterStats.py
- VirsorterStatsTest.py
    - also depends on MyTestCase.py
- PhageContentPlots.py
- MyTestCase.py
    - parent class of all test case classes (contains generic code for logging start of function)

Code for learning Pandas and python:
- IctvStats.py
    - code for getting statistics from: 
    - ICTV Master Species List 2018b.v1.csv
        - downloaded from https://talk.ictvonline.org/
    - prokaryotes.cv
- IctvStatsTest.py
    - accompanying tests


