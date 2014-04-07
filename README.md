# Supplemental scripts for analysis of the planarian regeneration transcriptome

Publication: [Damian Kao, Daniel Felix, Aziz Aboobaker. The planarian regeneration transcriptome reveals a shared but temporally shifted regulatory program between opposing head and tail scenarios. BMC Genomics 2013, 14:797. doi:10.1186/1471-2164-14-797](http://www.biomedcentral.com/1471-2164/14/797)

9 scripts are included in this repository for various processes of transcriptome consolidation and tag count analysis:  

* **cluster_hierarchical.py** - Performs hierarchical clustering on libraries using complete linkage and correlation distance.  
* **consolidated_fusionScore.py** - Calculates fusion scores from tab delimited BLAST output. More information on BLAST output format is commented in the script.  
* **consolidated_getRep.py** - Using cluster information along with nucleotide and longest ORF length data, get the sequence with the longest ORF.  
* **consolidated_multiSource.py** - from parsed CAP3 output, get clusters that came from at least 2 different sources.  
* **consolidated_parseCAP3.py** - parse CAP3 stdout.  
* **counts_filterLow.py** - Filter a tab delimited file of tag counts so all transcripts with less than 20 tag counts in all libraries are taken out.    
* **counts_topExpression.py** - List transcripts that take up more than 1% of total reads in more than 2 libraries.  
* **fa_lengths.py** - Generate a tab delimited list of nucleotide lengths where column 1 is transcript id and column 2 is transcript length.  
* **fa_longestORF.py** - Takes EMBOSS' getorf .fasta output and generates a tab delimited list of longest ORF lengths where column 1 is transcript id and column2 is longest ORF length.  
