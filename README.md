# TumorCNV
# Introduction 
TumorCNV is a tool designed to jointly detecting of germline and somatic copy number events from WGS data of the matched tumor-normal sample pair.
# Installation
The easiest way to get TumorCNV is to download the binary distribution from the TumorCNV github release page. Alternatively, you can build TumorCNV from source with gradle.
1. git clone --recursive https://github.com/yongzhuang/TumorCNV.git
2. Install gradle build tool (https://gradle.org/)
3. cd TumorCNV 
4. gradle build   
You'll find the executable jar file in TumorCNV/build/libs/. 

If you want to run TumorCNV, you'll need:
1. Install Java SE Development Kit 8
2. Install R (Rscript exectuable must be on the path)
3. Install Runiversal (https://cran.r-project.org/web/packages/Runiversal/index.html), VGAM (https://cran.r-project.org/web/packages/VGAM/index.html) and qcc(https://cran.r-project.org/web/packages/qcc/) package in R

# Running
usage: java -jar TumorCNV.jar [OPTIONS]
1. preprocess  
   This option is used to extract the information from the normal and tumor BAM files.

   usage: java -jar TumorCNV.jar preprocess [OPTIONS]

   -referenceSequenceFile  <FILE>   reference genome file (required)  
   -normalVCFFile <FILE>   normal sample's vcf file (optional)   
   -normalBAMFile <FILE>   normal sample's bam file (required)  
   -tumorBAMFile  <FILE>   tumor sample's bam file (required)  
   -mappabilityFile  <FILE>   mappability file (required)  
   -outputPrefix  <FILE>    prefix of output file (required)  
   -windowSize <INT> window size (optional, default 500)  
   -minMappingQuality   <INT> minimum mapping quality (optional, default 1)  
   -minBaseQuality   <INT> minimum base quality (optional, default 20)  

2. call  
   This option is used to jointly call germline and soamtic CNVs.  

   usage: java -jar TumorCNV.jar call [OPTIONS]  

   -rdFile  <FILE>   read depth file (required)  
   -afFile  <FILE>   allele frequency file (optional)  
   -mappabilityFile  <FILE>   mappability file (required)  
   -outputFile <FILE>   prefix of toutput file (required)  
   -exclude <FILE>   exclude regions  
   -transitionProb   <FLOAT>  transition probability of different states (optional, default 0.00001)  
   -minMappability   <FLOAT>  minimum mappability of window (optional, default 0.3)  
   -minDisatance  <INT> minimum distance to merge adjacent CNVs (optional, default 10000)  
   -purity  <FLOAT>  tumor purity (optional, default 1.0)  
   -ploidy  <INT> tumor ploidy (optional, default 2)  
   -outlier <FLOAT>  the percentage of outliers (optional, default 0.1)  
   -nt   <INT> number of threads (optional, default 1)  

# Example

The sample data sets and script can be found at TumorCNV/example, but the users need to download the reference genome file (ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz) and mappability file (http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig).

# Contact 
   yzhuangliu@gmail.com
