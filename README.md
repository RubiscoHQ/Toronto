## Input Data

### Reference Genome
- Hg19

### ExAC data source
- ExAC release 0.3.1
	
	We have removed individuals affected by severe pediatric disease, so this data set should serve as a useful reference set of allele frequencies for severe disease studies.[[http://exac.broadinstitute.org/about](http://exac.broadinstitute.org/about)]  
	[ftp://ftp.broadinstitute.org/pub/ExAC\_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz](ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz)
	use annovar convert2annovar.pl to convert annovar input  
	` > perl convert2annovar.pl -format vcf4 ExAC.r0.3.1.sites.vep.vcf > ExAC.r0.3.1.sites.vep.vcf.avinput  `
	` > NOTICE: Finished reading 9362538 lines from VCF file  `
	` > NOTICE: A total of 9362318 locus in VCF file passed QC threshold, representing 9415611 SNPs (6265499 transitions and 3150112 transversions) and 780261 indels/substitutions  `
	` > NOTICE: Finished writing 9415611 SNP genotypes (6265499 transitions and 3150112 transversions) and 780261 indels`
- ExAC on non-TCGA samples based on release 0.3  
	- [http://annovar.openbioinformatics.org/en/latest/user-guide/download/](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)
- ExAC on non-Psychiatric disease samples based on release 0.3  
	- [http://annovar.openbioinformatics.org/en/latest/user-guide/download/](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)
- ExAC non-both  
	- Merge both non-psy and non-tcga. 

### Clinvar data source
- Clinvar\_20160302  
	- [http://annovar.openbioinformatics.org/en/latest/user-guide/download/](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)

### Annotation data source
- refGene (gene based annotation)  
	- [http://annovar.openbioinformatics.org/en/latest/user-guide/download/](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)
- dbNSFP version 3.0a   
	- (whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, DANN, fitCons, PhyloP and SiPhy scores from )  
	- [https://sites.google.com/site/jpopgen/dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP)  
	- [http://annovar.openbioinformatics.org/en/latest/user-guide/download/](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)

### Annotation software
- ANNOVAR  
	- All annotation database also download form Annovar  
	- [http://annovar.openbioinformatics.org/en/latest/](http://annovar.openbioinformatics.org/en/latest/)  

### Data processing
1. Using annovar to annotate respective input databases with refGene, dbNSFP, Clinvar  
	`perl table_annovar.pl inputfile humandb/ -buildver hg19 -out exportfile -protocol refGene,dbnsfp30a,clinvar_20160302 -operation g,f,f -nastring NA -remove --thread 4 --maxgenethread --onetranscript`
2. Filtered out mutations  
	- by type, remain only non-syndromic mutations  
	- by Pathogenicity: for ExAC mutations exclude those annotated “Pathogenic” by Clinvar; for Clinvar only remains “Pathogenic” mutations
3. Annotate ExAC vcf information using ExAC.r0.3.1.sites.vep.vcf
4. ExAC HighQuality Filter  
	- Generate HQ group based on ExAC high quality definition
	1. they were given a PASS filter status by VQSR,   
	2. at least 80% of the individuals in the dataset had at least depth (DP) \>= 10 and genotype quality (GQ) \>= 20 (i.e. AN_Adj \>= 60706*0.8*2 or 97130),   
	3. there was at least one individual harboring the alternate allele with depth \>= 10 and GQ \>= 20,  
	4. the variant was not located in the 10 1-kb regions of the genome with the highest levels of multi-allelic (quad-allelic or higher) variation.  
	##    chrom bin_start   bin_end freq  
	## 1     14 106330000 106331000  109  
	## 2      2  89160000  89161000  103  
	## 3     14 106329000 106330000   76  
	## 4     14 107178000 107179000   63  
	## 5     17  18967000  18968000   46  
	## 6     22  23223000  23224000   44  
	## 7      1 152975000 152976000   42  
	## 8      2  89161000  89162000   39  
	## 9     14 107179000 107180000   38  
	## 10    17  19091000  19092000   38
5. ExAC at least one Hom mutation  
	Generate 1hom group based on AC\_Hom \>= 1, which means at least one individual carry a homozygous allele.
6. retrench and generate input file  
	- extract useful informations for analysis. Score values, MAF (using ExAC popmax), Gene and so on.  
	- Defining 6 Freq group according to MAF:  
	I:   singleton, AN\_All == 1;  
	II:  2-10 alleles, 2 \<= AN\_All \<= 10;  
	III: 0.0001 \<= AF\_All \<= 0.001  
	IV:  0.001 \<= AF\_All \<= 0.01  
	V:   0.01 \<= AF\_All \<= 0.1  
	VI:  0.1 \<= AF\_All  
	 
7. remove all variants contain NA value for all analyzed scores (SIFT, Polyphen2\_HDIV, Polyphen2\_HVAR, MutationTaster, MutationAssessor, PROVEAN, VEST3, CADD, DANN, fathmm-MKL\_coding, MetaSVM, MetaLR, GERP++\_RS, phyloP7way\_vertebrate, phyloP20way\_mammalian, phastCons7way\_vertebrate, phastCons20way\_mammalian, SiPhy\_29way\_logOdds, 18 scores. Discard LRT and FATHMM for too many NAs).

### Input File
1. ExAC.all.hg19\_multianno.new.vcf.hq.1hom.input.narm.txt
2. ExAC.all.hg19\_multianno.new.vcf.hq.input.narm.txt
3. ExAC.nonpat.hg19\_multianno.new.vcf.hq.1hom.input.narm.txt
4. ExAC.nonpat.hg19\_multianno.new.vcf.hq.input.narm.txt
5. clinvar.hg19\_multianno.filt.new.narm.txt

Currently, use 3 as control dataset and 5 as pathogenic dataset

## Data Analysis
1. Read input file as Patho and Control
2. Analysis data in different levels  
	Score [score\_list] ×   
	Group [group list] ×   
	select [genome, gene\_specific\_select (gss)] ×   
	contaminate [True, False] ×   
	contaminate\_rate [contaminate\_start to contaminate\_end, step = contaminate\_step]  
	  
	- Group is the frequent group defined in data processing;  
	- Select: genome is genome-wide randomly selecting control data; gss is selecting same number of control mutation for each gene according to the number in path dataset.  
	- contaminate is always False for gss selection  
	- contaminate preform like: if contaminate\_rate = 0.05 and contaminate\_select = 9000, then randomly select 9000 patho mutations as patho, (9000×0.05 patho + 9000×0.95 ctrl) as ctrl.  
	  
	For each analysis level, calculate their PRC, AUPRC and other data.
3. Write all AUPRC, PR90, PR10 data into out\_file
4. Draw plots
5. Analysis out\_file  

