##Manipulating VCF files

#1. Describe the format of the file and the data stored
VCF (Variant Calling Format) is a format that specifies a text file used in bioinformatics for storing gene sequence variations. 
Usually a DNA sample is sequenced through a next generation sequencing system,  producing a raw sequence file. 
That raw sequence data is then aligned, creating BAM/SAM files as a result. From there, variant calling identifies changes to a 
particular genome as compared to the reference genome. That output is stored in a variant call format, VCF for short. 
VCF is a text file format it is most likely stored in a compressed manner.It contains meta information lines, a header
line, and then data lines each containing information about a position in the genome. The format also has the ability
to contain genotype information on samples for each position.

#2. What does the header section of the file contain
The header contains an arbitrary number of meta-information lines, each starting with characters "##" and a TAB delimited field 
definition line, starting with a single "#"character. The meta-information header lines provide a standardized description of tags 
and annotations used in the data section. The use of meta-information allows the information stored within a VCF file to be tailored
to the dataset in question. It can be also used to provide information about the means of file creation, date of creation, version of 
the reference sequence, software used and any other information relevant to the history of the file. The field definition line names 
eight mandatory columns, corresponding to data columns representing the chromosome (CHROM), a 1-based position of the start of the variant 
(POS), unique identifiers of the variant (ID), the reference allele (REF), a comma separated list of alternate non-reference alleles (ALT),
a phred-scaled quality score (QUAL), site filtering information (FILTER) and a semicolon separated list of additional, user extensible annotation (INFO)

#3. How many samples are in the file
6 samples
bcftools query -l sample.vcf | wc -l
#Another option
vcf-query -l sample.vcf | wc -l

#4. How many variants are in the file
398246
#all these different command options give the same answer
bcftools view -H sample.vcf | wc -l
bcftools query -f '%ALT\n' sample.vcf | wc -l
grep -v '^#' sample.vcf | wc -l

#5. How would you extract the chromosome, position, QualByDepth and RMSMappingQuality fields? Save the output to a tab-delimited filePostion 1, 2, 

bcftools query -f '%CHROM\t%POS[\t%QD;%MQ]\n' sample.vcf > fields.txt

#6. Extract data that belongs to chromosomes 2,4 and MT
awk '$1=="2" || $1=="4" || $1=="MT"' sample.vcf > chromosome.txt

#7. Print out variants that do not belong to chr20:1-30000000

grep -v '^##' sample.vcf > chr20.txt

#8. Extract variants that belong to SRR13107019
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -s SRR13107019 sample.vcf > SRR7019.txt

#9. Filter out variants with a QualByDepth above 7
awk -F '\t' '{if ($6>=7) print $0}' sample.vcf > 9filtered.vcf

#10. How many contigs are referred to in the file. Check the header section
2211
grep -c "^##contig" sample.vcf 

#11. Comment on the eighth and ninth columns of the file

These  columns of a VCF (Variant Call Format) file typically contain information about the genotype of the samplebeing analyzed. 
The eighth column, also known as the "format" field, lists the data fields that are present in the ninth column, also known as the "sample" field.
These fields may include the genotype (GT), read depth (DP), and genotype quality (GQ) of the sample among others. The format of the fields in these 
columns can vary depending on the specific software used to generate the VCF file.The annotations contained in the INFO field are represented as 
tag-value pairs, where the tag and value are separated by an equal sign, ie =, and pairs are separated by colons.

#12. Extract data on the read depth of called variants for sample SRR13107018
bcftools querry -f '%DP\n' -s SRR13107018 sample.vcf > 12read_depth.vcf

#13. Extract data on the allele frequency of alternate alleles. Combine this data with the chromosome and position of the alternate allele
bcftools query -f '%CHROM\t%POS\t%AF\n' sample.vcf > alleles.vcf

#Manipulating SAM files
#1. Describe the format of the file and the data stored
The file begins with a header, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., 
this will change depending on the aligner being used. Following the header is the alignment section. Each line that follows corresponds to alignment 
information for a single read. Each alignment line has 11 mandatory fields for essential mapping information and a variable number of other fields for 
aligner specific information.

#2. What does the header section of the file contain

Each header line begins with the character ‘@’ followed by one of the two-letter header record type codes defined in this section.
In the header, each line is TAB-delimited and, apart from @CO lines, each data fieldfollows a format ‘TAG:VALUE’ where TAG is a two-character 
string that defines the format and content of VALUE. Thus header lines match /^@(HD|SQ|RG|PG)(\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ or /^@CO\t.*/.
Within each (non-@CO) header line, no field tag may appear more than once and the order in which the fields appear is not significant.

#3. How many samples are in the file
249
samtools view -H sample.sam | grep -c '@RG'

#4. How many alignments are in the file
35511
samtools view -c -F 4 sample.sam

#5. Get summary statistics for the alignments in the file

samtools flagstat sample.sam

  #summary statistics
36142 + 0 in total (QC-passed reads + QC-failed reads)
36142 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
35511 + 0 mapped (98.25% : N/A)
35511 + 0 primary mapped (98.25% : N/A)
27791 + 0 paired in sequencing
7407 + 0 read1
7409 + 0 read2
26134 + 0 properly paired (94.04% : N/A)
26529 + 0 with itself and mate mapped
631 + 0 singletons (2.27% : N/A)
228 + 0 with mate mapped to a different chr
126 + 0 with mate mapped to a different chr (mapQ>=5)

#6. Count the number of fields in the file

awk '{print NF}' sample.sam | sort -nu | tail -1
18 lines

#7. Print all lines in the file that have @SQ and sequence name tag beginning with NT_

grep '@SQ.*NT_' sample.sam > SQ.txt

#8. Print all lines in the file that have @RG and LB tag beginning with Solexa

grep '@RG.*LB:Solexa' sample.sam > solexa.txt

#9. Extract primarily aligned sequences and save them in another file.

samtools view -f 2 -b sample.sam > primary_alignment.bam

#10. Extract alignments that map to chromosomes 1 and 3. Save the output in BAM format.

samtools view -Sb sample.sam > sample.bam
samtools index sample.bam  #indexing the sam file first
samtools view -b sample.bam 1 3 > sample_chr1_3.bam

#11. How would you obtain unmapped reads from the file

samtools view -f 4 sample.bam> unmapped_reads_1.bam

#12. How many reads are aligned to chromosome 4
samtools index sample.bam
samtools view -f4 -c sample.bam
 631 reads
#13. Comment of the second and sixth column of the file
The second column of a SAM (Sequence Alignment/Map) file contains the name of the reference sequence (e.g. chromosome) that 
the read is aligned to. The sixth column, also known as the "CIGAR" (Compact Idiosyncratic Gapped Alignment Report) field, 
describes the alignment of the read to the reference sequence using a series of letters and numbers. The letters represent 
the type of alignment operation (e.g. M for match, I for insertion, and D for deletion) and the numbers represent the number of bases 
that the operation applies to. Together, this information allows one to infer the location, orientation, and quality of the alignment 
of the read to the reference.

#14. Extract all optional fields of the file and save them in “optional_fields.txt”

samtools view -h sample.bam | cut -f 12- | grep -Eo "RG:Z:[^\t]+" > Extracted_optional_fields.txt
