# DenPro
Fast software that calculates <b>Den</b>sity <b>Pro</b>file and precise mean density of aligned DNA sequence into inside and outside given regions.<br>
‘Density profile’ means a set of frequencies of the observed equal parts of the sequence with the same density. The program splits each given region over equal parts (called ‘windows’), and then counts the windows with the same density.<br>
‘Precise’ means that all the undefined regions in reference genome (i.e. regions filled by running ambiguous reference characters ‘N’) are excluded from consideration.<br>
If the input regions are not defined, the results will be calculated for the entire chromosome (each one in the sequence).

It runs on the command line under Windows, Linux and Mac OS X.

## Usage
```DenPro [options] -g|--gen <name> sequence```

## Help
```
Input:
  -g|--gen <name>       genome size file, or genome library, or single nucleotide sequence<br>
  --gap-len <int>       minimal length of undefined nucleotides region in genome
                        which is declared as a gap.
                        Ignored for genome size file [1000]
  -d|--dupl <OFF|ON>    accept duplicate reads [ON]
  --diff-sz <OFF|ON>    allow to ignore reads with different size [OFF]
Treatment:
  -c|--chr <chars>      treat stated chromosome only (all)
  --min-scr <int>       score threshold for treated reads (lack)
  --cons <int>          step of number of consolidated reads
  -f|--fbed <name>      'template' bed file which features define given regions
  -e|--exp-len <int>    length of expanding features in 'template' bed file [0]
  -s|--space <int>      resolution: span in bps in which reads will be counted
                        to define a density [50]
Ambig output:
  --alarm               output features ambiguities, if they exist
  --stat                output features ambiguities statistics, if they exist
Other:
  -t|--time             output run time
  -v|--version          print program's version and quit
  -h|--help             print usage information and quit
  ```

## Details

### Input
Aligned DNA sequence in BED format. Chromosomes in sequence can be unsorted, but reads within each chromosome should be sorted.<br>
Zipped files (.gz) are accepted too.

### Output
Mean density is measured in read per kilobase.<br>
Density profile is printed as a set of pairs \<number of read in window> – \<count of window>.<br>
The results are calculated for each chromosome separately, but the total mean density is printed as well.<br>
The results are printed on the screen and duplicated into a plain text file with a name ‘DenPro_out.txt’, located in current directory.<br>


### Options description
```-g|--gen <file>```<br>
Genome size file, or genome library, or single nucleotide sequence.<br>
Genome library is a directory contained nucleotide sequences for each chromosome in FASTA format.<br>
The difference between genome size file and genome library/file is that in the last case all the undefined regions in reference genome (gaps), will be excluded from calculation.<br>
Undefined regions are regions with only ambiguous reference characters ‘N’ in them.<br>
The minimal length of accounting gaps is managed by ```--gap-len``` option.<br>
For example, chromosome 1 from mm9 library contains 14 regions, separated by gaps with length more then 400 bps, and 10 regions, separated by gaps with length more then 1000.<br>
Indicating genome library has the same effect as ```-f|--fbed``` option, where ‘template’ is a set of defined regions.<br>
You can obtain genome library in  [UCSC](ftp://hgdownload.soe.ucsc.edu/goldenPath/) or in [Ensemble](ftp://ftp.ensembl.org/pub/release-73/fasta) storage. In the second case please copy genomic sequences with the same masked type only, f.e. unmasked (‘dna'), since program does not recognise mask’s types.<br>
Zipped .fa files can be mixed with unzipped.<br>
The single pointed FASTA file has the same effect as ```-c|--chr``` option.<br>
This option is required.

```--gap-len <int>```
Minimal length of undefined nucleotides region which is taken as a gap. For more details see ```-g|--gen``` option.<br>
Ignored for genome size file, pointed by ```-g|--gen```.<br>
Default: 1000

```-d|--dupl <OFF|ON>```
Accept or deny duplicated reads to participate in density calculation. Lower case value is appropriate as well.<br>
Default: ON

```--diff-sz <OFF|ON>```
Ignore reads with different length. Such reads are obtained for some alignments, especially in paired-end mode. They are scanty for Bowtie2, BWA, but can reach hundreds (MOSAIK) or even thousands (SMALT). This option allows to continue the calculation, otherwise the sequence is considered as incorrect. Issuance of information on such reads is regulated by options ```––alarm``` and ```––stat```.<br>
Default: OFF

```-c|--chr <chars>```
Treat stated chromosome only. Samples of option’s value: 1, 20, X.<br>
Reduces run time on 1.5-20 times depends of how far this chromosome is placed in an alignment. <br>
Default: all.

```--min-scr <int>```
Score threshold for treated reads. Reads with the score equal or less then stated will be ignored.<br>
Default: all reads are accepted.

```--cons <int>```
Step of number of consolidated windows. Allows to reduce output for very spread distribution by merging reads count for more than one ascending step. For example, setting the value of 5 will result in the issuance of pairs of \<1-5 reads> – \<N1 windows>, \<6-10 reads> – \<N2 windows> etc., in addition to ordinary output \<1 read> – \<X1 windows>, \<2 reads> – \<X2 windows> etc.<br>
Default: 1 (no consolidation).

```-f|--fbed <name>```
'Template' ordinary bed file with features that defines given regions.

```-e|--exp-len <int>```
Length of expanding regions in 'template' bed file.<br>
If set, all the features from 'template' bed file are be-expanding before processing: start positions are decreased for this value, end positions are increased.<br>
If expanded features become intersected, they are joined.<br>
This option is constructed for the special cases. For example, ChIP-seq treatment of transcription factor binding sites (TFBS). Such sites have a well-defined length (8-20 bps), while the length of corresponded enriched regions can reach 500-2000 pbs, but also uniform. If we know the TFBS coordinates, we can calculate read density within enriched regions immediately, by placing coordinates as 'template' bed file and setting this option.<br>
Default: 0.

```-s|--space <int>```
Resolution: the length of windows (span) in bs in which reads will be counted to define a density.<br>
Default: 50.

```--alarm```<br>
Output ambiguities, if they exist.
For alignments such ambiguities can be duplicated reads or reads with different length (size).<br>
In some circumstances you need to be aware about these issues. There are two ways to know about them: detailed and summary.
This option provides a detailed way. If it is pointed, information about type of ambiguity, number of line where it occurs, and resulting treatment would be printed each time when it happens.<br>
Summary way is managed by option ```–-stat```.<br>
Duplicated reads are not printed this way since they are considered as normal for alignment, but they are reported in summary way as well.

```--stat```<br>
Output features ambiguities statistics, if they exist. Prints number of all recognized certain type ambiguities, and appropriate treatment.<br>
In particular, this is a simple way to know the number of duplicated reads in alignments.<br>
For more details about ambiguities see ```--alarm``` option.

