# DenPro
Fast software that calculates <b>Den</b>sity <b>Pro</b>file and precise mean density of aligned DNA sequence into inside and outside given regions.<br>
*Density profile* means a set of frequencies of the observed equal parts of the sequence with the same density. The program splits each given region over equal parts (called ‘windows’), and then counts the windows with the same density.<br>
*Precise* means that all the undefined regions in reference genome (i.e. regions filled by running ambiguous reference characters ‘N’) are excluded from consideration.<br>
If the input regions are not defined, the results will be calculated for the entire chromosome (each one in the sequence).

The program runs on the command line under Linux and Windows.

## Installation
### Ready executable file
Try this way first. 

**Linux**<br>
Go to the desire directory and type commands:<br>
```wget -O DenPro.gz https://github.com/fnaumenko/DenPro/releases/download/1.0/DenPro-Linux-x64.gz```<br>
```gzip -d DenPro.gz```<br>
```chmod +x DenPro```

**Windows**<br>
Download archive from [here](https://github.com/fnaumenko/DenPro/releases/download/1.0/DenPro-Windows-x64.zip) and unzip by any archiver, for instance **WinRar**.

### Compiling in Linux
Required libraries:<br>
g++<br>
zlib (optionally)

Go to the desired directory and type commands:<br>
```wget -O DenPro.zip https://github.com/fnaumenko/DenPro/archive/1.0.zip```<br>
```unzip DenPro.zip```<br>
```cd DenPro-1.0```<br>
```make```

If **zlib** is not installed on your system, a message will be displayed from the linker.<br>
In that case you can compile the program without the ability to work with .gz files. 
To do this, open *makefile* in any text editor, uncomment last macro in the second line, comment third line, save *makefile*, and try again ```make```.<br>
To be sure about **zlib** on your system, type ```whereis zlib```.

## Usage
```DenPro [options] -g|--gen <name> sequence```

## Help
```
Input:
  -g|--gen <name>       chromosome sizes file, reference genome library, or single nucleotide sequence. Required
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
Output:
  -i|--info <NOTE|STAT> output summary information about feature ambiguities, if they exist:
                        NOTE - notice, STAT - statistics
  -w|--warn             output each feature ambiguity, if they exist
  -o|--out              duplicate standard output to DenPro_out.txt file
Other:
  -t|--time             output run time
  -v|--version          print program's version and quit
  -h|--help             print usage information and quit
  ```

## Details

### Input
Aligned DNA sequence in [BED](https://www.ensembl.org/info/website/upload/bed.html) format. 
Chromosomes in sequence can be unsorted, but reads within each chromosome should be sorted.<br>

Compressed files in gzip format (.gz) are acceptable.

### Output
Mean density is measured in read per kilobase.<br>
Density profile is printed as a set of pairs \<number of read in window> – \<count of window>.<br>
The results are calculated for each chromosome separately, and the total mean density is printed.<br>

### Options
Enumerable option values are case insensitive.
Enumerable option values are case insensitive.

```-g|--gen <file>```<br>
Chromosome sizes file, reference genome library, or single nucleotide sequence.<br>
Genome library is a directory contained nucleotide sequences for each chromosome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
If ```name``` is a .fa[.gz] file, **DenPro** accepts the corresponding chromosome as the only treated.<br>
Otherwise first the program searches for .fa files in the directory ```name```. 
If there are no such files in this directory, **DenPro** searches for .fa.gz files.<br>
If chromosome is stated by option ```–c|--chr```, the program searches for the corresponding .fa[.gz] file.<br>
The difference between chromosome sizes file and genome library/file is that in the latter all the undefined regions in the reference genome (gaps), will be excluded from calculation.<br>
Undefined regions are regions with only ambiguous reference characters ‘N’ in them.<br>
The minimal length of accounting gaps is managed by ```--gap-len``` option.<br>
For example, chromosome 1 from mm9 library contains 14 regions, separated by gaps with length more then 400 bps, and 10 regions, separated by gaps with length more then 1000.<br>
Indicating genome library has the same effect as ```-f|--fbed``` option, where ‘template’ is a set of defined regions.<br>
The program searches for gaps once. On subsequent calls with the same length of gaps, it uses the search results stored in the specified directory.<br>
The program also generates once a chromosome sizes file in the same directory.<br>
One can obtain a genome library in  UCSC: ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble: ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, f.e. unmasked (‘dna'), since program does not recognise mask’s types.<br>
This option is required.

```--gap-len <int>```<br>
Minimal length of undefined nucleotides region which is taken as a gap. For more details see ```-g|--gen``` option.<br>
Ignored for chromosome sizes file, specified by ```-g|--gen```.<br>
Default: 1000

```-d|--dupl <OFF|ON>```<br>
Accept or deny duplicated reads to participate in density calculation.<br>
Default: ```ON```

```--diff-sz <OFF|ON>```<br>
Ignore reads with different length. 
Such reads are obtained for some alignments, especially in paired end mode. 
They are scanty for **Bowtie2**, **BWA**, but can reach hundreds (**MOSAIK**) or even thousands (**SMALT**). 
This option allows to continue the calculation, otherwise the sequence is considered as incorrect. 
Issuance of information on such reads is regulated by options ```––i|info``` and ```–w|-–warn```.<br>
Default: ```OFF```

```-c|--chr <chars>```<br>
Treat stated chromosome only. Samples of option’s value: 1, 20, X.<br>
Reduces run time on 1.5-20 times depending on how far this chromosome is placed in an alignment.<br>
Default: all.

```--min-scr <int>```<br>
Score threshold for treated reads. 
Reads with the score equal or less then stated will be ignored.<br>
Default: all reads are accepted.

```--cons <int>```<br>
Step of number of consolidated windows. 
Allows to reduce output for very spread distribution by merging reads count for more than one ascending step. 
For example, setting the value of 5 will result in the issuance of pairs of \<1-5 reads> – \<N1 windows>, \<6-10 reads> – \<N2 windows> etc., in addition to ordinary output \<1 read> – \<X1 windows>, \<2 reads> – \<X2 windows> etc.<br>
Default: 1 (no consolidation).

```-f|--fbed <name>```<br>
'Template' ordinary bed file with features that defines given regions.

```-e|--exp-len <int>```<br>
Length of expanding regions in 'template' bed file.<br>
If set, all the features from 'template' bed file are be-expanding before processing: *start* positions are decreased for this value, *end* positions are increased.<br>
If expanded features become intersected, they are joined.<br>
This option is constructed for the special cases, for instance, ChIP-seq treatment of transcription factor binding sites (TFBS). 
Such binding sites have a well-defined length (8-20 bps), while the length of corresponded enriched regions can reach 500-2000 pb, but also uniform. 
If we know the TFBS coordinates, we can calculate read density within enriched regions immediately, by placing coordinates as 'template' bed file and setting this option.<br>
This option is only relevant in addition to the option ```-f|--fbed```.<br>
Default: 0.

```-s|--space <int>```<br>
Resolution: the length of windows (span) in bp in which reads will be counted to define a density.<br>
Default: 100.

```-i|--info <NOTE|STAT>```<br>
Output feature ambiguity statistics, if they exist.<br> 
For alignments, such ambiguities can be duplicated reads or reads with different length (size).<br>
In some circumstances you need to be aware of these issues. 
There are two methods used to identify them: detailed and summary.<br>
This option provides the summary method. 
It forces to display number of all recognized certain type ambiguities, and appropriate treatment.<br>
If ```STAT``` value is set, the total number of reads left after treatment is displayed additionally.<br>
In particular, this is a simple way to know the number of duplicated reads.<br>
The detailed method is managed by option ```-w|--warn```.

```-w|--warn```<br>
Output ambiguity for each feature, if they exist.<br>
If it is specified, information about type of ambiguity, number of line where it occurs, and resulting treatment would be printed each time it happens.<br>
Duplicated reads are not printed when using this method as they are considered normal for alignment, but they are reported in the summary method, see ```-i|--info``` option.

```-o|--out```<br>
Duplicate standard output to **DenPro_out.txt** file. 
It is an analogue of the *tee* Linux command and is rather useful by calling **DenPro** under Windows.


