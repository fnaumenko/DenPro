# DenPro
Fast software that calculates <b>Den</b>sity <b>Pro</b>file and precise mean density of aligned DNA sequence into inside and outside given regions.<br>
*Density profile* means a set of frequencies of the observed equal parts of the sequence with the same density. 
The program splits each given region into equal non-overlapping parts (windows), 
and then counts the number of windows with the same density.<br>
*Precise* means that all the undefined regions in the reference genome (i.e. regions filled by running ambiguous reference characters ‘N’) are excluded from consideration.<br>
If the input regions are not defined, then only mean density is calculated for each chromosome.

The program runs on the command line under Linux and Windows.

## Installation
### Executable file

**Linux**<br>
Go to the desire directory and type commands:<br>
```wget -O DenPro.gz https://github.com/fnaumenko/DenPro/releases/download/1.0/DenPro-Linux-x64.gz```<br>
```gzip -d DenPro.gz```<br>
```chmod +x DenPro```

**Windows**<br>
Download archive from [here](https://github.com/fnaumenko/DenPro/releases/download/1.0/DenPro-Windows-x64.zip) 
and unzip by any archiver, for instance [WinRar](https://www.win-rar.com/download.html?&L=0).

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
  --gap-len <int>       minimal length of undefined nucleotides region in genome which is declared as a gap.
                        Ignored for the genome size file [1000]
  -d|--dupl <OFF|ON>    accept duplicate reads [ON]
Treatment:
  -c|--chr <chars>      treat specified chromosome only
  --min-scr <int>       score threshold for treated reads
  --cons <int>          step of number of consolidated reads
  -f|--fbed <name>      'template' bed file which features define treated regions
  -e|--ext-len <int>    length by which the features in the 'template' bed file
                        will be extended in both directions before treatment [0]
  -s|--space <int>      resolution: span in bp by which reads will be counted to define a density [100]
Output:
  -i|--info <NM|CNT|STAT>       print information about file:
                        NM - name only, CNT - number of reads, STAT - statistics [CNT]
  -w|--warn             print each read ambiguity, if they exist
  -o|--out              duplicate standard output to DenPro_out.txt file
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit 
```

## Details

### Input
Aligned DNA sequence in [BED](https://www.ensembl.org/info/website/upload/bed.html) format, in which reads should be clustered in chromosomes.<br> 
For faster processing, reads belonging to the same chromosome also should be sorted in position ascending order.<br>
The simplest way to prepare bed files is to sort them, for example by using **sortBed** utility from [bedtools](http://bedtools.readthedocs.io/en/latest/) package 
(though for **DenPro** the order of the chromosomes themselves does not matter).<br>

Compressed files in gzip format (.gz) are acceptable.

### Output
Mean density is measured in read per kilobase.<br>
Density profile is printed as a set of a ‘key’-‘value’ pairs: \<number of read in window>\<count of window>.<br>
The results are calculated for each chromosome separately, though the total mean density is also printed.

### Options description
Enumerable option values are case insensitive.

```-g|--gen <file>```<br>
Chromosome sizes file, reference genome library, or single nucleotide sequence.<br>
Genome library is a directory contained nucleotide sequences for each chromosome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
If ```name``` is a *.fa[.gz]* file, **DenPro** accepts the corresponding chromosome as the only treated.<br>
Otherwise first the program searches for *.fa* files in the directory ```name```. 
If there are no such files in this directory, **DenPro** searches for *.fa.gz* files.<br>
If chromosome is specified by option ```–c|--chr```, the program searches for the corresponding *.fa[.gz]* file.<br>
The chromosome sizes file is recognized by *.sizes* extension.<br>
The difference between chromosome sizes file and genome library/file is that in the latter all the undefined regions in the reference genome (gaps), will be excluded from calculation.<br>
Undefined regions are regions with only ambiguous reference characters ‘N’ in them.<br>
The minimal length of accounting gaps is managed by ```--gap-len``` option.<br>
For example, chromosome 1 from mm9 library contains 14 regions, separated by gaps with length more then 400 bps, and 10 regions, separated by gaps with length more then 1000.<br>
Indicating genome library has the same effect as ```-f|--fbed``` option, where 'template' is a set of defined regions.<br>
The program searches for gaps once. On subsequent calls with the same length of gaps, it uses the search results stored in the specified directory.<br>
The program also generates once a chromosome sizes file in the same directory.<br>
One can obtain a genome library in  UCSC: ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble: ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, f.e. unmasked (‘dna'), since program does not recognise mask’s types.<br>
This option is required.

```--gap-len <int>```<br>
Minimal length of undefined nucleotides region which is taken as a gap. 
For more details see ```-g|--gen``` option.<br>
Ignored for chromosome sizes file, specified by ```-g|--gen```.<br>
Default: 1000

```-d|--dupl <OFF|ON>```<br>
Accept or deny duplicated reads to participate in density calculation.<br>
Default: ```ON```

```-c|--chr <chars>```<br>
Treat specified chromosome only. Samples of option’s value: 1, 20, X.<br>
Reduces run time on 1.5-20 times depending on how far this chromosome is placed in an alignment.<br>

```--min-scr <int>```<br>
Score threshold for treated reads. 
Reads with the score equal or less then stated will be ignored.<br>
Default: all reads are accepted.

```--cons <int>```<br>
Step of number of consolidated windows. 
Allows to reduce output for very spread distribution by merging reads count for more than one ascending step. 
For example, setting the value of 5 will result in the issuance of pairs of *\<1-5 reads>\<N1 windows>*, *\<6-10 reads>\<N2 windows>* etc., in addition to ordinary output *\<1 read>\<X1 windows>*, *\<2 reads>\<X2 windows>* etc.<br>
Default: 1 (no consolidation).

```-f|--fbed <name>```<br>
'Template' ordinary bed file with features that defines treated regions. 
The density profile will be constructed on the region resulting from the merger of all regions specified by the 'template' regions.<br>
This option abolishes the merge of defined regions specified by ```-g|--gen <name>``` option.

```-e|--ext-len <int>```<br>
Value by which all features in 'template' bed should be stretched in both directions before the density count.<br>
If set, all the features from 'template' bed file will be stretched before processing: *start* positions are decreased for this value, *end* positions are increased.<br>
If stretched features become intersected, they are joined.<br>
This option is constructed for the special cases, for instance, ChIP-seq treatment of transcription factor binding sites (TFBS). 
Such binding sites have a well-defined uniform length (8-20 bp), while the length of corresponded enriched regions can reach 500-2000 pb. 
If we know the TFBS coordinates, we can calculate read density within enriched regions immediately, by placing coordinates as 'template' bed file and setting this option.<br>
This option is only relevant in addition to the option ```-f|--fbed```.<br>
Default: 0.

```-s|--space <int>```<br>
Resolution: the length of windows (span) in bp by which reads will be counted to define a density.<br>
Read is considered belonging to span if it`s centre is placed within the span. 
Thus each read is counted once. 
Then, while using the reference genome from the input sequences, and if 'template' is not specified, the reference undefined regions (gaps) are preliminarily removed from the compared sequences. 
This means that the regions separated by a gap are merged.<br>
As a result, the program compares the actual read density distributions.<br>
Default: 100.

```-i|--info <NM|CNT|STAT>```<br>
Output information about number of items (features/reads/intervals).<br>
```NM```:  brief output. Prints file names without any additional information.<br>
```CNT```:  prints file names and number of all and accepted items, if they are different.<br>
```STAT```: prints item ambiguities statistics, if they exist.<br>
For the alignments, such ambiguities can be duplicated reads or reads with different length (size).<br>
Thus, not all reads present in the file can be accepted.<br>
In some circumstances you need to be aware of these issues. 
There are two methods used to identify them: detailed and summary.<br>
The ```STAT``` value provides the summary method. 
It forces to display number of all recognized certain type ambiguities, and appropriate treatment.<br>
In particular, this is a simple way to know the number of duplicated reads.<br>
The detailed method is managed by option ```-w|--warn```.

```-w|--warn```<br>
Output ambiguity for each read, if they exist.<br>
If it is specified, information about type of ambiguity, number of line where it occurs, and resulting treatment would be printed each time it happens.<br>
Duplicated reads are not printed when using this method as they are considered normal for the alignment, 
but they are reported in the summary method, see ```-i|--info``` option.

```-o|--out```<br>
Duplicate standard output to **DenPro_out.txt** file. 
It is an analogue of the *tee* Linux command and is rather useful by calling **DenPro** under Windows.

##
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write me on fedor.naumenko@gmail.com
