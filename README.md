# minibar
Dual barcode and primer demultiplexing for MinION sequenced reads

Minibar is developed for and accompanies the paper:
\
Long-read nanopore sequencing of ribosomal DNA: a portable, cost-effective, phylogenetically robust approach for biodiversity assessments across broad taxonomic scale. Authors: Henrik Krehenwinkel, Aaron Pomerantz, James B. Henderson, Susan R. Kennedy, Jun Ying Lim, Varun Swamy, Juan Diego Shoobridge, Nipam H. Patel, Rosemary G. Gillespie, Stefan Prost

    Usage: minibar.py barcode_file sequence_file [-pct <pct> | -e <int> -E <int>] [-l <int>]
                                                 [-F [-P <prefix>]] [-M 1|2|3]
                                                 [-S | -T | -C | -CC | -D]
                                                 [-cols <int_list>] [-info fwd|rev|primer]
                                                 [-w] [-fh | -nh] [-n <num_seqs> | -n <first_seq>,<num_seqs>]

        Identify MinION sequence by dual barcode indexes and primers.
        The sequence file can be in Fasta or Fastq format, gzipped or plain text.
        Sample ID is placed at end of header comment with match hit info before it.
        (minibar.py version 0.18)

        Example: ./minibar.py -C -F Demultiplex.txt example.fq

        -h display this with all option's descriptions     -v displays version
        -p <pct> percentage match (.75)
        -e <int> barcode edit distance value, overrides -p (4)
        -E <int> primer edit distance value (11)
        -l <int> length to search for index and primer at start and end of sequence (80)

        -F create individual sample files for sequences with -S or -C output (default: False)
        -f write to stdout instead of creating files (default: True)
        -P <str> if -F, <str> is prefix for individual files, followed by sample ID. (default: sample_)

        -M 1|2|3 Method to identify sample types using the barcodes (default: 3)
                 1 requires approximate match of barcode and primer at sequence start, this
                   and barcodes matched at the other end are used to identify sample IDs
                 2 finds matched barcodes on both ends of sequence, identifies pairs that match a sample ID
                 3 uses Method 1 and if it does not succeed, uses Method 2

        -S outputs sequence record in fasta or fastq format of input (default output)
        -T trims barcode and primer from each end of the sequence, then outputs record
        -C similar to S but uses upper/lower case to show found barcode indexes and primers
        -CC also colors found barcode blue, primer green if found, primer red otherwise
        
        -cols <int_list> column position in barcode_file for: sample, fwd index, fwd primer, rev index, rev primer
                 (default: 1,2,3,4,5 if 5 cols; 1,3,4,5,6 if 6 cols; 1,3,4,6,7 if 7 cols; 1,3,5,8,10 if 10 or more cols)
                        
        -info fwd|rev|primer display barcode index or primer info, including edit distances
        -w  treat duplicates in barcode_file as warning, not error
        -fh first line of barcode file considered a header (default: auto detect header)
        -nh first line of barcode file is not a header (default: auto detect header)

        -n <num_seqs> number of sequences to read from file (ex: -n 100)
        -n <first_seq>,<num_seqs> (ex: -n 102,3)
                 
### Requirements
**minibar.py** is written in Python version 2.7 and is compatible with Python version 3. It imports the edlib library which
you can install using **pip install edlib** and, if you are interested, its source can be found at https://github.com/Martinsos/edlib. With the requirements of a typical Python installation, the minibar.py source and the edlib module installed, minibar should run on MacOS, Linux and Windows.

### Examples
This simply adds the Sample ID and hit quality info to the comment field of each record and saves the results in another Fasta file. Also shown are the headers for the first two Fasta records with the added information in **bold** which you'll need to scroll to the end of the lines to see.
<pre>
$ minibar.py IndexCombinationPeperomonia.txt PeperomiaTestSet.fasta >PeperomiaTestSet_SampleIDs.fa
IndexCombinationPeperomonia.txt PeperomiaTestSet.fasta Index edit dist 4, Primer edit dist 11, Search Len 80, Search Method 3, Output Type S
750 seqs: H 750 HH 679 Hh 62 hh 0 IDs 741 Mult_IDs 0 (0.1245s)

$ grep "^>" PeperomiaTestSet_SampleIDs.fa -m 2
>5fcfed05-3207-4f44-bbc8-dd8d62042384 runid=8eea9dccb71575d312e59de22190819144152897 read=24 ch=325 start_time=2018-04-02T23:38:51Z <b>H-(0,6),H+(0,2) Jun_40</b>
>8d745d8e-3d06-4778-81cc-d25dd7f3b3a7 runid=8eea9dccb71575d312e59de22190819144152897 read=25 ch=262 start_time=2018-04-02T23:39:03Z <b>H+(1,1),H-(1,1) Jun_38</b>
</pre>
