# MoncktonWaterfall


# Usage Instructions

Usage: MoncktonWaterfall.py [-h] -i,--inFastx INFASTX -o,--out OUT [-m,--motifs MOTIFS] [-L,--locus LOCUS]
                             [-q,--plotQV] [-Q,--quality MIN_QUALITY_THRESHOLD] [-t,--quality-type QUALITY_TYPE]
                             [-S,--sort-type SORT_TYPE] [-R,--reverse-complement FASTQ_RC]
                             [-a,--allele_length ALLELE_LENGTH] [-C,--colours COLOURS]
                             [-M,--max-read-length MAX_LENGTH] [-F,--flanks FLANKS] [-y,--ylabel YLABEL]
                             [-f,--format FORMAT] [-d,--dpi DPI]

Quick waterfall plot from sequencing reads in fasta/fastq format. Some options incompatible with fasta format.

optional arguments:
  -h, --help            show this help message and exit
  -i,--inFastx INFASTX  Input Fastx file. (required)
  -o,--out OUT          Output file. (required)
  -m,--motifs MOTIFS    Search motifs, comma separated, most frequent first, e.g. 'CGG,AGG'. Either motifs (-m) or
                        locus (-L) must be specified.
  -L,--locus LOCUS      By specifying the locus, the motifs and colours are defined. Additional motifs to be coloured
                        can be specified using the -m flag. Either motifs (-m) or locus (-L) must be specified.
                        (Options: HTT, DMPK1)
  -q,--plotQV           Plot additional QV waterfall (visualising base quality of reads). Default False
  -Q,--quality MIN_QUALITY_THRESHOLD
                        Set average read quality filter (Q score), see -t for other quality scoring options. Default 0
  -t,--quality-type QUALITY_TYPE
                        Set the read quality scoring method. Options: Average, First200. Default Average
  -S,--sort-type SORT_TYPE
                        Set how reads are sorted. Options: repeat_length, repeat_count, flank_position, flank_presence, length. Can combine options by comma-separating them, the order sets the hierarchy of which option is considered first (e.g. 'flank_presence,repeat_length').
                        Default Read length
  -R,--reverse-complement FASTQ_RC
                        Set whether reads should be reverse complemented based on sequence. Options: repeats, flanks,
                        both. Default False
  -a,--allele_length ALLELE_LENGTH
                        Set whether reads are split between two graphs based on allele. Set with the number of repeats
                        in the short allele. Not compatible with -q. Default None
  -C,--colours COLOURS  Specify colours for motifs using comma-separated HEX codes. Should match length and order of
                        motifs. Default colours from the 'tab10' colour palette.
  -M,--max-read-length MAX_LENGTH
                        Do not plot reads above this length. Default None
  -F,--flanks FLANKS    Comma separated upstream and downstream flanking motifs. Can be used to colour-code flanks,
                        remove reads without either or both flanks, or sort by downstream flank position. Default None
  -y,--ylabel YLABEL    Y-axis label for plot. Default Reads
  -f,--format FORMAT    Image format. Default png
  -d,--dpi DPI          Image resolution. Default 400

Minimum usage:
	python MoncktonWaterfall.py -i reads.fastq -m CAG,CCG,CAA,CCA,CCT -o reads_output.png
	
Specifying locus and plotting a read quality plot:
	python MoncktonWaterfall.py -i reads.fastq -L HTT -q -o reads_output.png
	
Specifying locus with flank as additional motifs:
	python MoncktonWaterfall.py -i reads.fastq -L HTT -m TTCCGATC -o reads_output.png
	
Reverse complement reverse reads and exclude reads longer than 1000 bases:
	python MoncktonWaterfall.py -i reads.fastq -L HTT -R -M 1000 -o reads_output.png

Filtering on read quality score (average > 30) and sorting by continuous CAG length:
	python MoncktonWaterfall.py -i reads.fastq -L HTT -Q 30 -t Average -S repeat_length -o reads_output.png

Specifying motifs with custom colours:
	python MoncktonWaterfall.py -i reads.fastq -m CAG,CCG,CAA,CCA,CCT -C #933b33,#9dbb61,#4bacc6,#ffc000,#ab9bc3 -o reads_output.png

Specifying locus and flanking sequences, and sorting by downstream flank position (works well as an alternative to sorting by repeat length):
	python MoncktonWaterfall.py -i reads.fastq -L HTT -F CGACCCT,AGCTTC -S flank_position -o reads_output.png
	

More details on options: 
-L/--locus -> you can simply choose the locus, which contains preset options for the motifs and colours; currently two options: 
	- HTT: using motifs: CAG,CAA,CCG,CCA,CCT and colours: #933b33,#9dbb61,#4bacc6,#ffc000,#ab9bc3, which translates to red, green, blue, yellow, purple
	- DMPK1: using motifs: CTG,CCG,CTC and colours: #933b33,#9dbb61,#4bacc6, translating to red, green and blue

-C/--colours ->  can be used to manually set colours for motifs specified under -m/—-motifs flag, the number of colours provided should match the number of motifs provided (if less colours are provided, will use default colours to colour any motifs provided). The colours should be provided in the same order as the motifs, to get the desired colour-coding. 
	- e.g. providing -m CAG,CAA,CCG -C #933b33,#9dbb61,#4bacc6, will colour CAG in #933b33 (red), CAA in #9dbb61 (green), and CCG in #4bacc6 (blue)
	- the -C flag currently only accepts HEX code (e.g. #933b33) and not names of colours (e.g. ‘red’), the HEX code for a wide range of colours can be found here: https://www.ditig.com/publications/256-colors-cheat-sheet or to find an appropriate palette and corresponding HEX codes you can use this interactive tool: https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=9

-S -> Can provide one or multiple options to be used to sort reads in plot. The hierarchy is set by the order, with the first option considered first. 
	- For example -S flank_presence,repeat_length will first order reads based on presence of both, one or no flanks in the sequence and consequently order by repeat length within those groups. 
	- Note that using flank_position or flank_presence requires the flanks to be specified using the -F flag.

-R -> Useful option for ONT reads. Based on the sequence determines which strand the read has come from and reverse complements as necessary 

-F -> These flanks are not the same as the flanking sequence in the amplicon. These should be limited to a short sequence upstream and downstream of the repeat tract



Additional changes from v4:
	- changed ordering in -S flank_presence
	- changed -R repeats to only count most frequent motif (the first one specified)
	- some error catching, defining choices for arguments etc.

Additional changes from v5:
	- when using -F and -R, the need to reverse complement will be determined based on the sequence between the flanks. No longer need to specify an option after -R 
