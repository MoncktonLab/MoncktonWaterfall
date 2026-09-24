import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
from collections import OrderedDict
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm
import sys,re,pysam
from operator import attrgetter

COLORMAP = list(matplotlib.colormaps['tab10'].colors)
UNKNOWN  = (0.75,)*3
BLANK    = (1,)*3
FLANK_COLOR = (0,)*3
XLABEL   = 'Position'
YLABEL   = 'Reads'
MAXQV    = 40
MINQV    = 0
CENTERQV = 20
QVCOLOR  = matplotlib.colormaps['RdYlBu']


def main(parser):
    global COLORMAP # setting COLORMAP as global because need to modify it
    args = parser.parse_args()
    
    min_quality_threshold = args.min_quality_threshold
    quality_type = args.quality_type
    sort_type = args.sort_type
    # Keep the legacy return value defined even when --split-direction is used.
    # In split mode the individual forward/reverse rasters are local to the
    # plotting branch, so there is no single `raster` object to return.
    raster = None
    
    infastx = args.inFastx
    #motifs = args.motifs.split(',') # change from v4 (to merge with v5)
    fastq_RC = args.fastq_RC
    if args.split_direction and not fastq_RC:
        raise Waterfall_Exception('--split-direction requires -R/--reverse-complement so that read orientation can be identified.')
    max_length = args.max_length
    min_length = args.min_length
 
    if args.allele_length != None:
      if '+' in args.allele_length:
          allele_length = int(args.allele_length.split('+')[0])
          difference = int(args.allele_length.split('+')[1])
      else:
          allele_length = int(args.allele_length)
          difference = 5

    # change from v4 (to merge with v5)
    if args.motifs == None and args.locus == None:
        raise Waterfall_Exception('Must specify locus (-L) or motifs (-m).')
    if args.motifs:
        motifs = args.motifs.split(',')
    else:
        motifs = None
    if args.colours:
        colours = args.colours.split(',')
    else:
        colours = None

    if args.flanks:
        flanks = args.flanks.split(',')
        upflank = flanks[0]
        downflank = flanks[1]

    # if locus specified, set motifs_locus and merge with any additional motifs specified by user
    if args.locus:
        if args.locus == 'HTT':
            motifs_locus = ['CAG','CCG','CAA','CCA','CCT']
            colours_locus = ['#DB6968','#4D97CD','#46C1BE','#E8C559','#B279B4']
        elif args.locus == 'DMPK1':
            motifs_locus = ['CTG','CCG','CTC']
            colours_locus = ['#DB6968', '#4D97CD','#46C1BE']
        elif args.locus == 'ERDA1':
            motifs_locus = ['CAG','CCG','CAC','CAT']
            colours_locus = ['#DB6968', '#4D97CD','#46C1BE','#E8C559']
        COLORMAP_locus = [mcolors.to_rgb(hex_color) for hex_color in colours_locus]
        if args.motifs:
            motifs_all = list(OrderedDict.fromkeys(motifs_locus+motifs))
            if colours:
                if len(motifs) <= len(colours):
                    if all(colour.startswith('#') for colour in colours):
                            colors_rgb = [mcolors.to_rgb(hex_color) for hex_color in colours]
                            COLORMAP_final = COLORMAP_locus+colors_rgb
                    else:
                        print('Colours must be in HEX code (e.g. #933b33). Using predefined colours.')
                        COLORMAP_final = COLORMAP_locus+COLORMAP
                else:
                    print('More motifs than colours specified. Using predefined colours.')
                    COLORMAP_final = COLORMAP_locus+COLORMAP
            else:
                COLORMAP_final = COLORMAP_locus+COLORMAP
        else:
            COLORMAP_final = COLORMAP_locus
            motifs_all = motifs_locus
    else:
        motifs_locus = None
        motifs_all = motifs
        if args.colours:
            if len(motifs) <= len(colours):
                if all(colour.startswith('#') for colour in colours):
                    colors_rgb = [mcolors.to_rgb(hex_color) for hex_color in colours]
                    COLORMAP_final = colors_rgb
                else:
                    print('Colours must be in HEX code (e.g. #933b33). Using predefined colours.')
                    COLORMAP_final = COLORMAP
            else:
                print('More motifs than colours specified. Using predefined colours.')
                COLORMAP_final = COLORMAP
        else:
            COLORMAP_final = COLORMAP

    # append the motifs with the code for the flanks to be coloured 
    if args.flanks:
        motifs_all.append('U')
        motifs_all.append('D')
        COLORMAP_final.append(FLANK_COLOR)
        COLORMAP_final.append(FLANK_COLOR)

    if fastq_RC:
        # Handle fastq_RC and generate output filename
        infastx_file = infastx.split('/')[-1]
        infastx_path = infastx.replace(infastx_file, '')
        outfastx = f"{infastx_path}out_{infastx_file}" if not infastx.startswith("out_") else f"out_{infastx}"
        if args.flanks:
            assignFlanks(infastx, upflank, downflank)
            fastq_to_process = "temp.fastq"
            reverseComplementFastq(fastq_to_process, outfastx, motifs_all)
            sortedRecs = sorted(list(pysam.FastxFile(outfastx, 'r')), key=sortFunc(args.sort_type, motifs_all[0], args.locus))
        else:
            reverseComplementFastq(infastx, outfastx, motifs_all)
            sortedRecs = sorted(list(pysam.FastxFile(outfastx, 'r')), key=sortFunc(args.sort_type, motifs_all[0], args.locus))
    else:
        if args.flanks:
            assignFlanks(infastx, upflank, downflank)
            fastq_to_process = "temp.fastq"
            sortedRecs = sorted(list(pysam.FastxFile(fastq_to_process, 'r')), key=sortFunc(args.sort_type, motifs_all[0], args.locus))
        else:
            sortedRecs = sorted(list(pysam.FastxFile(infastx, 'r')), key=sortFunc(args.sort_type, motifs_all[0], args.locus))

    if not sortedRecs:
        print(f'No records in {infastx}')
        return None

    if max_length is not None:
        sortedRecs = [rec for rec in sortedRecs if len(rec.sequence) <= max_length]
    if min_length is not None:
        sortedRecs = [rec for rec in sortedRecs if len(rec.sequence) >= min_length]
		
    if quality_type is not None:
        filteredRecs = filterByQuality(sortedRecs, min_quality_threshold, quality_type)
    else:
        filteredRecs = sortedRecs



    # change from v3
    if args.allele_length:
        short_reads, long_reads = split_reads_by_allele_length(filteredRecs, allele_length, difference, motifs_all[0], args.locus)
        colors  = OrderedDict([(m,COLORMAP_final[i]) for i,m in enumerate(motifs_all)])
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        raster = motifRaster(short_reads, motifs_all, colors)
        # replace 'U' and 'D' keys to 'flank' in legend
        legend_list = []
        for motif,color in list(colors.items()):
            if motif not in ['U','D']:
                legend_list.append((motif,color))
            elif motif == 'U':
                motif = 'flank'
                legend_list.append((motif,color))
        patches = [ mpatches.Patch(color=color, label=motif ) for motif,color in legend_list]
        n_rows = raster.shape[0]
        plotWaterfall(raster, XLABEL, args.ylabel, labels=patches, ax=ax1)
        ax1.set_title(f'Short Allele Reads, n= {n_rows}')
        raster = motifRaster(long_reads, motifs_all, colors)
        n_rows = raster.shape[0]
        plotWaterfall(raster, XLABEL, args.ylabel, labels=patches, ax=ax2)
        ax2.set_title(f'Long Allele Reads, n= {n_rows}')
        out = f"{args.out}.{args.format}"
        plt.tight_layout()
        fig.savefig(out, dpi=args.dpi, format=args.format)
        if args.plotQV:
            # Generate the QV raster and plot it, now with the filtered records
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
            qvraster = qvRaster(short_reads, 'Base QV')
            norm     = TwoSlopeNorm(CENTERQV, vmin=MINQV, vmax=MAXQV)
            plotWaterfall(qvraster, XLABEL, args.ylabel, norm=norm, cmap=QVCOLOR, colorbar='QV', ax=ax1)
            ax1.set_title('Short Allele Reads - Quality')
            qvraster = qvRaster(long_reads, 'Base QV')
            norm     = TwoSlopeNorm(CENTERQV, vmin=MINQV, vmax=MAXQV)
            plotWaterfall(qvraster, XLABEL, args.ylabel, norm=norm, cmap=QVCOLOR, colorbar='QV', ax=ax2)
            ax2.set_title('Long Allele Reads - Quality')
            name, ext = out.rsplit('.', 1)
            plt.tight_layout()
            fig.savefig(f'{name}.QV.{ext}', dpi=args.dpi, format=args.format)  
    else:
        colors  = OrderedDict([(m,COLORMAP_final[i]) for i,m in enumerate(motifs_all)])

        # replace 'U' and 'D' keys to 'flank' in legend
        legend_list = []
        for motif,color in list(colors.items()):
            if motif not in ['U','D']:
                legend_list.append((motif,color))
            elif motif == 'U':
                motif = 'flank'
                legend_list.append((motif,color))
        patches = [ mpatches.Patch(color=color, label=motif ) for motif,color in legend_list]

        out = args.out if args.out.endswith(args.format) else '%s.%s' % (args.out,args.format)

        if args.split_direction:
            # The orientation has to have been identified before plotting.
            # With -R, reverseComplementFastq records the original orientation
            # in the FASTQ header, even though reverse reads may have been
            # reverse-complemented for plotting.
            forward_reads = [rec for rec in filteredRecs if getReadOrientation(rec) == "forward"]
            reverse_reads = [rec for rec in filteredRecs if getReadOrientation(rec) == "reverse"]

            f, ax = plotSplitWaterfall(
                forward_reads,
                reverse_reads,
                XLABEL,
                args.ylabel,
                motifs_all,
                colors,
                labels=patches
            )
            ax.set_title(f'{out}, forward n={len(forward_reads)}, reverse n={len(reverse_reads)}')
            plt.tight_layout()
            f.savefig(out, dpi=args.dpi, format=args.format)

            if args.plotQV:
                # Plot QV with the same orientation and row ordering.
                forward_sorted = sorted(forward_reads, key=lambda rec: len(rec.sequence))
                reverse_sorted = sorted(reverse_reads, key=lambda rec: len(rec.sequence), reverse=True)
                all_reads = forward_sorted + reverse_sorted
                width = max(len(rec.sequence) for rec in all_reads) if all_reads else 1

                def padded_qv_raster(recs):
                    if not recs:
                        return np.ma.masked_all((0, width))
                    raster = qvRaster(recs, 'Base QV')
                    if raster.shape[1] < width:
                        pad = np.ma.masked_all((raster.shape[0], width - raster.shape[1]))
                        raster = np.ma.concatenate((raster, pad), axis=1)
                    elif raster.shape[1] > width:
                        raster = raster[:, :width]
                    return raster

                qv_forward = padded_qv_raster(forward_sorted)
                qv_reverse = padded_qv_raster(reverse_sorted)
                qv_reverse_for_plot = qv_reverse[::-1]

                fig, qvax = plt.subplots()
                if len(forward_sorted):
                    norm = TwoSlopeNorm(CENTERQV, vmin=MINQV, vmax=MAXQV)
                    qvax.imshow(qv_forward, origin='lower', aspect='auto',
                                interpolation='nearest',
                                extent=(0, width, 0, len(forward_sorted)),
                                norm=norm, cmap=QVCOLOR)
                if len(reverse_sorted):
                    norm = TwoSlopeNorm(CENTERQV, vmin=MINQV, vmax=MAXQV)
                    qvax.imshow(qv_reverse_for_plot, origin='lower', aspect='auto',
                                interpolation='nearest',
                                extent=(0, width, -len(reverse_sorted), 0),
                                norm=norm, cmap=QVCOLOR)
                qvax.axhline(0, color='black', linewidth=0.8)
                qvax.set_xlabel(XLABEL)
                qvax.set_ylabel(args.ylabel)
                qvax.set_title(f'{out}, forward n={len(forward_sorted)}, reverse n={len(reverse_sorted)} - Quality')
                fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=QVCOLOR), ax=qvax, shrink=0.6, label='QV')
                plt.tight_layout()
                name, ext = out.rsplit('.', 1)
                fig.savefig(f'{name}.QV.{ext}', dpi=args.dpi, format=args.format)
        else:
            raster  = motifRaster(filteredRecs,motifs_all,colors)
            n_rows = raster.shape[0]
            f,ax    = plotWaterfall(raster,XLABEL,args.ylabel,labels=patches)
            ax.set_title(f'{out}, n= {n_rows}')
            plt.tight_layout()
            f.savefig(out,dpi=args.dpi,format=args.format)
            if args.plotQV:
                # Generate the QV raster and plot it, now with the filtered records
                qvraster = qvRaster(filteredRecs, 'Base QV')
                norm     = TwoSlopeNorm(CENTERQV, vmin=MINQV, vmax=MAXQV)
                f, ax    = plotWaterfall(qvraster, XLABEL, args.ylabel, norm=norm, cmap=QVCOLOR, colorbar='QV')
                name, ext = out.rsplit('.', 1)
                f.savefig(f'{name}.QV.{ext}', format=args.format) 
                
    print("Done")
    return raster



def assignFlanks(infastx, upflank, downflank):
    with pysam.FastxFile(infastx) as input_fastq, open("temp.fastq", "w") as output_fastq:
        up_rc = reverseComplementSequence(upflank)
        down_rc = reverseComplementSequence(downflank)
        up_replacement = "U" * len(upflank)
        down_replacement = "D" * len(downflank)

        for rec in input_fastq:
            sequence = rec.sequence
            quality = rec.quality
            up_pos = sequence.find(upflank)
            down_pos = sequence.find(downflank)

            if up_pos != -1 and down_pos != -1:
                # Forward flanks
                modified_sequence = sequence.replace(upflank, up_replacement).replace(downflank, down_replacement)
                flank_info = f"upflank_{up_pos}:downflank_{down_pos}"
            else:
                # Try reverse complement flanks
                up_rc_pos = sequence.find(up_rc)
                down_rc_pos = sequence.find(down_rc)
                if up_rc_pos != -1 and down_rc_pos != -1:
                    modified_sequence = sequence.replace(up_rc, up_replacement).replace(down_rc, down_replacement)
                    flank_info = f"upflank_{up_rc_pos}_rc:downflank_{down_rc_pos}_rc"
                else:
                    # No matching pair found
                    modified_sequence = sequence
                    flank_info = "upflank_-1:downflank_-1"
            print(flank_info)
            modified_header = f"{rec.name}:{flank_info}"
            output_fastq.write(f"@{modified_header}\n")
            output_fastq.write(f"{modified_sequence}\n")
            output_fastq.write(f"+\n")
            output_fastq.write(f"{quality}\n")



def countCAG(rec, motif):
    return rec.count(motif)

def getFlankPos(rec):
    name_split = rec.split(":")
    for part in name_split:
        if "downflank" in part:
            downflank_pos = part.split("_")[1] 
            return int(downflank_pos)
    return -1 

def getFlankPresence(rec):
    downflank = None
    upflank = None
    name_split = rec.split(":")
    for part in name_split:
        if "downflank" in part:
            downflank = part
        if "upflank" in part:
            upflank = part

    downflank_pos = int(downflank.split("_")[1]) if downflank else -1
    upflank_pos = int(upflank.split("_")[1]) if upflank else -1
    if downflank_pos == -1 and upflank_pos == -1:
        presence = 0  
    elif downflank_pos != -1 and upflank_pos != -1:
        presence = 2  
    elif downflank_pos == -1 and upflank_pos != -1:
        presence = 1  
    elif downflank_pos != -1 and upflank_pos == -1:
        presence = 3 
    return presence
    
                              
def countCountinuousCAG(rec, motif, locus):
    count = 0
    length = len(motif)
    i = 0
    while i < len(rec) - length + 1:
        # Case 1: Continuous "CAG" motifs
        if rec[i:i+length] == motif and rec[i+length:i+2*length] == motif:
                count += 1
                i += length  # Move forward by one "CAG" to continue checking
        # Case 2: Single non-"CAG" sandwiched between "CAG"s
        elif rec[i:i+length] == motif and rec[i+2*length:i+3*length] == motif:
                count += 1
                i += length  # Skip past the non-"CAG" triplet
                continue
        else:
            i += 1
    #if locus == "HTT" and "CAACAG" in rec:
    #    count -= 2 
    return count


def checkReverseComplement(header, sequence, motifs):
    up_pos, down_pos = -1, -1
    is_rc = False
    name_split = header.split(":")
    for part in name_split:
        if "upflank" in part:
            up_pos = int(part.split("_")[1])
            if "rc" in part:
                is_rc = True
        elif "downflank" in part:
            down_pos = int(part.split("_")[1])
            if "rc" in part:
                is_rc = True
        
    # Determine flanked region
    if up_pos != -1 and down_pos != -1:
        region = sequence[min(up_pos, down_pos):max(up_pos, down_pos)]
    elif up_pos != -1:
        region = sequence[:up_pos] if is_rc else sequence[up_pos:]
    elif down_pos != -1:
        region = sequence[down_pos:] if is_rc else sequence[:down_pos]
    else:
        region = sequence

    # Count motifs
    for_count = sum(region.count(motif) for motif in motifs)
    rev_count = sum(region.count(reverseComplementSequence(motif)) for motif in motifs)

    if rev_count > for_count:
        return reverseComplementSequence(sequence)
    return sequence



def reverseComplementSequence(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U':'U', 'D':'D'}
    sequence = ''.join(complement[base] for base in reversed(sequence))
    return sequence
	
def reverseComplementFastq(input_file, output_file, motifs):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = []
        for i, line in enumerate(infile):
            lines.append(line.strip())
            if (i + 1) % 4 == 0:
                header, seq, plus, qual = lines
                new_seq = checkReverseComplement(header, seq, motifs)

                # Keep track of the original read orientation.  This is needed
                # later when --split-direction is requested.  The sequence can
                # still be reverse-complemented for plotting, but the original
                # orientation is retained in the header.
                if new_seq != seq:
                    orientation = "reverse"
                    # Reverse quality score if sequence was reversed
                    qual = qual[::-1]
                else:
                    orientation = "forward"

                header = f"{header}:orientation_{orientation}"
                outfile.write(f"{header}\n{new_seq}\n{plus}\n{qual}\n")
                lines = []


def getReadOrientation(rec):
    """Return the original read orientation recorded in the FASTQ header."""
    if ":orientation_reverse" in rec.name:
        return "reverse"
    if ":orientation_forward" in rec.name:
        return "forward"

    # When flanks were assigned, their *_rc annotation already identifies
    # the original orientation.
    if "_rc" in rec.name:
        return "reverse"

    return "forward"
                
def filterByQuality(recs, min_quality_threshold, quality_type):
    filtered_recs = []
    if quality_type == "Average":
        for rec in recs:
            qv = phred2QV(rec)  # Get the quality values for this record
            if np.mean(qv) >= min_quality_threshold:  # Check if all quality values are above the threshold
                filtered_recs.append(rec)
    if quality_type == "First200": 
        for rec in recs:
            qv = phred2QV(rec)  # Get the quality values for this record
            if np.all(qv[:200] >= min_quality_threshold):  # Check if all quality values are above the threshold
                filtered_recs.append(rec)
    return filtered_recs
    
# change from v3    
def split_reads_by_allele_length(recs, allele_length, difference, motif, locus):
    short_reads = [rec for rec in recs if countCountinuousCAG(rec.sequence, motif, locus) < allele_length + difference]
    long_reads = [rec for rec in recs if countCountinuousCAG(rec.sequence, motif, locus) >= allele_length + difference]
    return short_reads, long_reads

def sortFunc(sort_type, motif, locus):
    if sort_type == None:
        return lambda rec: -len(rec.sequence)
    elif sort_type == "repeat_length":
        return lambda rec: -countCountinuousCAG(rec.sequence, motif, locus)
    elif sort_type == "repeat_count":
        return lambda rec: -countCAG(rec.sequence, motif)
    elif sort_type =="flank_position":
        return lambda rec: -getFlankPos(rec.name)
    elif sort_type =="flank_presence":
        return lambda rec: getFlankPresence(rec.name)


def plotSplitWaterfall(forward_reads, reverse_reads, xlabel, ylabel, motifs, colors, labels=None, ax=None, **kwargs):
    """
    Plot forward reads above the x-axis and reverse reads below it.

    Forward reads are ordered shortest -> longest, so the shortest read is
    closest to the x-axis. Reverse reads are ordered longest -> shortest,
    so the longest reverse read is closest to the x-axis.
    """
    if ax is None:
        f, ax = plt.subplots()
    else:
        f = plt.gcf()

    # Keep a common x dimension. motifRaster already uses the first sequence
    # length as the raster width, so use the maximum read length here and pad
    # each half before plotting.
    all_reads = forward_reads + reverse_reads
    if not all_reads:
        raise ValueError("No reads available for split-direction waterfall")

    width = max(len(rec.sequence) for rec in all_reads)

    def padded_motif_raster(recs):
        if not recs:
            return np.empty((0, width, 3))
        raster = motifRaster(recs, motifs, colors)
        if raster.shape[1] < width:
            pad = np.ones((raster.shape[0], width - raster.shape[1], 3))
            raster = np.concatenate((raster, pad), axis=1)
        elif raster.shape[1] > width:
            raster = raster[:, :width, :]
        return raster

    # Forward: shortest -> longest, with shortest adjacent to y=0.
    forward_sorted = sorted(forward_reads, key=lambda rec: len(rec.sequence), reverse=True)

    # Reverse: longest -> shortest, with longest adjacent to y=0.
    reverse_sorted = sorted(reverse_reads, key=lambda rec: len(rec.sequence))

    forward_raster = padded_motif_raster(forward_sorted)
    reverse_raster = padded_motif_raster(reverse_sorted)

    # For the lower half, the last raster row is closest to y=0. Reverse the
    # raster so the longest reverse read (first in reverse_sorted) is nearest
    # the x-axis and the shortest is furthest below it.
    #reverse_raster_for_plot = reverse_raster[::-1]

    if len(forward_sorted):
        ax.imshow(
            forward_raster,
            origin='lower',
            aspect='auto',
            interpolation='nearest',
            extent=(0, width, 0, len(forward_sorted)),
            **kwargs
        )

    if len(reverse_sorted):
        ax.imshow(
            reverse_raster,
            origin='lower',
            aspect='auto',
            interpolation='nearest',
            extent=(0, width, -len(reverse_sorted), 0),
            **kwargs
        )
        
    ax.set_xlim(0, width)
    ax.set_ylim(-len(reverse_sorted), len(forward_sorted))
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Labels for orientation are useful when the two halves have very
    # different numbers of reads.
    ax.text(0.1, 0.99, f"Forward (n={len(forward_sorted)})",
            transform=ax.transAxes, ha='left', va='top')
    ax.text(0.1, 0.01, f"Reverse (n={len(reverse_sorted)})",
            transform=ax.transAxes, ha='left', va='bottom')

    if labels:
        ax.legend(handles=labels, bbox_to_anchor=(1.25, 0.6), loc='best', frameon=False)

    return f, ax


def plotWaterfall(array, xlabel, ylabel, labels=None, colorbar=False, ax=None, **kwargs):
    if ax is None:
        f, ax = plt.subplots()
    else:
        f = plt.gcf()
    image = ax.imshow(array, origin='lower', aspect='auto', interpolation='nearest', **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if labels:
        ax.legend(handles=labels, bbox_to_anchor=(1.25, 0.6), loc='best', frameon=False)
    if colorbar:
        f.colorbar(image, ax=ax, shrink=0.6, label=colorbar)
    return f, ax
    
def motifRaster(recs,motifs,colors):
    if not recs:
        return np.empty((0, 0, 3))

    width = max(len(rec.sequence) for rec in recs)
    raster = np.ones((len(recs), width, 3))
    patt = re.compile('|'.join(['(%s)'%m for m in motifs]))
    colors['other'] = UNKNOWN

    for i,rec in enumerate(recs):
        for j in patt.finditer(rec.sequence):
            raster[i,j.start():j.end(),:] = colors[j.group()]

        # Fill in unknown cells only within the read.
        blank = np.all(raster[i, :len(rec.sequence)] == BLANK, axis=1)
        raster[i, :len(rec.sequence)][blank, :] = colors['other']
    return raster

def qvRaster(recs,title):
    raster  = -np.ones((len(recs),len(recs[0].sequence)))
    for i,rec in enumerate(recs):
        raster[i,:len(rec.sequence)] = phred2QV(rec)
    return np.ma.masked_where(raster == -1,raster)

def phred2QV(rec):
    return np.array([ord(q)-33 for q in rec.quality])

class Waterfall_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='waterfall_Monckton.py', description='Quick waterfall plot from sequencing reads in fasta/fastq format. Some options incompatible with fasta format.')
    parser.add_argument('-i,--inFastx', dest='inFastx', type=str, default=None, required=True, 
                    help='Input Fastx file. (required)')
    parser.add_argument('-o,--out', dest='out', type=str, default=None, required=True,
                    help='Output file. (required)')
    parser.add_argument('-M,--motifs', dest='motifs', type=str, default=None,
                    help='Search motifs, comma separated, most frequent first, e.g. \'CGG,AGG\'. Either motifs (-m) or locus (-L) must be specified.')
    parser.add_argument('-L,--locus', dest='locus', type=str, default=None,
                    help='By specifying the locus, the motifs and colours are defined. Additional motifs to be coloured can be specified using the -m flag. Either motifs (-m) or locus (-L) must be specified. (Options: HTT, DMPK1)')    
    parser.add_argument('-q,--plotQV', dest='plotQV', action='store_true', default=False,
                    help='Plot additional QV waterfall (visualising base quality of reads).  Default False')
    parser.add_argument('-Q,--quality', dest='min_quality_threshold', type=int, default=0,
                    help='Set average read quality filter (Q score), see -t for other quality scoring options. Default 0')
    parser.add_argument('-t,--quality-type', dest='quality_type', type=str, default='Average',
                    help='Set the read quality scoring method. Options: Average, First200. Default Average')
    parser.add_argument('-S,--sort-type', dest='sort_type', type=str, default=None,
                    help='Set how reads are sorted. Options: repeat_length, repeat_count. Default Read length')
    parser.add_argument('-R,--reverse-complement', dest='fastq_RC', action='store_true',
                    help='Set whether reads should be reverse complemented based on motif counts. Default False')
    parser.add_argument('--split-direction', dest='split_direction', action='store_true', default=False,
                    help='Plot forward reads above the x-axis and reverse reads below it. Forward reads are sorted shortest-to-longest; reverse reads longest-to-shortest. Requires -R to identify reads by motif orientation.')
    parser.add_argument('-a,--allele_length', dest='allele_length', type=str, default=None,
                    help='Set whether reads are split between two graphs based on allele. Set with the number of repeats in the short allele. Not compatible with -q. Default None')
    parser.add_argument('-C,--colours', dest='colours', type=str, default=None,
                    help='Specify colours for motifs using comma-separated HEX codes. Should match length and order of motifs. Default colours from the \'tab10\' colour palette.')
    parser.add_argument('-max,--max-read-length', dest='max_length', type=int, default=None, 
                        help='Do not plot reads above this length. Default None')
    parser.add_argument('-min,--min-read-length', dest='min_length', type=int, default=None, 
                        help='Do not plot reads above this length. Default None')
    parser.add_argument('-F,--flanks', dest='flanks', type=str, default=None,
                    help='Comma separated upstream and downstream flanking motifs. Can be used to colour-code flanks') ##Don't think functionality beyond here is implemented ##, remove reads without either or both flanks, or sort by downstream flank position. Default None')
    #parser.add_argument('-f,--flank-filtering', dest='flank-filtering', type=str, default=None,
    #                help='Keep reads based Apply filters. Options: up, down, both')
    parser.add_argument('-y,--ylabel', dest='ylabel', type=str, default=YLABEL, required=False,
                    help='Y-axis label for plot. Default %s'%YLABEL)
    parser.add_argument('-f,--format', dest='format', type=str, default='png',
                    help='Image format.  Default png')
    parser.add_argument('-d,--dpi', dest='dpi', type=int, default=400,
                    help='Image resolution.  Default 400')
    
    try:
        main(parser)
    except Waterfall_Exception as e:
        print(f'ERROR:{e}')
        sys.exit(1) 
