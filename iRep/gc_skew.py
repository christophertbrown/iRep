#!/usr/bin/env python3

"""
script for calculating gc skew

Chris Brown
ctb@berkeley.edu
"""

# python modules
import os
import sys
import argparse
import numpy as np
from scipy import signal
from itertools import cycle, product

# ctb
from iRep.fasta import iterate_fasta as parse_fasta

def plot_two(title, subtitle, A, B, labels, legend, vert = False):
    """
    plot with differnt y axes
    title = title for chart
    A = data for left axis [[x], [y]]
    B = data for right axis
    lables = [left label, right label, x label]
    legend = [[left legend], [right legend]]
    """
    import matplotlib.pyplot as plt
    plt.rcParams['pdf.fonttype'] = 42
    from matplotlib.backends.backend_pdf import PdfPages

    fig, ax1 = plt.subplots()
    colors = ['0.75', 'b', 'r', 'c', 'y', 'm', 'k', 'g']
    a_colors = cycle(colors)
    b_colors = cycle(colors[::-1])
    a_label = cycle(legend[0])
    b_label = cycle(legend[1])
    # plot left axis and x - axis
    for a in A:
        x, y = a
        ax1.set_ylabel(labels[0])
        ax1.set_xlabel(labels[-1])
        ax1.plot(x, y, c = next(a_colors), marker = 'o', ms = 4, label = next(a_label))
    # add vertical lines
    if vert is not False:
        for i in vert:
            x, c = i
            ax1.axvline(x = x, c = c, label = next(a_label), linewidth = 2)
    # plot right axis
    ax2 = ax1.twinx()
    for b in B:
        x, y = b
        ax2.set_ylabel(labels[1])
        ax2.plot(x, y, c = next(b_colors), linewidth = 2, label = next(b_label))
    xmin = min([min(i[0]) for i in A] + [min(i[0]) for i in B])
    xmax = max([max(i[0]) for i in A] + [max(i[0]) for i in B])
    ax2.set_xlim(xmin, xmax)
    # title
    plt.suptitle(title, fontsize = 16)
    plt.title(subtitle, fontsize = 10)
    # legend
    ax1.legend(loc = 'upper left', bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    plt.legend(loc = 'upper right', bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    # save
    pdf = PdfPages('%s.pdf' % title.replace(' ', '_'))
    pdf.savefig(bbox_inches = 'tight')
    plt.close()
    pdf.close()

def check_peaks(peaks, length):
    """
    select pair of min and max that are not too close or
    too far apart and have greatest y distance between one another
    """
    # if ori/ter peaks are too close or too far apart, they are probably wrong
    closest, farthest = int(length * float(0.45)), int(length * float(0.55))
    pairs = []
    for pair in list(product(*peaks)):
        tr, pk = pair # trough and peak
        a = (tr[0] - pk[0]) % length
        b = (pk[0] - tr[0]) % length
        pt = abs(tr[1] - pk[1]) # distance between values
        if (a <= farthest and a >= closest) or (b <=farthest and b >= closest):
            pairs.append([pt, tr, pk])
    if len(pairs) == 0:
        return [False, False]
    pt, tr, pk = sorted(pairs, reverse = True)[0]
    return [tr[0], pk[0]]

def find_ori_ter(c_skew, length):
    """
    find origin and terminus of replication based on 
    cumulative GC Skew
    """
    # find origin and terminus of replication based on 
    # cumulative gc skew min and max peaks
    c_skew_min = signal.argrelextrema(np.asarray(c_skew[1]), np.less, order = 1)[0].tolist()
    c_skew_max = signal.argrelextrema(np.asarray(c_skew[1]), np.greater, order = 1)[0].tolist()
    # return False if no peaks were detected
    if len(c_skew_min) == 0 or len(c_skew_min) == 0:
        return [False, False]
    else:
        c_skew_min = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_min]
        c_skew_max = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_max]
        ori, ter = check_peaks([c_skew_min, c_skew_max], length)
    return ori, ter

def gc_skew(name, length, seq, window, slide, plot_skew):
    """
    calculate gc skew and cumulative sum of gc skew over sequence windows
     gc skew = ((G - C) / (G + C)) * window size * genome length
    """
    # convert to G - C
    replacements = {'G':1, 'C':-1, 'A':0, 'T':0, 'N':0}
    gmc = [] # G - C
    for base in seq:
        try:
            gmc.append(replacements[base])
        except:
            gmc.append(0)
    # convert to G + C
    gpc = [abs(i) for i in gmc] # G + C
    # calculate sliding windows for (G - C) and (G + C)
    weights = np.ones(window)/window
    gmc = [[i, c] for i, c in enumerate(signal.fftconvolve(gmc, weights, 'same').tolist())]
    gpc = [[i, c] for i, c in enumerate(signal.fftconvolve(gpc, weights, 'same').tolist())]
    # calculate gc skew and cummulative gc skew sum
    skew = [[], []] # x and y for gc skew
    c_skew = [[], []] # x and y for gc skew cummulative sums
    cs = 0 # cummulative sum
    # select windows to use based on slide
    for i, m in gmc[0::slide]:
        p = gpc[i][1]
        if p == 0:
            gcs = 0
        else:
            gcs = m/p
        cs += gcs
        skew[0].append(i)
        c_skew[0].append(i)
        skew[1].append(gcs)
        c_skew[1].append(cs)
    ori, ter = find_ori_ter(c_skew, length)
    # plot data
    if plot_skew is True:
        title = '%s GC Skew' % (name)
        subtitle = '(window = %s, slide = %s)' % (window, slide)
        labels = ['GC Skew', 'Cumulative GC Skew', 'Position on Genome (bp)']
        # remove some points for plotting (approx. 1,000 datapoints)
        N = int(len(skew[0])/1000)
        if N != 0:
            skew = [skew[0][0::N], skew[1][0::N]]
        if ori is False:
            plot_two(title, subtitle, [skew], [c_skew], labels, \
                [[labels[0]], [labels[1]]])
        else:
            plot_two(title, subtitle, [skew], [c_skew], labels, \
                [[labels[0], 'Ori:%s' % ('{:,}'.format(ori)), \
                'Ter:%s' % ('{:,}'.format(ter))], [labels[1]]], \
                vert = [(ori, 'r'), (ter, 'b')])
    return ori, ter, skew, c_skew

def parse_genomes(fastas, single):
    """
    generator for parsing fastas
    if single is True, combine sequences in multifasta file
    """
    if single is True:
        for genome in fastas:
            sequence = []
            for seq in parse_fasta(genome): 
                sequence.extend(list(seq[1].upper()))
            yield (genome.name.rsplit('.', 1)[0], len(sequence), sequence)
    else:
        for genome in fastas:
            for seq in parse_fasta(genome):
                ID = seq[0].split('>', 1)[1].split()[0]
                yield (ID, len(seq[1]), list(seq[1].upper()))

def open_files(files):
    """
    open files in list, use stdin if first 
    item in list is '-'
    """
    if files is None:
        return files
    if files[0] == '-':
        return (sys.stdin)
    return (open(i) for i in files)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = \
            '# calculate gc skew and find Ori and Ter of replication')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', required = True, \
            help = 'fasta(s)')
    parser.add_argument(\
            '-l', default = False, type = int, \
            help = 'minimum contig length (default = 10 x window)')
    parser.add_argument(\
            '-w', default = 1000, type = int, \
            help = 'window length (default = 1000)')
    parser.add_argument(\
            '-s', default = 10, type = int, \
            help = 'slide length (default = 10)')
    parser.add_argument(\
            '--single', action = 'store_true', \
            help = 'combine multi-fasta sequences into single genome')
    parser.add_argument(\
            '--no-plot', action = 'store_false', \
            help = 'do not generate plots, print GC Skew to stdout')
    args = vars(parser.parse_args())
    fastas = open_files(args['f'])
    single, plot_skew = args['single'], args['no_plot']
    window, slide = args['w'], args['s']
    min_len = args['l']
    if min_len is False:
        min_len = 10 * window
    for name, length, seq in parse_genomes(fastas, single):
        if length < min_len:
            print('%s: Too Short' % (name), file=sys.stderr)
            continue
        ori, ter, skew, c_skew = gc_skew(name, length, seq, window, slide, plot_skew)
        if ori == False:
            ori, ter = 'n/a', 'n/a'
        else:
            ori, ter = '{:,}'.format(ori), '{:,}'.format(ter)
        print('%s -> Origin: %s Terminus: %s' \
                % (name, ori, ter), file=sys.stderr)
        if plot_skew is False:
            print('\t'.join(['# Name', 'Position', 'GC Skew', 'Cumulative GC Skew']))
            for i, pos in enumerate(skew[0]):
                out = [name, pos, skew[1][i], c_skew[1][i]]
                print('\t'.join([str(i) for i in out]))
