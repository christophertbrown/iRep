#!/usr/bin/env python3

"""
script for estimating growth rate from peak-to-trough coverage ratios
based on method from Korem et al. 2015

Chris Brown
ctb@berkeley.edu
"""

# python modules
import os
import sys
import lmfit
import argparse
import numpy as np
from scipy import signal
from scipy import ndimage
import pickle as pickle
from itertools import product
import matplotlib.pyplot as plt
from multiprocessing import Pool

# plotting modules
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['pdf.fonttype'] = 42

# ctb
from iRep.mapped import get_reads as mapped
from iRep.gc_skew import gc_skew as gc_skew
from iRep.fasta import iterate_fasta as parse_fasta

def plot_coverage(coverage, cov_label, avg_cov, fit, c_skew, ori, ter, m_filter, ptr, title):
    """
    plot coverage profiles with fit, ori, and ter
    """
    # remove some points for plotting (approx. 1,000 datapoints)
    N = int(len(coverage[0])/1000)
    if N != 0:
        X, Y = [coverage[0][0::N], coverage[1][0::N]]
    else:
        X, Y = coverage
    # plot
    fig, ax1 = plt.subplots()
    # plot coverage data
    ax1.plot(X, Y, label = 'coverage', c = '0.60', marker = 'o', ms = 2)
    # add origin and terminus
    ax1.axvline(x = ori, c = 'r', \
            label = 'consensus Ori: %s' % ('{:,}'.format(ori)), linewidth = 2)
    ax1.axvline(x = ter, c = 'b', \
            label = 'consensus Ter: %s' % ('{:,}'.format(ter)), linewidth = 2)
    # plot median filter and estimated ori and ter
    if m_filter is not None and m_filter is not False:
        ax1.plot(m_filter[0], m_filter[1], \
            label = 'median filter', c = 'k', linewidth = 3)
    # plot fit
    if fit is not False:
        ax1.plot(coverage[0], fit, label = 'least squares fitting', \
            c = 'm', alpha = 0.75, linewidth = 2)
    # plot cumulative gc skew
    if c_skew is not False:
        ax2 = ax1.twinx()
        ax2.set_ylabel('cumulative GC skew')
        ax2.plot(c_skew[0], c_skew[1], label = 'cumulative GC skew', \
            c = 'g', alpha = 0.75, linewidth = 2)
    # title
    plt.suptitle(title, fontsize = 12)
    if ptr != 'n/a':
        ptr = '%.2f' % (ptr)
    plt.title('ptr: %s     avg. cov: %.2f' % (ptr, avg_cov), fontsize = 10)
    # label axes
    ylab = cov_label
    xlab = 'position on genome (bp)'
    ax1.set_ylabel(ylab)
    ax1.set_xlabel(xlab)
    ax1.set_xlim(min(X), max(X))
    # legend
    ax1.legend(loc = 'upper right', \
            bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    plt.legend(loc = 'upper left', \
            bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    # save
    plt.close()
    return fig

def plot_genomes(genomes, plot_name):
    """
    plot coverage and fitting data for each genome and sample pair
    """
    # PDF for plost
    if '.pdf' not in plot_name:
        plot_name = '%s.pdf' % (plot_name)
    pdf = PdfPages(plot_name)
    for g_name, genome in list(genomes.items()):
        if 'c_skew' in genome:
            c_skew = genome['c_skew']
        else:
            c_skew = False
        for s_name, sample in sorted(genome['samples'].items()):
            if sample is False:
                continue
            if 's_ptr' in sample:
                fit = sample['s_ptr'][-1]
            else: 
                fit = False
            # plot coverage
            y = sample['cov']
            x = list(range(0, len(y)))
            fig = plot_coverage([x, y], 'coverage', \
                    sample['avg_cov'], False, c_skew, \
                    genome['ORI'], genome['TER'], None, \
                    sample['ptr'], 'genome: %s sample: %s' % \
                    (g_name.rsplit('.', 1)[0], s_name.rsplit('.', 1)[0]))
            pdf.savefig(fig, bbox_inches = 'tight')
            # plot filtered/transformed coverage and fitting
            if sample['filtered'] is not False:
                fig = plot_coverage(sample['filtered'], 'coverage (log2, filtered)', \
                        sample['avg_cov'], fit, c_skew, \
                        genome['ORI'], genome['TER'], sample['m_filter'], \
                        sample['ptr'], 'genome: %s sample: %s' % \
                        (g_name.rsplit('.', 1)[0], s_name.rsplit('.', 1)[0]))
                pdf.savefig(fig, bbox_inches = 'tight')
    # save PDF
    pdf.close()

def simple_plot(x, y, horiz = [], title = 'n/a'):
    # plot every 100th coverage window
    plt.plot(x, y, label = 'coverage', c = '0.75') # , marker = 'o')
    for i in horiz:
        plt.axvline(x = i)
    plt.title(title)
    plt.show()

def calc_coverage(genomes, mappings, id2g):
    """
    for each sample:
        count number of read starts at each position in genome
    # genomes[genome]['samples'][sample]['contigs'][ID] = cov
    """
    for sample, reads in mappings:
        for read in reads:
            c = read[2] # contig
            # left-most position of read on contig
            start, length = int(read[3]), len(read[9])
            end = start + length - 1
            for i in range(start - 1, end):
                try: 
                    genomes[id2g[c]]['samples'][sample]\
                            ['contigs'][c][i] += 1
                except:
                    continue
    # combine coverage data for contigs
    for genome in list(genomes.values()):
        order, samples = genome['order'], genome['samples'] 
        for sample in list(samples.values()):
            for contig in order:
                sample['cov'].extend(sample['contigs'][contig])
                del sample['contigs'][contig]
            sample['avg_cov'] = np.average(sample['cov'])
    return genomes

def coverage_function(pars, X, data = None, printPs = False): 
    """
    piecewise linear function representing 
    coverage profile across genome
    """
    results = []
    x1 = pars['x1'].value
    x2 = pars['x2'].value
    y1, y2 = pars['y1'].value, pars['y2'].value
    if y1 > y2: # y1 ~ ori
        a = float(y2 - y1) / float(x2 - x1)
    else:
        a = float(y1 - y2) / float(x1 - x2)
    if printPs is True:
        print('x1: %s x2: %s y1: %s y2: %s a:%s' \
                % ('{:,}'.format(int(x1)), '{:,}'.format(int(x2)), y1, y2, a))
    for x in X:
        if x <= x1:
            results.append(-1*a*x + y1 + a*x1)
        elif x1 < x < x2:
            results.append(a*x + y1 - a*x1)
        elif x >= x2:
            results.append(-1*a*x + y2 + a*x2)
    if data is None:
        return np.asarray(results)
    return np.asarray([y - data[i] for i, y in enumerate(results)]) # model - data

def check_peaks(peaks, length):
    """
    select pair from peaks and troughs that are not too close or
    too far apart and have greatest y distance between one another
    """
    # if ori/ter peaks are too close or too far apart, they are probably wrong
    closest, farthest = int(length * float(0.45)), int(length * float(0.55))
    pairs = []
    for pair in list(product(*peaks)):
        pk, tr = pair # peak and trough
        a = (tr[0] - pk[0]) % length
        b = (pk[0] - tr[0]) % length
        if pk[1] < tr[1]:
            continue
        peak_dist = abs(pk[1] - tr[1]) # distance between values
        if (a <= farthest and a >= closest) or (b <=farthest and b >= closest):
            pairs.append([peak_dist, pk, tr])
    if len(pairs) == 0:
        return False, False
    peak_dist, pk, tr = sorted(pairs, reverse = True)[0]
    return pk, tr

def estimate_pars(pars, window = 999):
    """
    estimate parameters for curve fitting
    """
    genome, sample, xy, length = pars
    if xy is False:
        return False
    x, y = xy
    y_med = [x, median_filter(y)]
    # find indexes of peaks and troughs
    pks = signal.find_peaks_cwt(y, np.arange(100,1000,10000))
    trs = signal.find_peaks_cwt([-i for i in y], np.arange(100,1000,10000))
    # find positions on genome for peaks and troughs
    pks = [[y_med[0][i], y_med[1][i]] for i in pks] 
    trs = [[y_med[0][i], y_med[1][i]] for i in trs]
    # find best pk/tr pair based on greatest distance in coverage 
    # and position on genome
    ori, ter = check_peaks([pks, trs], length)
    x1, x2 = ori[0], ter[0]
    y1, y2 = ori[1], ter[1]
    if genome is not None:
        return genome, sample, (x1, x2, y1, y2, y_med)
    else:
        return x1, x2, y1, y2, y_med

def fit_coverage(pars, window = 999, est_pars = False):
    """
    use coverage data to estimate parameters for
    coverage function
    """
    x, y, length = pars
    # estimate variables using median filter
    if est_pars is True:
        x1_est, x2_est, y1_est, y2_est, y_med = \
            estimate_pars((None, None, (x, y), length), window = window)
    else:
        x1_est = int(length * float(0.25))
        x2_est = int(length * float(0.75))
        y1_est, y2_est = min(y), max(y)
        y_med = [x, median_filter(y)]
    # how close can origin and terminus of 
    # replication be to one another?
    closest, farthest = \
            int(length * float(0.45)), int(length * float(0.55))
    # Parameters
    Pars = lmfit.Parameters()
    ## x1 and x2
    Pars.add('length', value = length, vary = False)
    Pars.add('x2', value = x2_est, min = 0, max = length)
    Pars.add('xdist', value = abs(x1_est - x2_est), min = closest, max = farthest)
    Pars.add('x1', expr = '(x2 + xdist) % length')
    ## y1 and y2
    Pars.add('y1', value = y1_est)
    Pars.add('y2', value = y2_est)
    # fit data to model
    mi = lmfit.minimize(coverage_function, Pars, args = (x,), \
            kws = {'data':y, 'printPs':False}, method = 'leastsq')
    # get fitted values
    fit = coverage_function(mi.params, x)
    return (x1_est, x2_est, y_med, mi.params, fit, mi.redchi)

def generate_permutations(x, y, perms):
    """
    create generator for permutations
    """
    for i in range(perms):
        i += 1
        if i > perms:
            break
        yield (x, np.random.permutation(y))

def check_sig(obs, perms, pval):
    """
    obs = observed value
    perms = list of values from permutation analysis
     check that obs is more extreme than perms
     len(extremes)/len(perms) <= pval
    """
    perms.append(obs)
    extremes = [p for p in perms if p <= obs]
    p = float(len(extremes))/float(len(perms))
    if p <= pval:
        return p, True
    return p, False

def permutation_analysis(x, y, CHI, length, perms, threads, pval = 0.05):
    """
    compare value of chi to values from permutation analysis
    """
    pool = Pool(threads)
    chi = []
    for result in pool.map(fit_coverage, [[px, py, length] \
            for px, py in generate_permutations(x, y, perms)]):
        chi.append(result[-1])
    return check_sig(CHI, chi, pval)

def log_trans(array):
    """
    log transform elements in array
    - leave 0 as 0
    """
    lt = []
    eps = 1e-50
    for i in array:
        if i < eps:
            lt.append(np.log2(eps))
        else:
            lt.append(np.log2(i))
    return lt

def zero_to_one(array):
    """
    scale array 0 -> 1
    """
    t = []
    mini = min(array)
    maxi = max(array)
    for i in array:
        if mini == maxi:
            t.append(0)
        else:
            t.append((i - mini)/(maxi - mini)*100)
    return t

def standardize(array):
    """
    standardize to have mean = 0 and stdev 1
    """
    mean = np.average(array)
    stdev = np.std(array)
    return [(i - mean)/stdev for i in array]

def find_y(X, x, y):
    """
    find y value for item in x closest to X
    """
    return y[sorted([[abs(X-p), i] for i, p in enumerate(x)])[0][1]]

def median_filter(y, window = 999):
    """
    return median filtered data
    """
    return ndimage.filters.median_filter(y, size = window, mode = 'reflect').tolist()

def median_filter_names(pars):
    """
    return names with median filter
    """
    g, s, sample = pars
    if 'filtered' not in sample:
        return g, s, False
    xy = sample['filtered']
    if xy is False:
        return g, s, False
    x, y = xy
    return g, s, [x, median_filter(y)]

def ori_from_cov(sample, error_threshold = 20000):
    """
    find x, y for ORI and TER based on
    result of fitting
    """
    g, s, length, sample = sample
    if sample['filtered'] is False:
        y = sample['cov']
        x = list(range(0, len(y)))
        sample['m_filter'] = [x, median_filter(y)]
        return g, s, sample
    x, y = sample['filtered']
    # find origin and terminus from coverage
    x1_est, x2_est, y_med, pars, fit, chi = \
            fit_coverage([x, y, length], est_pars = False)
    x1, x2 = int(pars['x1'].value), int(pars['x2'].value)
    x1_err, x2_err = pars['x1'].stderr, pars['x2'].stderr
    y1 = find_y(x1, y_med[0], y_med[1])
    y2 = find_y(x2, y_med[0], y_med[1])
    if y1 > y2:
        ORI = (x1, y1, x1_err)
        TER = (x2, y2, x2_err)
    else:
        ORI = (x2, y2, x2_err)
        TER = (x1, y1, x1_err)
    # calculate PTR from coverage at Ori and Ter
    # exclude if x1 or x2 error > error_threshold
    if x1_err > error_threshold or x2_err > error_threshold:
        ptr = False
    elif TER[1] == 0:
        ptr = False
    else:
        ptr = (2**ORI[1])/(2**TER[1])
    sample['s_ptr'] = (ptr, ORI, TER, fit)
    sample['m_filter'] = y_med
    return g, s, sample

#def calc_dist(a, b, length):
#    """
#    calculate circular distance between points
#    """
#    return min([(a - b) % length, (b - a) % length])

def sample_ptr_from_coverage(genomes, mappings, perms, threads, ptr_threshold = 1.1):
    """
    threshold = 1.1 in Korem et al 2015
    calculate ptr ratio
    """
    # calculate cumulative gc skew for each genome
    pool = Pool(threads)
    for name, ori, ter, c_skew in \
        pool.map(calc_gc_skew, [g for g in list(genomes.items())]):
            genomes[name]['c_skew'] = c_skew
    # generate list of all genomes and samples
    pairs = product(list(genomes.keys()), [i[0] for i in mappings])
    # remove samples with no coverage, also get genome length
    pairs = [(g, s, genomes[g]['len']) for g, s in pairs \
            if genomes[g]['samples'][s] is not False]
    # calculate ptr for each genome / sample pair
    pool = Pool(threads)
    # calculate ptr for each genome and sample pair
    for sample in \
            pool.map(ori_from_cov, \
            [(g, s, l, genomes[g]['samples'][s]) \
                for g, s, l in pairs]):
        g, s, sample = sample
        genomes[g]['samples'][s] = sample
    for genome in list(genomes.values()):
        length = genome['len']
        for sample in list(genome['samples'].values()):
            if sample['filtered'] is False:
                continue
            x, y = sample['filtered']
            ptr, ORI, TER, fit = sample['s_ptr']
            # make sure ptr > threshold
            if ptr is False or ptr < ptr_threshold:
                sample['s_ptr'] = (False, False, ptr, ORI, TER, fit)
            # calculate significane by comparing reduced chi square value
            # to values calculated from permutation analysis
            # only save significant Ori/Ter estimations to genome[ORI] / genome[TER]
            elif perms is not False: 
                p_vals, sig = permutation_analysis(x, y, chi, length, perms, threads)
                sample['s_ptr'] = (sig, p_vals, ptr, ORI, TER, fit)
                if sig is True:
                    genome['ORI'].append(ORI[0])
                    genome['TER'].append(TER[0])
            else:
                sample['s_ptr'] = (None, None, ptr, ORI, TER, fit) 
                genome['ORI'].append(ORI[0])
                genome['TER'].append(TER[0])
    return genomes

def calc_gc_skew(genome, window = 1000, slide = 10):
    """
    calculate gc skew
    """
    name, genome = genome
    ori, ter, skew, c_skew = gc_skew(name, genome['len'], genome['seq'], window, slide, False)
    return name, ori, ter, c_skew

def ori_from_gc_skew(genomes, mapping, threads):
    """
    find ori and ter based on gc skew
    """
    # find ori and ter for each genome based on gc skew
    pool = Pool(threads)
    for name, ori, ter, c_skew in \
        pool.map(calc_gc_skew, [g for g in list(genomes.items())]):
            genomes[name]['c_skew'] = c_skew
            genomes[name]['ORI'] = [ori]
            genomes[name]['TER'] = [ter]
    # generate list of all genomes and samples
    pool = Pool(threads)
    pairs = product(list(genomes.keys()), [i[0] for i in mappings])
    # calculate median coverage filter for each genome and sample pair
    for sample in \
        pool.map(median_filter_names, \
            [[g, s, genomes[g]['samples'][s]] for g, s in pairs]):
        g, s, y_med = sample
        if y_med is False:
            genomes[g]['samples'][s]['m_filter'] = False
        else:
            genomes[g]['samples'][s]['m_filter'] = y_med
    return genomes

def circular_median(P, G):
    """
    calculate circular median 
    see Korem et al 2015
    P = list of positions
    G = genome length
    """
    if len(P) == 1:
        return P[0]
    TmP = []
    for t in P:
        dists = []
        for p in P:
            if t != p:
                dists.append((p - t) % G)
        TmP.append([max(dists) - min(dists), t])
    Tm = sorted(TmP)[0][1]
    return (np.median([(p - Tm) % G for p in P]) + Tm) % G

def calc_ptr(genomes, min_samples = 1):
    """
    calculate PTR for each sample based on coverage at Ori and Ter
     locations re-calculated based on the circular median of Ori and Ter 
     positions determined across samples
    """
    for genome in list(genomes.values()):
        # define consensus Ori and Ter based on circular median of
        # Ori and Ter determined from each sample
        ORI, TER = genome['ORI'], genome['TER']
        if len(ORI) < min_samples or len(TER) < min_samples:
            ORIx, TERx = False, False
        else:
            ORIx = circular_median(ORI, genome['len'])
            TERx = circular_median(TER, genome['len'])
        genome['ORI'] = ORIx
        genome['TER'] = TERx
        # calculate ptr for each genome and sample based
        # on consensus Ori/Ter locations
        for name, sample in list(genome['samples'].items()):
            if sample is False:
                continue
            # check for passing significance threshold
            if ORIx is False or TERx is False or ('m_filter' not in sample) or sample['m_filter'] is False:
                sample['ptr'] = 'n/a'
            else:
                x, y = sample['m_filter']
                ORIc = find_y(ORIx, x, y)
                TERc = find_y(TERx, x, y)
                if TERc == 0:
                    sample['ptr'] = 'n/a'
                else:
                    sample['ptr'] = (2**ORIc)/(2**TERc)
    return genomes

def filter_windows(X, Y, window = 1000, inc_threshold = 0.60, \
        mdiff = float(8), wdiff = float(1.5)):
    """
    filter low and high coverage bins based on running average
    and difference from median
    # genomes[genome]['samples'][sample]['contigs'][ID] = cov
    """
    filtered = [[], []] # store x and y
    weights = np.ones(window)/window
    med = np.median(Y)
    avgs = signal.fftconvolve(Y, weights, 'same').tolist()
    for xi, avi, yi in zip(X, avgs, Y):
        # skip if zero
        if yi <= 0 or avi <= 0 or med <= 0:
            continue
        # skip if >wdiff different from average including 1000 coverage windows
        if abs(float(max([yi, avi]))/float(min([yi, avi]))) > wdiff:
            continue
        # skip if >mdiff different from median
        if abs(float(max([yi, med]))/float(min([yi, med]))) > mdiff:
            continue
        filtered[0].append(xi)
        filtered[1].append(yi)
    # > inc_threashold windows must remain
    if float(len(filtered[0])) / float(len(X)) < inc_threshold:
        return False
    return filtered

def coverage_windows(sample, window = 10000, slide = 100):
    """
    sliding window smoothing of coverage data
    # genomes[genome]['samples'][sample]['contigs'][ID] = cov
    """
    g, s, sample = sample
    # calculate coverage windows for sample
    cov = sample['cov']
    weights = np.ones(window)
    if len(cov) < len(weights):
        sample['filtered'] = False
        return (g, s, sample)
    windows = [[], []] # x and y values for windows
    i = 0
    for c in signal.fftconvolve(cov, weights, 'valid').tolist()[0::slide]:
        windows[0].append(i)
        windows[1].append(c/window)
        i += slide
    # filter high and low coverage windows
    filtered = filter_windows(windows[0], windows[1])
    if filtered is False:
        sample['filtered'] = False
    else:
        # log transform
        sample['filtered'] = [filtered[0], log_trans(filtered[1])]
    return (g, s, sample)

def calc_cov_windows(genomes, mappings, threads):
    """
    calculate coverage windows for all pairs
    of genomes and samples
    """
    # generate list of all genomes and samples
    pairs = product(list(genomes.keys()), [i[0] for i in mappings])
    pool = Pool(threads)
    # calculate coverage windows for each pair
    for sample in \
            pool.map(coverage_windows, \
            [[g, s, genomes[g]['samples'][s]] for g, s in pairs]):
        g, s, sample = sample
        genomes[g]['samples'][s] = sample
    return genomes

def parse_genomes(fastas, mappings):
    """
    read fastas into dictionary:
     genomes[genome] = {order: [contig order], {samples}}
      samples[sample] = {cov, {contigs}, window_sum, sliding_sum, sliding_average}
      contigs[contig][sample] = [cov]
    """
    id2g = {} # contig ID to genome lookup 
    genomes = {} # dictionary for saving genome info
    for genome in fastas:
        sequence = []
        name = genome.name
        samples = {s[0]:{'contigs':{}, 'cov':[]} for s in mappings}
        g = genomes[name] = {'order':[], 'samples':samples, 'ORI':[], 'TER':[]}
        g['len'] = 0
        for seq in parse_fasta(genome): 
            sequence.extend(list(seq[1].upper()))
            ID = seq[0].split('>', 1)[1].split()[0]
            g['order'].append(ID)
            id2g[ID] = name
            length = len(seq[1])
            g['len'] += length
            cov = [0 for i in range(0, length)]
            for sample in list(samples.keys()):
                g['samples'][sample]['contigs'][ID] = \
                [0 for i in range(0, length)]
        g['seq'] = sequence
    return genomes, id2g

def growth_from_ptr(fastas, mappings, out, pickle_in, pickle_out, \
        method, plot, perms, threads):
    """
    est. growth rate from peak-to-trough coverage ratio 
     1) calculate coverage
     2) smooth/normalize/filter coverage calculations
     3) find ori and ter of replication based on either:
        - fitting coverage function* to smoothed/normalized/filtered coverage 
           using Levenberg-Marquardt algorithm 
           * = piecewise linear function from Korem 2015
        - gc skew
     4) calculate PTR for each genome and sample based on consesus location of Ori/Ter
    """
    if pickle_in is False:
        # get genome info from fastas
        genomes, id2g = parse_genomes(fastas, mappings)
        # get coverage from sam files
        genomes = calc_coverage(genomes, mappings, id2g)
        # calc coverage windows
        genomes = calc_cov_windows(genomes, mappings, threads)
        # calculate per sample ptr
        if method == 'coverage':
            genomes = sample_ptr_from_coverage(genomes, mappings, perms, threads)
        elif method == 'gc_skew':
            genomes = ori_from_gc_skew(genomes, mappings, threads)
        else: 
            print('# specify method: gc_skew or coverage', file=sys.stderr)
            exit()
        # calculate ptr based on ori, ter calculated across samples
        genomes = calc_ptr(genomes)
    else:
        genomes = pickle.load(open(pickle_in, 'rb'))
    if pickle_out is not False:
        pickle.dump(genomes, open(pickle_out, 'wb'))
    if out is not False: 
        print_table(genomes, out)
    if args['plot'] is not False:
        plot_genomes(genomes, args['plot'])
    return genomes

def print_table(genomes, out):
    """
    print table of results
    """
    samples = []
    for g in list(genomes.values()):
        for s in list(g['samples'].keys()):
            samples.append(s)
    samples = sorted([i for i in set(samples)])
    header = ['# genome', 'ORI', 'TER'] + samples
    print('\t'.join(header), file=out)
    for gn, genome in list(genomes.items()):
        ptr = [gn, '{:,}'.format(genome['ORI']), '{:,}'.format(genome['TER'])]
        for sample in samples:
            sample = genome['samples'][sample]
            if sample is False:
                ptr.append(False)
            else:
                ptr.append(sample['ptr'])
        print('\t'.join([str(i) for i in ptr]), file=out)

def open_files(files):
    """
    open files in list, use stdin if first 
    item in list is '-'
    """
    if files is None:
        return files
    if files[0] == '-':
        return [sys.stdin]
    return [open(i) for i in files]

def filter_mapping(sam, mismatches, sort_sam, sbuffer):
    """
    create generator for filtering mapping
    """
    for type, read in mapped(sam, False, mismatches, \
            'both', sort_sam, False, False, sbuffer):
        if type == 10 or type == 20:
            yield read

def validate_args(args):
    """
    check that arguments are supplied correctly
    """
    if args['p'] is not False and args['p'] <= 10:
        print('# -p must be >= 10 or False (no permutation analysis)', file=sys.stderr)
        exit()
    if args['p'] is not False and type(args['p']) is not int and args['p'].lower() == 'false':
        args['p'] = False
    elif args['p'] is not False:
        args['p'] = int(args['p'])
    if args['c'] is False:
        if args['f'] is None or args['s'] is None:
            print('# -f and -s are required, or -c (-h for help)', file=sys.stderr)
            exit()
    if (args['m'] == 'gc_skew' or args['m'] == 'coverage') is False:
        print('# method may be \'gc_skew\' or \'coverage\'', file=sys.stderr)
        exit()
    # check if files already exist
    found = []
    for i in args['o'], args['pickle'], args['plot']:
        if i is not False and os.path.isfile(i) and args['ff'] is False:
            found.append(i)
    if len(found) > 0:
        print('# file(s): %s already exist. Use -ff to overwrite.' \
                % (', '.join(found)), file=sys.stderr)
        exit()
    return args

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = \
            '# est. growth rate from peak-to-trough coverage ratio')
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', required = False, \
            help = 'fasta(s)')
    parser.add_argument(\
            '-s', nargs = '*', action = 'store', required = False, \
            help = 'sorted sam file(s) for each sample (e.g.: bowtie2 --reorder)')
    parser.add_argument(\
            '-m', required = True, type = str, \
            help = 'method for detecting Ori/Ter of replication: gc_skew or coverage')
    parser.add_argument(\
            '-c', required = False, default = False, \
            help = 'pre-computed data from growth_ptr.py (optional: pickle file)')
    parser.add_argument(\
            '-o', required = True, type = str, \
            help = 'filename for output table')
    parser.add_argument(\
            '-pickle', required = False, default = False, \
            help = 'filename for output pickle file (optional)')
    parser.add_argument(\
            '-plot', required = True, type = str, \
            help = 'filename for coverage profile plots (default: no plots)')
    parser.add_argument(\
            '-mm', required = False, default = False, type = int, \
            help = 'maximum number of mapping mismatches allowed (default: no limit)')
    parser.add_argument(\
            '-p', required = False, default = False, \
            help = 'number of permutations to perform (default: None)')
    parser.add_argument(\
            '--sort', action = 'store_true', help = 'sort the sam file')
    parser.add_argument(\
            '-b', default = '100', help = 'max memory (GB) for sorting sam (default: 100)')
    parser.add_argument(\
            '-ff', action = 'store_true', help = 'overwrite files')
    parser.add_argument(\
            '-t', required = False, default = 6, type = int, \
            help = 'threads (default: 6)')
    args = vars(parser.parse_args())
    args = validate_args(args)
    fastas = open_files(args['f'])
    sams, mm, sort, sort_b = args['s'], args['mm'], args['sort'], args['b']
    out, pickle_in, pickle_out = open(args['o'], 'w'), args['c'], args['pickle']
    if sams is not None:
        mappings = [[s, filter_mapping(s, mm, sort, sort_b)] for s in sams]
    else:
        mappings = False
    genomes = growth_from_ptr(fastas, mappings, out, pickle_in, pickle_out, \
            args['m'], args['plot'], args['p'], args['t'])
