#!/usr/bin/env python3

"""
script for estimating microbial population replication rates (iRep)
from slope of coverage across complete or draft-quality genomes

Chris Brown
ctb@berkeley.edu
"""

# python modules
import os
import sys
import lmfit
import random
import argparse
import numpy as np
import pickle as Pickle
from scipy import signal
from itertools import product
from multiprocessing import Pool

# plotting modules
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
from matplotlib.backends.backend_pdf import PdfPages

# ctb
from iRep.mapped import get_reads as mapped
from iRep.fasta import iterate_fasta as parse_fasta

def plot_coverage(cov, trimmed, avg_cov, length, fit, iRep, uiRep, r2, kept_windows, title):
    """
    plot coverage profile with fit
    """
    if len(cov[0]) < 100:
        return False
    fig, ax1 = plt.subplots()
    # plot coverage data
    ax1.plot(cov[0], cov[1], label = 'filtered', c = '0.75', alpha = 0.5)
    # plot sorted coverage data
    ax1.plot(cov[0], sorted(cov[1]), label = 'sorted', c = 'k', linewidth = 3)
    # plot trimmed coverage data
    if trimmed is not False:
        ax1.plot(trimmed[0], trimmed[1], label = 'trimmed', c = 'g', linewidth = 4)
    # plot fit
    if fit is not False:
        fit, m, b = fit
        ax1.plot(fit[0], fit[1], label = 'least squares fit to trimmed', \
            c = 'm', alpha = 0.75, linewidth = 4)
        ax1.axhline(m * length + b, c = 'r', label = 'Ori/Ter')
        ax1.axhline(b, c = 'r')
    # title
    subtitle = ['iRep: %s' % (iRep), \
                'un-filtered iRep: %s' % (uiRep), \
                'r^2: %s' % (r2), \
                'avg. cov: %s' % (avg_cov), \
                '%' + ' windows: %s' % (kept_windows)]
    title = title + '\n' + '    '.join(subtitle)
    plt.suptitle(title, fontsize = 10)
    # label axes
    xlab = 'position (bp)'
    ylab = 'coverage (log2)'
    ax1.set_ylabel(ylab)
    ax1.set_xlabel(xlab)
    if cov[0] != []:
        ax1.set_xlim(min(cov[0]), max(cov[0]))
    # legend
    ax1.legend(loc = 'upper right', \
            bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    plt.legend(loc = 'upper left', \
            bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    # save
    plt.close()
    return fig

def plot_coverage_gc(sample, avg_cov, length, kept_windows, title, show = False):
    """
    plot gc bias: cov vs GC fraction
    """
    if 'GC' not in sample:
        return False
    if 'Cwindows' in sample:
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
    else:
        fig, ax1 = plt.subplots()
    GC = sample['GC']
    COV = sample['Fwindows']
    if len(GC[0]) < 100:
        return False
    # plot cov. vs. gc
    ax1.plot(GC[1][1::100], COV[1][1::100], \
            marker = 'o', \
            label = 'filtered coverage', \
            c = '0.75', alpha = 1, markersize = 5, linewidth = 0)
    # plot fit
    fit = sample['GC fit']
    if fit is not False:
        ax1.plot(fit[0], fit[1], \
                label = 'least squares fit to GC filtered data', \
                c = 'm', alpha = 0.75, linewidth = 4)
    # label axes
    xlab = 'GC content'
    ylab = 'coverage'
    ax1.set_ylabel(ylab)
    ax1.set_xlabel(xlab)
    # plot cov. vs gc for corrected data
    if 'Cwindows' in sample:
        corrected = sample['Cwindows']
        ax2 = fig.add_subplot(122)
        ax2.plot(GC[1][1::100], corrected[1][1::100], \
                    marker = 'o', label = 'corrected', \
                    c = 'g', alpha = 1, markersize = 5, linewidth = 0)
        ax2.set_xlabel(xlab)
        # plot fit to corrected data 
        Cfit = sample['Cfit']
        if Cfit is not False:
            ax2.plot(Cfit[0], Cfit[1], \
                    label = 'least squares fit to corrected data', \
                    c = 'm', alpha = 0.75, linewidth = 4)
    # title
    iRep = sample['fiRep']
    if iRep != 'n/a':
        iRep = '%.2f' % (iRep)
    if 'raw iRep' in sample:
        riRep = '%.2f' % (sample['raw iRep'])
    else:
        riRep = 'n/a'
    subtitle = ['iRep: %s' % (iRep), \
                'raw iRep: %s' % (riRep), \
                'GC r^2: %.2f' % (sample['GC r2']), \
                'GC bias: %.2f' % (sample['GC bias']), \
                'avg. cov: %s' % (avg_cov), \
                '%' + ' windows: %s' % (kept_windows)]
    title = title + '\n' + '    '.join(subtitle)
    plt.suptitle(title, fontsize = 10)
    # legend
    ax1.legend(loc = 'upper center', \
            bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    plt.legend(loc = 'upper center', \
            bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    if show is True:
        plt.show()
    # save
    plt.close()
    return fig

def plot_coverage_hist(cov, avg_cov, iRep, r2, kept_windows, title, bins = 50):
    """
    plot coverage profile with fit
    """
    # histogram of coverage data
    fig, ax1 = plt.subplots()
    data = plt.hist(cov, bins = bins, color = 'green', alpha = 0.75)
    plt.ylabel('number of windows')
    plt.xlabel('coverage')
    # title
    subtitle = ['iRep: %s' % (iRep), \
                'r^2: %s' % (r2), \
                'avg. cov: %s' % (avg_cov), \
                '%' + ' windows: %s' % (kept_windows)]
    title = title + '\n' + '    '.join(subtitle)
    plt.suptitle(title, fontsize = 10)
    # save
    plt.close()
    return fig

def plot_genomes(genomes, mappings, plot_name):
    """
    plot coverage and fitting data for each genome and sample pair
    """
    print('# plotting data', file=sys.stderr)
    sys.stderr.flush()
    # PDF for plots
    pdf = PdfPages(plot_name)
    sns.set_style('whitegrid')
    for g_name, genome in list(genomes.items()):
        for s_name in [i[0] for i in mappings]:
            if s_name not in genome['samples']:
                sample = False
            else:
                sample = genome['samples'][s_name]
            if sample is False or sample['windows'] is False:
                continue
            title = 'genome: %s sample: %s' % \
                    (g_name.rsplit('.', 1)[0], s_name.rsplit('.', 1)[0])
            x, y = sample['windows'][0], sorted(sample['windows'][1])
            # plot coverage histogram
            stats = [sample['fiRep'], sample['iRep'], sample['r2'], sample['avg_cov']]
            stats = ['%.2f' % i if i != 'n/a' else i for i in stats]
            iRep, uiRep, r2, avg_cov = stats
            kept_windows = sample['kept_windows']
            if kept_windows != 'n/a':
                kept_windows = '%.2f' % (kept_windows * 100)
            fig = plot_coverage_hist(
                    y, avg_cov, iRep, \
                    r2, kept_windows, title)
            pdf.savefig(fig, bbox_inches = 'tight')
            # plot gc bias with corrected data
            fig = plot_coverage_gc(sample, avg_cov, genome['len'], kept_windows, title)
            if fig is not False:
                pdf.savefig(fig, bbox_inches = 'tight')
            # plot filtered coverage rank with fitting
            fig = plot_coverage(
                    sample['LTwindows'], sample['trimmed'], \
                    avg_cov, genome['len'], sample['fit'], \
                    iRep, uiRep, r2, kept_windows, title)
            if fig is not False:
                pdf.savefig(fig, bbox_inches = 'tight')
    # save PDF
    pdf.close()

def simple_plot(xy, xy2 = False, horiz = [], title = 'n/a'):
    # plot every 100th coverage window
    x, y = xy
    plt.plot(x, y, label = 'coverage', c = 'red') # , marker = 'o')
    if xy2 is not False:
        x2, y2 = xy2
        plt.plot(x2, y2, label = 'coverage', c = 'blue') # , marker = 'o')
    for i in horiz:
        plt.axvline(x = i)
    plt.title(title)
    plt.savefig('simple_plot.pdf')
    plt.show()

def calc_coverage(genomes, mappings, id2g, mask_edges = True):
    """
    for each sample:
        calcualte coverage at each position in genome
    # genomes[genome]['samples'][sample]['contigs'][ID] = cov
    """
    print("# parsing mapping files", file=sys.stderr)
    sys.stderr.flush()
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
    # combine coverage and gc data for contigs
    for genome in list(genomes.values()):
        order, samples = genome['order'], genome['samples'] 
        for sample in list(samples.values()):
            for contig in order:
                try:
                    seq = sample['contigs'][contig]
                    if mask_edges is True:
                        seq = seq[100:len(seq)-100]
                    sample['cov'].extend(seq)
                except:
                    continue
            sample['avg_cov'] = np.average(sample['cov'])
    return genomes

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

def coverage_function(pars, X, data = None, printPs = False): 
    """
    linear function for sorted coverage profile
    y = mx + b
    """
    m = pars['m'].value
    b = pars['b'].value
    if printPs is True:
        print('m: %s b: %s' % \
            ('{:,}'.format(int(m)), '{:,}'.format(int(b))))
    results = [float(m * x) + b for x in X]
    if data is None:
        return np.asarray(results)
    return np.asarray([y - data[i] for i, y in enumerate(results)]) # model - data

def fit_coverage(pars):
    """
    fit line to sorted coverage values to get slope
    """
    x, y, info, return_fit = pars
    if len(x) <= 2: # make sure there is more data than parameters
        if return_fit is False:
            return (False, False, False, info)
        else:
            return (False, False, False, False, info)
    # Parameters
    Pars = lmfit.Parameters()
    ## y = mx + b
    Pars.add('m', value = 1, vary = True)
    Pars.add('b', value = 1)
    # fit data to model
    mi = lmfit.minimize(coverage_function, Pars, args = (x,), \
            kws = {'data':y, 'printPs':False}, method = 'leastsq')
    # calculate r-squared 
    r2 = 1 - (mi.residual.var() / np.var(y))
    if return_fit is False:
        return (mi.params['m'].value, mi.params['b'].value, r2, info)
    # get fitted values
    fit = [x, coverage_function(mi.params, x)]
    return (mi.params['m'].value, mi.params['b'].value, fit, r2, info)

def fragments_test(sample, threads, n = 100, \
        lengths = [100000,200000,300000,400000,500000,1000000]):
    """
    calculate slope of coverage for random fragments of different lengths
    """
    g, s, length, sample = sample
    X, Y = sample['LTwindows']
    sample['test'] = [[], [], []]
    random_fragments = []
    # length to percent of genome
    l2p = {l:str(int(float(l)/length*100)) + '%' for l in lengths}
    for l in lengths:
        if l >= length:
            continue
        possible_starts = list(range(0, length - l - 1))
        for i in range(0, n):
            start = random.choice(possible_starts)
            end = start + l
            x, y = [], []
            for xi, yi in zip(X, Y):
                if xi >= start and xi <= end:
                    x.append(xi)
                    y.append(yi)
            if len(x) > 1:
                random_fragments.append((l, x, y))
    pool = Pool(threads)
    for m, b, r2, l in pool.map(fit_coverage, \
            [(x, y, l, False) for l, x, y in random_fragments]):
        if m < 0:
            m = -m
        sample['test'][0].append('%s (%s)' % (l, l2p[l]))
        sample['test'][1].append(False)
        sample['test'][2].append(m)
    return sample

def calc_n50(sequences):
    """
    calculate n50 for list of sequences
    """
    lengths = sorted([float(len(i)) for i in sequences], reverse = True)
    total = float(sum(lengths))
    n = total * float(0.50)
    n50 = running = lengths[0]
    for length in lengths:
        if running >= n:
            return int(n50)
        else:
            n50 = length
            running += n50

def randomly_fragment(sequence, max_pieces, \
        alpha = 0.1, beta = 100000, \
        min_length = 1000, max_length = 200000):
    """
    randomly fragment genome and return
    random subset of fragments
    """
    shuffled = []
    # break into pieces of random length
    while sequence is not False:
        s = int(random.gammavariate(alpha, beta))
        if s <= min_length or s >= max_length:
            continue
        if len(sequence) < s:
            seq = sequence[0:]
        else:
            seq = sequence[0:s]
        sequence = sequence[s:]
        shuffled.append(seq)
        if sequence == []:
            break
    # shuffle pieces
    random.shuffle(shuffled)
    # subset fragments
    subset, total = [], 0
    for fragment in shuffled:
        length = len(fragment)
        if total + length <= max_pieces:
            subset.append(fragment)
            total += length
        else:
            diff = max_pieces - total
            subset.append(fragment[0:diff])
            break
    return subset

def trim_data(data, xy, p = 0.1):
    """
    remove data from ends of sorted list
    """
    if xy is False:
        length = len(data)
        num = int(length * (p/2))
        return data[num:length - num]
    X, Y = data
    length = len(X)
    num = int(length * (p/2))
    return X[num:length - num], Y[num:length - num]


def windows2iRep(windows, L, thresholds):
    """
    calculate iRep from slide window coverage calculations
    """
    # filter zero coverage windows
    Fwindows = filter_windows(windows)
    total_windows = len(windows[0])
    total_Fwindows = len(Fwindows[0])
    kept_windows = float(total_Fwindows)/float(total_windows)
    if kept_windows < thresholds['min_wins']:
        return 'n/a'
    # log transform
    x, y = [Fwindows[0], log_trans(Fwindows[1])]
    y = sorted(y)
    x, y = trim_data([x, y], xy = True)
    m, b, r2, info = fit_coverage((x, y, False, False))
    if r2 < thresholds['min_r2']:
        return 'n/a'
    return 2**(m * L) # iRep

def iRep_from_fragments(pars):
    """
    calculate irep from random genome fragments
    """

    # parameters and coverage of complete genome
    cov, length, p2l, test = pars
    p, method, window, slide, min_length, max_length, alpha_beta, mask_edges = test
    alpha, beta = alpha_beta
    test = 'method:%s window:%s slide:%s min_len:%s max_len:%s a:%s b:%s' \
                % (method, window, slide, min_length, max_length, alpha, beta)
    percent = '%s (%s)' % (int(p * 100), p2l[p])
    results = {'iRep':'n/a', 'n50':None, 'fraction':percent, \
            'fragments':None, 'method':method, \
            'window':window, 'slide':slide, 'min_length':min_length, \
            'mask_edges':mask_edges, 'test':test, 'range':'n/a'}

    # randomly fragment genome
    L = int(length * p)
    fragments = randomly_fragment(cov, L, \
                    min_length = min_length, max_length = max_length, \
                    alpha = alpha, beta = beta)

    # mask edges
    if mask_edges is True:
        fragments = [i[100:len(i)-100] for i in fragments]

    # calc n50 for random fragments
    results['n50'] = calc_n50(fragments)

    # report number of fragments
    results['fragments'] = len(fragments)

    # coverage methods

    ## cat - combine coverage values from fragments, then calculate windows
    if method == 'iRep':
        windows = [[], []] # x and y values for windows
        weights = np.ones(window)
        fragments = [base for fragment in fragments for base in fragment]
        if len(fragments) < len(weights):
            return results
        x = 1
        for y in signal.fftconvolve(fragments, weights, 'valid').tolist()[0::slide]:
            windows[0].append(x)
            windows[1].append(y/window)
            x += slide
        results['iRep'] = windows2iRep(windows, L, thresholds)
        return results

    ## median - median of 10 iReps using cat method
    elif method == 'iRep_median':
        iReps = []
        for i in range(0, 10):
            windows = [[], []] # x and y values for windows
            weights = np.ones(window)
            random.shuffle(fragments)
            combined = [base for fragment in fragments for base in fragment]
            if len(combined) < len(weights):
                return results
            x = 1
            for y in signal.fftconvolve(combined, weights, 'valid').tolist()[0::slide]:
                windows[0].append(x)
                windows[1].append(y/window)
                x += slide
            iReps.append(windows2iRep(windows, L, thresholds))
        iReps = [i for i in iReps if i != 'n/a']
        if len(iReps) == 0:
            return results
        results['iRep'] = np.median(iReps)
        results['range'] = max(iReps) - min(iReps)
        return results

    ## valid - calculate coverage windows for each fragment, then combine
    elif method == 'scaffold_windows':
        windows = [[], []] # x and y values for windows
        weights = np.ones(window)
        x = 1
        for contigCov in fragments:
            try:
                for y in signal.fftconvolve(contigCov, weights, 'valid').tolist()[0::slide]:
                    windows[0].append(x)
                    windows[1].append(y/window)
                    x += slide
            except:
                windows[0].append(x)
                windows[1].append(np.average(contigCov))
                x += slide
        results['iRep'] = windows2iRep(windows, L, thresholds) 
        return results

    ## other methods?
    else:
        print(sys.stderr, '# methods: cat or valid')
        exit()

def iRep_test(sample, thresholds, threads, n = 100, \
        fraction = [0.75], \
        method = ['iRep'], \
        window = [5000], slide = [100], \
        min_length = [1000], \
        max_length = [10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000], \
        alpha_beta = [(0.1, 21000)], \
        mask_edges = [True]):
    """
    calculate iRep for random subsets of genome
    representing different fractions of the entire genome
    as an approximation for incomplete or draft genomes
    """
    g, s, length, sample = sample
    p2l = {p:'{:,}'.format(int(float(p)*length/1000)) \
            for p in fraction} # % of genome -> length
    tests = []
    for test in product(fraction, method, window, slide, \
                            min_length, max_length, alpha_beta, mask_edges):
        for i in range(0, n):
            tests.append(test)
    sample['test'] = {'iRep':[], 'n50':[], 'fraction':[], 'fragments':[], 'method':[], \
            'window':[], 'slide':[], 'min_length':[], 'max_length':[], \
            'alpha, beta':[], 'mask_edges':[], 'test':[], 'range':[]}
    # make sure complete genome provided
    if len(sample['contigs']) != 1:
        print(sys.stderr, '# complete genome required when running tests')
        exit()
    cov = [contig for contig in list(sample['contigs'].values())][0]
    pool = Pool(threads)
    for test_results in pool.map(iRep_from_fragments, \
            [(cov, length, p2l, test) for test in tests]):
        for i, result in list(test_results.items()):
            sample['test'][i].append(result)
    return sample

def plot_tests(genomes, pairs, out, plot, cats, y_lab, normalize = False):
    """
    plot test data
    """
    lengths = []
    slopes = []
    samples = []
    n50s = []
    for g, s in pairs:
        sample = genomes[g]['samples'][s]
        s = s.rsplit('.', 1)[0].replace('_', ' ')
        l, n50, m = sample['test']
        lengths.extend(l)
        slopes.extend(m) 
        samples.extend([s for i in m])
        n50s.extend(n50)
    if normalize == 'log2':
        slopes = log_trans(slopes)
    slope_fs = pd.DataFrame({cats:lengths, y_lab:slopes, 'sample':samples, 'n50':n50s})
    slope_fs.to_csv(out, sep = '\t')
    slope_fs = slope_fs[slope_fs[y_lab] != False]
    sns.set_style('whitegrid')
    sns.set_context('poster')
    sns_plot = sns.boxplot(x = cats, y = y_lab, data = slope_fs, \
            hue = 'sample', palette = 'deep')
    sns.stripplot(x = cats, y = y_lab, data = slope_fs, \
            hue = 'sample', palette = 'deep', \
            jitter = True, size = 5, edgecolor = 'gray')
    plt.legend(loc = 'upper right', bbox_to_anchor=(1.05, 1))
    sns_plot.figure.savefig('%s' % (plot), bbox_inches = 'tight')

def print_tests(genomes, pairs, out, cats, y_lab, normalize = False):
    """
    plot test data
    """
    lengths = []
    slopes = []
    samples = []
    n50s = []
    fragments = []
    tests = []
    ranges = []
    for g, s in pairs:
        sample = genomes[g]['samples'][s]
        s = s.rsplit('.', 1)[0]
        test_results = sample['test']
        l, n50, num_fragments, m, test, rg = \
                test_results['fraction'], test_results['n50'], \
                test_results['fragments'], \
                test_results['iRep'], test_results['test'], \
                test_results['range']
        lengths.extend(l)
        slopes.extend(m) 
        samples.extend([s for i in m])
        n50s.extend(n50)
        fragments.extend(num_fragments)
        tests.extend(test)
        ranges.extend(rg)
    if normalize == 'log2':
        slopes = log_trans(slopes)
    slope_fs = pd.DataFrame({cats:lengths, y_lab:slopes, 'sample':samples, 'n50':n50s, \
                                'num. fragments':fragments, 'test':tests, 'range':ranges})
    slope_fs.to_csv(out, sep = '\t')

def test_slopes(genomes, pairs, out, plot, test, thresholds, threads):
    """
    test methods for using slope to approximate growth:
    1) 'fragments': calculate slope of genome fragemnts 
        - test min. size fragment required for reliable results
    2) 'iRep': calculate iRep based on sorted coverage windows from simulated genome fragments
        - test percent of genome required for reliable results
    """
    if test == 'fragments':
        print('# calculating coverage slope of random fragments', file=sys.stderr)
        for g, s in pairs:
            genomes[g]['samples'][s] = \
                fragments_test((g, s, genomes[g]['len'], genomes[g]['samples'][s]), threads)
        plot_tests(genomes, pairs, out, plot, \
                'fragment length (percent of genome)', 'log2(slope)', normalize = 'log2')
    if test == 'iRep':
        print('# calculating iRep for random genome subsets', file=sys.stderr)
        for g, s in pairs:
            genomes[g]['samples'][s] = \
                iRep_test((g, s, genomes[g]['len'], genomes[g]['samples'][s]), \
                    thresholds, threads)
        print_tests(genomes, pairs, out, \
                'percent of genome (length in kbp)', 'index of replication (iRep)')
    return genomes

def iRep_from_windows(X, Y, length, sample, thresholds):
    """
    calculate iRep from coverage windows
    """
    # thresholds
    min_coverage, min_windows, min_r2, maxFragMbp = \
            thresholds['min_cov'], thresholds['min_wins'], \
            thresholds['min_r2'], thresholds['fragMbp']

    # sort coverage windows
    Ys = sorted(Y)
    windows = len(Ys)
    if windows == 0:
        m = False
        sample['trimmed'] = False
    else:
        dif = float(length)/float(windows)
        Xs = [int(i * dif) + 1 for i, value in enumerate(Ys, 0)]

        # trim ends of sorted data
        Xt, Yt = trim_data((Xs, Ys), xy = True)
        sample['trimmed'] = (Xt, Yt)
        m, b, fit, r2, info = fit_coverage((Xt, Yt, None, True))

    # calculate iRep
    if m is False:
        sample['fit'] = False
        sample['r2'] = 'n/a'
        sample['iRep'] = 'n/a'
        sample['fiRep'] = 'n/a'
        return sample

    iRep = 2**(m * length)
    sample['iRep'] = iRep
    sample['fit'] = fit, m, b
    sample['r2'] = r2

    # filter iRep based on window inclusion and coverage thresholds
    if sample['kept_windows'] < min_windows or \
       sample['avg_cov'] < min_coverage or \
       sample['r2'] < min_r2 or \
       sample['fragMbp'] > maxFragMbp:
            sample['fiRep'] = 'n/a'
    else:
        sample['fiRep'] = sample['iRep']

    return sample

def iRep_calc(sample):
    """
    calculate iRep based on slope of sorted coverage values
    """
    g, s, length, sample, thresholds = sample

    # calculate fragments / Mbp
    sample['fragMbp'] = len(list(sample['contigs'].keys()))/(float(length)/1000000)

    # calculate iRep for un-corrected data
    if 'OLTwindows' in sample:
        X, Y = sample['OLTwindows']
        sample['raw iRep'] = iRep_from_windows(X, Y, length, sample, thresholds)['iRep']

    # calculate iRep
    X, Y = sample['LTwindows']
    sample = iRep_from_windows(X, Y, length, sample, thresholds)

    return g, s, sample

def calc_growth(genomes, pairs, thresholds, threads):
    """
    estimate growth from slope of sorted (filtered) coverage windows
    """
    print('# calculating coverage slope and index of replication (iRep)', file=sys.stderr)
    sys.stderr.flush()
    pool = Pool(threads)
    for g, s, sample in pool.map(iRep_calc, \
            [(g, s, genomes[g]['len'], genomes[g]['samples'][s], thresholds)
                for g, s in pairs]):
        genomes[g]['samples'][s] = sample
    pool.terminate()
    return genomes

def filter_windows(win, other = False, mdif = float(8)):
    """
    filter windows based on difference from median
    """
    X, Y = [], []
    med = np.median(win[1])
    if other is False:
        for x, y in zip(win[0], win[1]):
            if y <= 0 or med <= 0:
                continue
            if abs(float(max([y, med])) / float(min([y, med]))) > mdif:
                continue
            X.append(x)
            Y.append(y)
        return X, Y
    else:
        oX, oY = [], []
        for x, y, ox, oy in zip(win[0], win[1], other[0], other[1]):
            if y <= 0 or med <= 0:
                continue
            if abs(float(max([y, med])) / float(min([y, med]))) > mdif:
                continue
            X.append(x)
            Y.append(y)
            oX.append(ox)
            oY.append(oy)
        return [X, Y], [oX, oY]

def coverage_windows(sample, window = 5000, slide = 100):
    """
    sliding window smoothing of coverage data
    # genomes[genome]['samples'][sample]['contigs'][ID] = cov
    """
    g, gc, s, sample, gc_correction, thresholds = sample
    # calculate coverage windows for sample
    cov = sample['cov']
    del sample['cov']
    weights = np.ones(window)
    windows = [[], []] # x and y values for windows
    i = 0
    if len(cov) < len(weights):
        sample['windows'] = False
        sample['kept_windows'] = False
        sample['LTwindows'] = False
        sample['iRep'] = 'n/a'
        sample['GC bias'] = 'n/a'
        sample['fiRep'] = 'n/a'
        sample['fit'] = False
        return (g, s, sample)
    for c in signal.fftconvolve(cov, weights, 'valid').tolist()[0::slide]:
        windows[0].append(i)
        windows[1].append(c/window)
        i += slide
    # filter zero coverage windows
    sample['windows'] = windows
    if gc is not False:
        Fwindows, Fgc = filter_windows(windows, gc)
        sample['GC'] = Fgc
        sample['Fwindows'] = Fwindows # only need this for plotting GC data
    else:
        Fwindows = filter_windows(windows)
    total_windows = len(windows[0])
    total_Fwindows = len(Fwindows[0])
    sample['kept_windows'] = float(total_Fwindows)/float(total_windows)
    # calculate gc bias
    if gc is not False:
        fit, filtered, corrected_windows, Cfit = \
                gc_bias(Fgc, Fwindows, thresholds['GC_min_r2'])
        fit, m, b, r2 = fit
        sample['GC fit'], sample['GC r2'] = fit, r2
        sample['GC bias'] = m * r2
        if gc_correction is True and corrected_windows is not False:
            sample['Cwindows'] = corrected_windows
            sample['Cfit'] = Cfit
            sample['OLTwindows'] = [Fwindows[0], log_trans(Fwindows[1])]
            Fwindows = corrected_windows
    # log transform
    sample['LTwindows'] = [Fwindows[0], log_trans(Fwindows[1])]
    return (g, s, sample)

def calc_cov_windows(genomes, pairs, mappings, gc_correction, thresholds, threads):
    """
    calculate coverage windows for all pairs
    of genomes and samples
    """
    print('# calculating coverage over sliding windows', file=sys.stderr)
    sys.stderr.flush()
    # filter out any genome -> sample pairs not passing thresholds
    pairs = [(g, s) for g, s in pairs if s in genomes[g]['samples']]
    pool = Pool(threads)
    # calculate coverage windows for each pair
    for sample in \
            pool.map(coverage_windows, \
            [[g, genomes[g]['gc'], s, genomes[g]['samples'][s], gc_correction, thresholds]\
                for g, s in pairs]):
        g, s, sample = sample
        genomes[g]['samples'][s] = sample
    pool.terminate()
    # filter out any genome -> sample pairs not passing thresholds
    pairs = [(g, s) for g, s in pairs if genomes[g]['samples'][s]['windows'] is not False]
    return genomes, pairs

def gc_bias(GC, COV, correction_threshold):
    """
    calculate GC bias by fitting line to cov vs %s gc
    """
    # fit line to cov vs gc
    m, b, fit, r2, l = fit_coverage((GC[1], COV[1], False, True))
    # remove outliers
    fGC, fCOV, error = [], [], []
    for gc, cov in zip(GC[1], COV[1]):
        est = m * gc + b
        error.append(abs(cov - est))
    try:
        cutoff = sorted(error, reverse = True)[int(len(error)*0.01)]
    except:
        cutoff = 0
    for gc, cov, e in zip(GC[1], COV[1], error):
        if e >= cutoff:
            continue
        fGC.append(gc)
        fCOV.append(cov)
    # re-fit with filtered data
    m, b, fit, r2, l = fit_coverage((fGC, fCOV, False, True))
    if r2 < correction_threshold:
        return (fit, m, b, r2), (fGC, fCOV), False, False
    corrected = [[], []]
    av = np.average(COV[1])
    for p, gc, cov in zip(GC[0], GC[1], COV[1]):
        est = m * gc + b
        correction = av - est
        corrected[0].append(p)
        corrected[1].append(cov + correction)
    Cm, Cb, Cfit, Cr2, Cl = fit_coverage((fGC, corrected[1], False, True))
    return (fit, m, b, r2), [fGC, fCOV], corrected, Cfit

def gc_content(seq, window = 5000, slide = 100):
    """
    calculate gc content over sequence windows
    """
    # convert GC
    replacements = {'G':1, 'C':1, 'A':0, 'T':0, 'N':0}
    GC = [] # G - C
    for base in seq:
        try:
            GC.append(replacements[base.upper()])
        except:
            GC.append(0)
    # calculate gc content over sliding windows
    i = 0
    weights = np.ones(window)
    GCwindows = [[], []]
    for gc in signal.fftconvolve(GC, weights, 'valid').tolist()[0::slide]:
        GCwindows[0].append(i)
        GCwindows[1].append(gc/window)
        i += slide
    return GCwindows

def parse_genomes_fa(fastas, mappings, mask_edges = True):
    """
    genomes[genome name] = {order: [contig order], samples: {}}
        samples[sample name] = {cov: [coverage by position], contigs: {}}
            contigs[contig name] = [coverage by position]
    """
    id2g = {} # contig ID to genome lookup 
    genomes = {} # dictionary for saving genome info
    for genome in fastas:
        name = genome.name
        samples = {s[0]:{'contigs':{}, 'cov':[]} for s in mappings}
        g = genomes[name] = {'order':[], 'samples':samples}
        g['len'] = 0
        genome_seq = []
        for seq in parse_fasta(genome): 
            ID = seq[0].split('>', 1)[1].split()[0]
            g['order'].append(ID)
            id2g[ID] = name
            length = len(seq[1])
            g['len'] += length
            if mask_edges is True:
                contig_seq = seq[1][100:len(seq[1])-100]
            else:
                contig_seq = seq[1]
            genome_seq.extend(contig_seq)
            for sample in list(samples.keys()):
                g['samples'][sample]['contigs'][ID] = \
                    [0 for i in range(0, length)]
        g['gc'] = gc_content(genome_seq)
    return genomes, id2g

def parse_genomes_sam(id2g, mappings):
    """
    id2g = {} # contig ID to genome lookup 
    genomes[genome name] = {order: [contig order], samples: {}}
        samples[sample name] = {cov: [coverage by position], contigs: {}}
            contigs[contig name] = [coverage by position]
    """
    genomes = {} # dictionary for saving genome info
    for sam in [i[0] for i in mappings]:
        for line in open(sam):
            if line.startswith('@') is False:
                break
            if line.startswith('@SQ') is False:
                continue
            line = line.strip().split()
            ID, length = line[1].split(':', 1)[1], int(line[2].split(':', 1)[1])
            if ID not in id2g:
                continue
            name = id2g[ID]
            if name not in genomes:
                genomes[name] = {'order':[], 'samples':{}, 'len': 0, 'gc':False}
            if ID not in genomes[name]['order']:
                genomes[name]['order'].append(ID)
                genomes[name]['len'] += length
            if sam not in genomes[name]['samples']:
                genomes[name]['samples'][sam] = {'contigs':{}, 'cov':[]}
            # add cov array for each contig for each sample
            genomes[name]['samples'][sam]['contigs'][ID] = [0 for i in range(0, length)]
    return genomes

def iRep(fastas, mappings, out, pickle, plot, \
            thresholds, gc_correction, threads):
    """
    est. growth from slope of coverage
     1) calculate coverage over windows
     2) sort and filter coverage window calculations
     3) calculate slope of sorted coverage values
        - fitting line to sorted/filtered coverage
           using Levenberg-Marquardt algorithm (least squares)
     4) approximate growth from length-normalized slope
    """
    # get genome info from fastas
    if fastas is not None:
        genomes, id2g = parse_genomes_fa(fastas, mappings)
    else:
        genomes = parse_genomes_sam(id2g, mappings)

    # get coverage from sam files
    genomes = calc_coverage(genomes, mappings, id2g)

    # generate list of all genomes and samples
    pairs = [i for i in product(list(genomes.keys()), [i[0] for i in mappings])]

    # filter out any genome -> sample pairs not passing thresholds
    pairs = [(g, s) for g, s in pairs if s in genomes[g]['samples']]

    # calc coverage windows
    genomes, pairs = \
            calc_cov_windows(genomes, pairs, mappings, gc_correction, thresholds, threads)

    # calculate per sample slope and estimate growth
    genomes = calc_growth(genomes, pairs, thresholds, threads)

    # save results
    if out is not False:
        print_table(genomes, mappings, out, thresholds)
    if pickle is not False:
        pickle_out = '%s.pickle' % (out.rsplit('.', 1)[0])
        pickle_out = open(pickle_out, 'wb')
        Pickle.dump(genomes, pickle_out)
        pickle_out.close()
    if plot is not False:
        plot_genomes(genomes, mappings, plot)

    return genomes

def print_table(genomes, mappings, out, thresholds):
    """
    print table of results
    """
    print('# saving results', file=sys.stderr)
    samples = [i[0] for i in mappings]
    iRep, fiRep, riRep, r2, cov, kept, fragMbp, gcb, gcr2 = \
            [], [], [], [], [], [], [], [], []
    for name, genome in list(genomes.items()):
        iRep.append([])
        fiRep.append([])
        riRep.append([])
        r2.append([])
        cov.append([])
        kept.append([])
        fragMbp.append([])
        gcb.append([])
        gcr2.append([])
        iRep[-1].append(name)
        fiRep[-1].append(name)
        riRep[-1].append(name)
        r2[-1].append(name)
        cov[-1].append(name)
        kept[-1].append(name)
        fragMbp[-1].append(name)
        gcb[-1].append(name)
        gcr2[-1].append(name)
        for sample in samples:
            if sample not in genome['samples']:
                sample = False
            else:
                sample = genome['samples'][sample]
            if sample is False:
                iRep[-1].append('n/a')
                fiRep[-1].append('n/a')
                riRep[-1].append('n/a')
                r2[-1].append('n/a')
                cov[-1].append('n/a')
                kept[-1].append('n/a')
                fragMbp[-1].append('n/a')
                gcb[-1].append('n/a')
                gcr2[-1].append('n/a')
            else:
                iRep[-1].append(sample['iRep'])
                fiRep[-1].append(sample['fiRep'])
                try:
                    riRep[-1].append(sample['raw iRep'])
                except:
                    riRep[-1].append('n/a')
                r2[-1].append(sample['r2'])
                cov[-1].append(sample['avg_cov'])
                kept[-1].append('%.2f' % (100 * sample['kept_windows']))
                fragMbp[-1].append('%.0f' % (sample['fragMbp']))
                gcb[-1].append(sample['GC bias'])
                try:
                    gcr2[-1].append(sample['GC r2'])
                except:
                    gcr2[-1].append('n/a')
    out = open(out, 'w')
    header = ['# genome'] + samples
    thresholds = 'min cov. = %s, min wins. = %s, min r^2 = %s, max fragments/Mbp = %s, GC correction min r^2 = %s' % \
            (thresholds['min_cov'], thresholds['min_wins'], thresholds['min_r2'], \
            thresholds['fragMbp'], thresholds['GC_min_r2'])
    for vals, desc in \
            (fiRep, '## index of replication (iRep) - thresholds: %s' % (thresholds)), \
            (iRep,  '## un-filtered index of replication (iRep)'), \
            (riRep, '## raw index of replication (no GC bias correction)'), \
            (r2,   '## r^2'), \
            (cov,  '## coverage'), \
            (kept, '## % windows passing filter'), \
            (fragMbp, '## fragments/Mbp'), \
            (gcb, '## GC bias'), \
            (gcr2, '## GC r^2'):
        print(desc, file=out)
        print('\t'.join(header), file=out)
        for i in vals:
            print('\t'.join([str(j) for j in i]), file=out)
        print('#', file=out)
    out.close()
    return out

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
    if args['f'] is None and args['b'] is None:
        print('# -f or -b is required (-h for help)', file=sys.stderr)
        exit()
    # check if files already exist
    args['plot'] = '%s.pdf' % (args['o'])
    args['table'] = '%s.tsv' % (args['o'])
    found = []
    for i in args['table'], args['plot']:
        if i is not False and os.path.isfile(i) and args['ff'] is False:
            found.append(i)
    if len(found) > 0:
        print('# file(s): %s already exist. Use -ff to overwrite.' \
                % (', '.join(found)), file=sys.stderr)
        exit()
    if args['f'] is None and args['gc_correction'] is True:
        print('# -f requred for % GC correction', file=sys.stderr)
        exit()
    return args

if __name__ == '__main__':
    desc = '# calculate the Index of Replication (iRep)'
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument(\
            '-f', nargs = '*', action = 'store', required = True, \
            help = 'fasta(s)')
    parser.add_argument(\
            '-s', nargs = '*', action = 'store', required = True, \
            help = 'sorted sam file(s) for each sample (e.g.: bowtie2 --reorder)')
    parser.add_argument(\
            '-o', required = True, type = str, \
            help = 'prefix for output files (table and plots)')
    parser.add_argument(\
            '--pickle', action = 'store_true', \
            help = 'save pickle file (optional)')
    parser.add_argument(\
            '-mm', required = False, default = 1, type = int, \
            help = 'max. # of read mismatches allowed (default: 1)')
    parser.add_argument(\
            '--sort', action = 'store_true', \
            help = 'optional - sort the sam file')
    parser.add_argument(\
            '-M', default = '100', \
            help = 'max. memory (GB) for sorting sam (default: 100)')
    parser.add_argument(\
            '--no-plot', action = 'store_true', \
            help = 'do not plot output')
    parser.add_argument(\
            '--no-gc-correction', action = 'store_false', \
            help = 'do not correct coverage for GC bias before calculating iRep')
    parser.add_argument(\
            '-ff', action = 'store_true', \
            help = 'overwrite files')
    parser.add_argument(\
            '-t', required = False, default = 6, type = int, \
            help = 'threads (default: 6)')
    args = vars(parser.parse_args())
    args = validate_args(args)
    fastas = open_files(args['f'])
    sams, mm, sort, sort_b = args['s'], args['mm'], args['sort'], args['M']
    # generator for mapping
    mappings = [[s, filter_mapping(s, mm, sort, sort_b)] for s in sams]
    # cancel plotting
    if args['no_plot'] is True:
        args['plot'] = False
    # thresholds
    thresholds = {'min_cov':5, 'min_wins':0.98, 'min_r2':0.90, \
                    'fragMbp':175, 'GC_min_r2':0.0}
    # calculate iRep
    genomes = iRep(\
                fastas, mappings, \
                args['table'], args['pickle'], args['plot'],
                thresholds, args['no_gc_correction'], args['t'])
