#!/usr/bin/env python3

"""
script for combining and filtering iRep.py output
"""

# python modules
import os
import sys
import argparse
import pandas as pd

def parse_tables(tables):
    """
    parse and combine iRep output returning
    iRep, coverage, and % of windows maintained
    """
    iRep = {}
    for table in tables:
        switch = False
        for line in open(table):
            line = line.strip().split('\t')
            if line[0].startswith('## index of replication'):
                switch = True
            if switch is False:
                continue
            if line[0] == '# genome':
                samples = line[1:]
                continue
            if line[0].startswith('##'):
                if '## un-filtered index of replication' in line[0]:
                    metric = 'iRep'
                elif '## raw index of replication' in line[0]:
                    metric = 'riRep'
                elif '## r^2' in line[0]:
                    metric = 'r2'
                elif 'coverage' in line[0]:
                    metric = 'cov'
                elif 'windows' in line[0]:
                    metric = 'wins'
                elif 'fragments/Mbp' in line[0]:
                    metric = 'fragMbp'
                elif 'GC bias' in line[0]:
                    metric = 'GCB'
                elif '## GC r^2' in line[0]:
                    metric = 'GC r2'
                else:
                    continue
                continue
            genome = line[0]
            stats = []
            for i in line[1:]:
                if i == 'False' or i == 'n/a':
                    stats.append('n/a')
                else:
                    stats.append(float(i))
            for sample, stat in zip(samples, stats):
                if genome not in iRep:
                    iRep[genome] = {}
                if sample not in iRep[genome]:
                    iRep[genome][sample] = {}
                iRep[genome][sample][metric] = stat
    return iRep

def filter_iRep(iRep, min_cov = 5, min_wins = 98, min_r2 = 0.90, max_fragMbp = 175, max_GCB = False):
    """
    filter iRep values based on coverage and % of windows
    passing coverage filters
    """
    # sum coverage across genomes from the same sample in order
    # to calculate relative abundance 
    sample2covSum = {}
    for genome, samples in list(iRep.items()):
        for sample, stats in list(samples.items()):
            if sample not in sample2covSum:
                sample2covSum[sample] = 0
            if stats['cov'] != 'n/a':
                sample2covSum[sample] += stats['cov']
    # filter iRep and add relative abundance calculation
    for genome, samples in list(iRep.items()):
        for sample, stats in list(samples.items()):
            cov, wins, r2, = stats['cov'], stats['wins'], stats['r2']
            try:
                fragMbp = stats['fragMbp']
            except:
                fragMbp = False
            try:
                GCB = stats['GCB']
            except:
                GCB = False
            if cov != 'n/a' and sample2covSum[sample] != 0:
                iRep[genome][sample]['ra'] = cov/sample2covSum[sample] * 100
            else:
                iRep[genome][sample]['ra'] = 'n/a'
            iRep[genome][sample]['fiRep'] = iRep[genome][sample]['iRep']
            if \
                cov < min_cov or \
                wins < min_wins or \
                r2 < min_r2:
                    iRep[genome][sample]['fiRep'] = 'n/a'
            if fragMbp is not False:
                if fragMbp > max_fragMbp:
                    iRep[genome][sample]['fiRep'] = 'n/a'
            if GCB is not False and max_GCB is not False:
                if GCB > max_GCB:
                    iRep[genome][sample]['fiRep'] = 'n/a'
    return iRep

def print_table(ptr, metric):
    """
    print ptr values
    """
    genomes = sorted(ptr.keys())
    samples = []
    for genome in genomes:
        samples.extend(list(ptr[genome].keys()))
    samples = sorted([i for i in set(samples)])
    yield ['# genome'] + samples
    for genome, ss in list(ptr.items()):
        out = [genome]
        for sample in samples:
            if sample not in ss:
                out.append('n/a')
            else:
                try:
                    out.append(ss[sample][metric])
                except:
                    continue
        yield out

def print_short(iRep):
    """
    print short form tables
    """
    metrics = [('fiRep', 'index of replication (iRep) - thresholds: %s' % (thresholds)), \
              ('iRep', 'un-filtered index of replication (iRep)'), \
              ('riRep', 'raw index of replication (no GC bias correction)'), \
              ('r2', 'r^2'), \
              ('cov', 'coverage'), \
              ('ra', 'relative abundance'), \
              ('wins', '% windows passing filter'), \
              ('fragMbp', 'fragments/Mbp'), \
              ('GCB', 'GC bias'), \
              ('GC r2', 'GC r^2')]
    for metric in metrics:
        yield('## %s' % (metric[1]))
        for line in print_table(iRep, metric[0]):
            yield(line)

def convert_to_PD(iRep):
    """
    convert to PD long format
    """
    Gs, Ss = [], []
    iReps, fiReps, riReps, R2s, COVs, RAs, WINs, FRAGs, GCBs, GCR2s = \
            [], [], [], [], [], [], [], [], [], []
    for genome, samples in list(iRep.items()):
        for sample, stats in list(samples.items()):
            Gs.append(genome)
            Ss.append(sample)
            iReps.append(stats['iRep'])
            fiReps.append(stats['fiRep'])
            riReps.append(stats['riRep'])
            R2s.append(stats['r2'])
            COVs.append(stats['cov'])
            RAs.append(stats['ra'])
            WINs.append(stats['wins'])
            try:
                FRAGs.append(stats['fragMbp'])
            except:
                FRAGs.append('n/a')
            try:
                GCBs.append(stats['GCB'])
                GCR2s.append(stats['GC r2'])
            except:
                GCBs.append('n/a')
                GCR2s.append('n/a')
    df = pd.DataFrame({\
            'genome':Gs, 'sample': Ss, 'iRep': fiReps, \
            'un-filtered iRep': iReps, 'raw iRep (no GC bias correction)':riReps, \
            'r^2':R2s, 'coverage':COVs, \
            'relative abundance':RAs, '% windows':WINs, \
            'fragments/Mbp':FRAGs, 'GC bias':GCBs, 'GC r^2':GCR2s})
    columns = ['sample', 'genome', 'coverage', 'relative abundance', \
                'iRep', 'un-filtered iRep', 'raw iRep (no GC bias correction)', \
                'r^2', '% windows', 'fragments/Mbp', 'GC bias', 'GC r^2']
    return df, columns

if __name__ == '__main__':
    parser = argparse.ArgumentParser(\
            description = '# combine and/or filter iRep.py output')
    parser.add_argument(\
            '-t', nargs = '*', action = 'store', required = True, \
            help = 'iRep table(s)')
    parser.add_argument(\
            '-c', default = float(5), type = float, \
            help = 'min. coverage (default = 5)')
    parser.add_argument(\
            '-w', default = float(98), type = float, \
            help = 'min. percent windows passing filter (default = 98)')
    parser.add_argument(\
            '-r', default = float(0.90), type = float, \
            help = 'min. r^2 value for fitting (default = 0.90)')
    parser.add_argument(\
            '-f', default = float(175), type = float, \
            help = 'max. fragments/Mbp (default = 175)')
    parser.add_argument(\
            '-g', default = False, type = float, \
            help = 'max. GC bias (default = no filter)')
    parser.add_argument(\
            '--long', action = 'store_true', \
            help = 'print in long format')
    args = vars(parser.parse_args())
    iRep = parse_tables(args['t']) # iRep[genome][sample] = {wins, cov, ri, iRep}
    iRep = filter_iRep(iRep, \
                min_cov = args['c'], min_wins = args['w'], \
                min_r2 = args['r'], max_fragMbp = args['f'], max_GCB = args['g'])
    thresholds = \
            'min cov. = %s, min wins. = %s, min r^2 = %s, max fragments/Mbp = %s' % \
            (args['c'], args['w'], args['r'], args['f'])
    if args['long'] is False:
        for line in print_short(iRep):
            if line[0] == '#':
                print(line)
            else:
                print('\t'.join([str(i) for i in line]))
    else:
        df, columns = convert_to_PD(iRep)
        print('# thresholds: %s' % (thresholds))
        print((df.to_csv(sep = '\t', columns = columns)))
