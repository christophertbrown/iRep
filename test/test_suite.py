#!/usr/bin/env python

'''
Testing suite for iRep
'''

from subprocess import call
import os
import shutil
import glob

import iRep.iRep_filter

def load_data(datum):
    '''
    return the system path to the datum requested
    '''

    if datum == 'test_outdir':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        'tmp/test_outdir')
        return loc

    elif datum =='test_genome':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        '../sample_data/l_gasseri.fna')
        return loc

    elif datum == 'test_sams':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        '../sample_data/*.sam')
        return glob.glob(loc)

    elif datum == 'solution_irep':
        loc = os.path.join(str(os.path.dirname(os.path.realpath(__file__))), \
        '../sample_output/test.iRep.tsv')
        return loc

    else:
        raise AttributeError("datum {0} is not found".format(datum))

def execute_cmd(cmd, dry=False, shell=False, stdout=None, stderr=None):
    '''
    just a wrapper to call commands
    '''

    devnull = open(os.devnull, 'w')
    if stdout == None:
        stdout = devnull

    if stderr == None:
        stderr = devnull

    if not shell:
        print(' '.join(cmd))
        #if not dry: call(cmd, stderr=stderr, stdout=stdout)
        if not dry: call(cmd)

    else:
        print(cmd)
        if not dry: call(cmd, shell=True, stderr=stderr, stdout=stdout)

    return

class TestiRep():
    '''
    class to do regression tests of iRep
    '''

    def setUp(self):
        self.outdir = load_data('test_outdir')
        if os.path.exists(self.outdir):
            shutil.rmtree(self.outdir)
        os.makedirs(self.outdir)

        self.genome = load_data('test_genome')
        self.sams = load_data('test_sams')
        self.irep_solution = load_data('solution_irep')

    def tearDown(self):
        pass

    def run_all(self):
        self.setUp()

        self.regression_test()

        self.tearDown()

    def regression_test(self):
        '''
        just call iRep as an external script, compare output to sample output
        '''

        # generate command
        test_out = os.path.join(self.outdir, 'test.iRep')
        cmd = ['iRep', '-f', self.genome, '-o', test_out, '-s'] + self.sams

        # run command
        execute_cmd(cmd)

        # verify output
        out_pdf = test_out + '.pdf'
        assert os.stat(out_pdf).st_size > 0

        out_tsv = test_out + '.tsv'
        assert os.stat(out_tsv).st_size > 0

        sol_tsv = self.irep_solution

        sol_ireps = self.extract_ireps(sol_tsv)
        test_ireps = self.extract_ireps(out_tsv)

        assert sorted(sol_ireps) == sorted(test_ireps)

    def extract_ireps(self, tsv):
        '''
        from the path to the .tsv file, return a list of iRep values
        '''
        values = []
        irep = iRep.iRep_filter.parse_tables([tsv])
        for genome in irep.keys():
            for sam in irep[genome].keys():
                val = irep[genome][sam]['iRep']
                values.append(float("{0:.2f}".format(val)))
        return values
        
def test_short():
    '''
    tests that shouldn't take very long
    '''

def test_long():
    '''
    full test suite (could take a while)
    '''

    t = TestiRep()
    t.run_all()
    return

if __name__ == "__main__":
    test_short()
    test_long()
    print("All tests passed!")
