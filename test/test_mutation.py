#!/usr/bin/env python

from bs import mutation

import argparse
import pysam

def testargs():
    args = argparse.Namespace(single        = False,
                              requirepaired = True,
                              maxdepth      = 1000,
                              bamFileName   = 'test_data/testregion.bam',
                              snvfrac       = 0.5,
                              ignorepileup  = False,
                              refFasta      = '/home/adam/dev/genomes/broad/Homo_sapiens_assembly19.fasta',
                              phasevcf      = 'test_data/input_vcf.vcf',
                              haplosize     = 100)

    return args

if __name__ == '__main__':

    # set up dummy argparse options

    args        = testargs()
    log         = open('testlog.txt', 'w')
    #mutpos_list = [34141180, 33926580, 33926590]
    mutpos_list = [34141180, 33926550, 33926580, 33926590]
    mutid_list  = ['foo', 'bar', 'baz']
    chrom       = '22'
    mutstart    = min(mutpos_list)
    mutend      = max(mutpos_list)+1
    bamfile     = pysam.Samfile(args.bamFileName, 'rb')
    reffile     = pysam.Fastafile(args.refFasta)

    testmut = mutation.Mutation(args, log, bamfile, chrom, mutstart, mutend, mutpos_list, reffile)

    testmut.collect_reads()
    numpaired = [rp for rp in testmut.readpairs if testmut.readpairs[rp].is_paired()]

    print "collected " + str(len(testmut.readpairs)) + " read names, " + str(len(numpaired)) + " are paired"

    print "fetching mates..."
    if args.requirepaired:
        testmut.find_hanging_mates()

        numpaired = [rp for rp in testmut.readpairs if testmut.readpairs[rp].is_paired()]
        print str(len(numpaired)) + " reads are now paired"

        print "allpaired result: " + str(testmut.allpaired())

    # test partitioning on existing mutations

    testmut.get_phase_sites(args.phasevcf)
