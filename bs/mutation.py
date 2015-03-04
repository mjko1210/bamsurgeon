#!/usr/bin/env python

from common import *

import subprocess
import networkx as nx

from collections import defaultdict as dd


def countBaseAtPos(bamfile,chrom,pos,mutid='null'):
    """ return list of bases at position chrom,pos
    """
    locstr = chrom + ":" + str(pos) + "-" + str(pos)
    args = ['samtools','mpileup',bamfile,'-r',locstr]

    p = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    pout = p.stdout.readlines()

    pileup = None 

    for line in pout:
        try:
            c = line.strip().split()
            assert len(c) > 5
            pileup = c[4].upper()
        except AssertionError:
            sys.stderr.write("INFO\t" + now() + "\t" + mutid + "\tmpileup failed, no coverage for base: " + chrom + ":" + str(pos) + "\n")
            return []
    bases = []
    if pileup:
        for b in pileup:
            if b in ['A','T','C','G']:
                bases.append(b)

    return bases


def makeins(read, start, ins, debug=False):
    assert len(read.seq) > len(ins) + 2
    
    if debug:
        print "DEBUG: INS: read.pos:", read.pos
        print "DEBUG: INS: start:   ", start
        print "DEBUG: INS: ins:     ", ins
        print "DEBUG: DEL: cigar:     ", read.cigarstring

    orig_len = len(read.seq)
    pos_in_read = start - read.pos + read.qstart


    if pos_in_read > 0: # insertion start in read span
        if debug:
            print "DEBUG: INS: pos_in_read:", pos_in_read
        left  = read.seq[:pos_in_read]
        right = read.seq[pos_in_read:]

        newseq = left + ins + right
        newseq = newseq[:orig_len]

    else: # insertion continues to the left of read
        right = read.seq[pos_in_read:]
        newseq = ins + right
        newseq = newseq[-orig_len:]

    if debug:
        print "DEBUG: INS: orig seq:", read.seq
        print "DEBUG: INS: newseq:  ", newseq
    return newseq


def makedel(read, chrom, start, end, ref, debug=False):
    assert len(read.seq) > end-start-2
    
    if debug:
        print "DEBUG: DEL: read.pos:", read.pos
        print "DEBUG: DEL: start:   ", start
        print "DEBUG: DEL: ins:     ", end
        print "DEBUG: DEL: cigar:     ", read.cigarstring
        print "DEBUG: DEL: orig seq:     ", read.seq

    orig_len = len(read.seq)
    #orig_end = read.pos + orig_len
    start_in_read = start - read.pos + read.qstart
    end_in_read = end - read.pos + read.qstart

    if debug:
        print "DEBUG: DEL: start_in_read:", start_in_read
        print "DEBUG: DEL: end_in_read:  ", end_in_read

    if start_in_read < 0: # deletion begins to the left of the read
        if debug:
            print "DEBUG: DEL: del begins to left of read." 

        assert end_in_read < orig_len
        right = read.seq[end_in_read:]
        left  = ref.fetch(chrom, start-(len(read.seq) - len(right)), start)

    elif end_in_read > orig_len: # deletion ends to the right of the read
        if debug:
            print "DEBUG: DEL: del ends to right of read."

        assert start_in_read > 0
        left  = read.seq[:start_in_read]
        right = ref.fetch(chrom, end, end+(len(read.seq) - len(left)))

    else:
        if debug:
            print "DEBUG: DEL: del starts and ends within read." 

        assert end_in_read <= orig_len and start_in_read >= 0 # deletion contained within the read
        left  = read.seq[:start_in_read]
        right = read.seq[end_in_read:]
        right += ref.fetch(chrom, read.pos+len(read.seq), read.pos+len(read.seq)+(len(read.seq)-len(left)-len(right)))

    if debug:
        print "DEBUG: DEL:  newseq:     ", left + right
    return left + right


class ReadPair:
    def __init__(self, read):
        assert not read.is_secondary and bin(read.flag & 2048) != bin(2048), "cannot process non-primary alignments!"
        self.read1 = None
        self.read2 = None

        if read.is_read1:
            self.read1 = read
        elif read.is_read2:
            self.read2 = read
        else:
            raise ValueError("encountered read " + self.name + " without is_read1 or is_read2 set!")

        self.name = read.qname

    def is_paired(self):
        return self.read1 is not None and self.read2 is not None

    def addmate(self, read):
        if read.is_read1 and self.read1 is None:
            self.read1 = read

        elif read.is_read2 and self.read2 is None:
            self.read2 = read


class VCFMut:
    def __init__(self, vcfline):
        self.vcfline = vcfline.strip()
        chrom, pos, siteid, alt, ref, alt, qual, filt = self.vcfline.split()[:8]
        
        self.pos    = int(pos)
        self.chrom  = chrom
        self.alt    = alt.split(',')
        self.filter = True

        if filt in ('PASS', '.'):
            self.filter=False

    def in_range(self, chrom, start, end):
        if self.chrom == chrom and start <= self.pos and end >= self.pos:
            return True
        else:
            return False

    def __lt__(self, other):
        if self.chrom == other.chrom:
            return self.pos < other.pos
        else:
            return self.chrom < other.chrom

    def __eq__(self, other):
        return self.vcfline == other.vcfline

    def __str__(self):
        return ':'.join(self.vcfline.split()[:5])



def nearest_site(chrom, pos, vcfmuts):
    mindist = None
    minsite = None

    for mut in vcfmuts:
        if chrom == mut.chrom:
            if minsite is None:
                mindist = abs(mut.pos - pos)
                minsite = mut
            else:
                if abs(mut.pos - pos) < mindist:
                    mindist = abs(mut.pos - pos)
                    minsite = mut

    return minsite


class Mutation:
    def __init__(self, args, log, bamfile, chrom, mutstart, mutend, mutpos_list, reffile):
        self.args        = args
        self.log         = log
        self.bamfile     = bamfile
        self.reffile     = reffile
        self.chrom       = chrom
        self.mutstart    = mutstart
        self.mutend      = mutend
        self.mutpos_list = mutpos_list # assume all on the same chromosome
        self.failed      = False       # flag for rejecting mutations

        assert mutend > mutstart, "mutation start must occur before mutation end: " + mutid

        self.region = 'haplo_' + chrom + '_' + str(self.mutstart) + '_' + str(self.mutend)

        self.readpairs = {}

        self.avoid        = None
        self.mutid_list   = None
        self.is_snv       = False
        self.mutbase_list = None
        self.is_insertion = False
        self.is_deletion  = False
        self.ins_seq      = None
        self.indel_start  = None
        self.indel_end    = None


    def collect_reads(self):
        ''' fetch all reads covering mutation region '''
        for pcol in self.bamfile.pileup(reference=self.chrom, start=self.mutstart, end=self.mutend):
            for pread in pcol.pileups:
                if not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048) and not pread.alignment.mate_is_unmapped: # only consider primary alignments

                    if pread.alignment.qname not in self.readpairs:
                        self.readpairs[pread.alignment.qname] = ReadPair(pread.alignment)

                    else:
                        self.readpairs[pread.alignment.qname].addmate(pread.alignment)


    def minpilepos(self):
        startposlist = []
        for name, rp in self.readpairs.iteritems():

            if rp.read1 is not None:
                startposlist.append(min(rp.read1.get_reference_positions()))

            if rp.read2 is not None:
                startposlist.append(min(rp.read2.get_reference_positions()))

        return min(startposlist)


    def maxpilepos(self):
        endposlist = []
        for name, rp in self.readpairs.iteritems():

            if rp.read1 is not None:
                endposlist.append(max(rp.read1.get_reference_positions()))

            if rp.read2 is not None:
                endposlist.append(max(rp.read2.get_reference_positions()))

        return max(endposlist)


    def reselect_by_mutation(self, invert=False):
        ''' retain entries in self.readpairs with(out) mutation in phasevcf: if >1 mutation in region choose the first one '''
        ''' invert: retain entries _without_ specified mutation '''
        pass


    def get_phase_sites(self, phasevcf):
        ''' find sites in phase vcf that are best for phasing '''
        possible_phase_sites = []

        with open(phasevcf, 'r') as vcf:
            for line in vcf:
                if not line.startswith('#'):
                    vcfmut = VCFMut(line)
                    if vcfmut.in_range(self.chrom, self.minpilepos(), self.maxpilepos()):
                        possible_phase_sites.append(vcfmut)

        print "number of possible phase sites: " + str(possible_phase_sites)

        nearest_phase_sites = {} # position in mutpos_list --> vcf record in phasing VCF
        
        # select the nearest sites to the mutation targets (1 site/target)
        for pos in self.mutpos_list:
            psite = nearest_site(self.chrom, pos, possible_phase_sites)
            if self.shared_reads(pos, psite.pos) > 0:
                nearest_phase_sites[str(psite)] = [pos, psite]

        print "nearest phase sites: "
        for rp,ns in nearest_phase_sites.iteritems():
            print "\t" + str(rp) + "-->" + str(ns)

        checked_pairs = dd(dict)
        readshare_graph = nx.DiGraph()

        phase_sites = {} # position in mutpos_list --> vcf record in phasing VCF

        for sitename1 in nearest_phase_sites:
            for sitename2 in nearest_phase_sites:

                site1 = nearest_phase_sites[sitename1][1]
                site2 = nearest_phase_sites[sitename2][1]

                if site1 != site2:
                    if (str(site1) not in checked_pairs or str(site2) not in checked_pairs[str(site1)]):
                        if (str(site2) not in checked_pairs or str(site1) not in checked_pairs[str(site2)]):
                            sr = self.shared_reads(site1.pos, site2.pos)
                            print str(site1), '--', str(site2), sr

                            checked_pairs[str(site1)][str(site2)] = sr
                            if sr > 0:
                                readshare_graph.add_edge(str(site1), str(site2), weight=sr)
                                readshare_graph.add_edge(str(site2), str(site1), weight=sr)

        # using strongly connected comp. instead of CC b/c if snpA --> snpB --> snpC
        # and not(snpA --> snpC) then reads aren't shared btwn A and C,
        # probably can't use one to to direct new mutation to the other

        scc = list(nx.strongly_connected_components(readshare_graph))
        print scc


    def find_hanging_mates(self):
        ''' clean up cases where mates were not in pileup '''
        for readname in self.readpairs:
            if not self.readpairs[readname].is_paired():

                if self.readpairs[readname].read1 is None:
                    self.readpairs[readname].addmate(self.bamfile.mate(self.readpairs[readname].read2))

                elif self.readpairs[readname].read2 is None:
                    self.readpairs[readname].addmate(self.bamfile.mate(self.readpairs[readname].read1))


    def allpaired(self):
        return (len(self.readpairs) - len([rp for rp in self.readpairs if self.readpairs[rp].is_paired()])) == 0


    def reads_with_position(self, pos):
        ''' return read names contatining aligned position pos '''
        readnames = {}
        for pcol in self.bamfile.pileup(reference=self.chrom, start=pos, end=pos+1):
            for pread in pcol.pileups:
                if pos in pread.alignment.get_reference_positions():
                    readnames[pread.alignment.qname] = True

        return readnames.keys()


    def shared_reads(self, pos1, pos2):
        ''' return number of reads shared by pair[0] and pair[1] '''
        reads_p1 = set(self.reads_with_position(pos1))
        reads_p2 = set(self.reads_with_position(pos2))

        return len(reads_p1.intersection(reads_p2))




def mutate(args, log, bamfile, bammate, chrom, mutstart, mutend, mutpos_list, avoid=None, mutid_list=None, is_snv=False, mutbase_list=None, is_insertion=False, is_deletion=False, ins_seq=None, reffile=None, indel_start=None, indel_end=None):
    assert mutend > mutstart, "mutation start must occur before mutation end: " + mutid

    outreads = {}
    mutreads = {}
    mutmates = {}

    region = 'haplo_' + chrom + '_' + str(mutstart) + '_' + str(mutend)

    for pcol in bamfile.pileup(reference=chrom, start=mutstart, end=mutend):
        if pcol.pos:
            refbase = reffile.fetch(chrom, pcol.pos-1, pcol.pos)
            basepile = ''
            for pread in pcol.pileups:
                if avoid is not None and pread.alignment.qname in avoid:
                    print "WARN\t" + now() + "\t" + region + "\tdropped mutation due to read in --avoidlist", pread.alignment.qname
                    return True, False, {}, {}, {}

                if not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048): # only consider primary alignments
                    basepile += pread.alignment.seq[pread.query_position-1]
                    pairname = 'F' # read is first in pair
                    if pread.alignment.is_read2:
                        pairname = 'S' # read is second in pair
                    if not pread.alignment.is_paired:
                        pairname = 'U' # read is unpaired

                    extqname = ','.join((pread.alignment.qname,str(pread.alignment.pos),pairname))

                    if pcol.pos in mutpos_list:
                        if not pread.alignment.is_secondary and bin(pread.alignment.flag & 2048) != bin(2048) and not pread.alignment.mate_is_unmapped:
                            outreads[extqname] = pread.alignment
                            mutid = mutid_list[mutpos_list.index(pcol.pos)]



                            if is_snv:
                                if extqname not in mutreads:
                                    mutreads[extqname] = pread.alignment.seq

                                mutbase = mutbase_list[mutpos_list.index(pcol.pos)]
                                mutbases = list(mutreads[extqname])
                                mutbases[pread.query_position-1] = mutbase
                                mutread = ''.join(mutbases)
                                mutreads[extqname] = mutread

                            if is_insertion:
                                mutreads[extqname] = makeins(pread.alignment, indel_start, ins_seq)

                            if is_deletion:
                                mutreads[extqname] = makedel(pread.alignment, chrom, indel_start, indel_end, reffile)

                            mate = None
                            if not args.single:
                                try:
                                    mate = bammate.mate(pread.alignment)
                                except:
                                    print "WARN\t" + now() + "\t" + mutid + "\twarning: no mate for", pread.alignment.qname
                                    if args.requirepaired:
                                        print "WARN\t" + now() + "\t" + mutid + "\tskipped mutation due to --requirepaired"
                                        return True, False, {}, {}, {}

                            if extqname not in mutmates:
                                mutmates[extqname] = mate

                            log.write(" ".join(('read',extqname,mutreads[extqname],"\n")))

                        if len(mutreads) > int(args.maxdepth):
                            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tdepth at site is greater than cutoff, aborting mutation.\n")
                            return True, False, {}, {}, {}

            # make sure region doesn't have any changes that are likely SNPs
            # (trying to avoid messing with haplotypes)
            maxfrac = 0.0
            hasSNP  = False

            basepile = countBaseAtPos(args.bamFileName,chrom,pcol.pos,mutid=region)
            if basepile:
                majb = majorbase(basepile)
                minb = minorbase(basepile)

                frac = float(minb[1])/(float(majb[1])+float(minb[1]))
                if minb[0] == majb[0]:
                    frac = 0.0
                if frac > maxfrac:
                    maxfrac = frac
                if frac > float(args.snvfrac):
                    sys.stderr.write("WARN\t" + now() + "\t" + region + "\tdropped for proximity to SNP, nearby SNP MAF: " + str(frac)  + " (max snv frac: " + args.snvfrac + ")\n")
                    hasSNP = True
            else:
                sys.stderr.write("WARN\t" + now() + "\t" + region + "\tcould not pileup for region: " + chrom + ":" + str(pcol.pos) + "\n")
                if not args.ignorepileup:
                    hasSNP = True

    return False, hasSNP, maxfrac, outreads, mutreads, mutmates # todo: convert to class
