import random
import numpy
from sys import getsizeof

#Constants
L = 50 #100?
K = 20 #30?
G = 10000000 #G
c = 10 #C: 20?
N = int(c*G/L)

class Read(object):
    """A read is a segment of DNA of length L that has been read.

    It contains information about the base pairs that were read (seq) as well
    as the k-mers that were constructed from the read (kmers).
    """

    #Dictionary of all reads by read base sequence
    reads = {}

    def __init__(self, bases):
        self.seq = bases
        self.kmers = [] #Array

kmers = {} #Dictionary of kmer_bases:kmer

class Kmer(object):
    """A K-mer is a segment of DNA of length K that has been pulled from a read.

    It contains information about the base pairs that were read (seq), the read
    that it was constructed from (read), and the K-mers which connect to it in
    the K-mer graph (in_k and out_k).
    """

    #Dictionary of all reads by read base sequence
    reads = {}

    def __init__(self, bases, prev_kmer):
        self.seq = bases
        self.reads = []
        self.in_k = []
        self.out_k = []

def link(srckmer, destkmer):
    """Constructs an edge from SRCKMER to DESTKMER.

    >>> s = make_kmer('ATCGGGTTACACCT')
    >>> d = make_kmer('ACCTAACTTACAA')
    >>> link(s, d)
    """
    #Constructs an edge between two K-mers.
    #If either doesn't exist, do nothing
    if not srckmer or not destkmer:
        return

    srcdat = kmer_bases(srckmer)
    destdat = kmer_bases(destkmer)

    #Find weight of edge
    weight = 0
    for w in range(1, min(len(srcdat), len(destdat))):
        if srcdat[-w:] == destdat[:w]:
            weight = w
    
    #Construct an edge
    edge = make_edge(weight, srckmer, destkmer)
    
    #Link source and destination K-mers
    in_edges(destkmer).append(edge)
    out_edges(srckmer).append(edge)

def add_read(read_bases):
    """Given a read READ_BASES as a string, construct a read and associated
    K-mers and add them to the database.
    """
    #If already in reads, do nothing.
    if read_bases in reads:
        return
    
    #Otherwise, make the read
    read = make_read(read_bases)
    rkmers = [] #K-mers from this read

    #Construct K-mers
    last_kmer = None
    for start_i in range(len(read_bases)-K):
        #Get bases of K-mer
        k_bases = read_bases[start_i:start_i+K]
        
        #If K-mer already present, use that
        if k_bases in kmers:
            kmer = kmers[k_bases]
        #Else construct it
        else:
            kmer = make_kmer(k_bases)
            kmers[k_bases] = kmer #Add to K-mer database
            
        #Link to previous K-mer
        link(last_kmer, kmer)
        
        #Make this K-mer the last K-mer
        last_kmer = kmer

    #Add read to database
    reads[read_bases] = read

def randDNA(length):
    """Generates an i.i.d. DNA sequence of length LENGTH.
    """
    dna = ''
    for _ in range(length):
        dna += random.choice('ACGT')
    return dna

def testdata():
    #Generate genome
    genome = randDNA(G)
    
    #Generate reads
    num_reads = 0
    while num_reads < N:
        for _ in range(10000):
            #Add a read
            read_i = random.randrange(G-L)
            add_read(genome[read_i:read_i+L])
        num_reads += 10000
        print(num_reads)
