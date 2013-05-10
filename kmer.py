import random
import numpy
from sys import getsizeof

#Constants
L = 100 #100?
K = 12 #30?
G = 1000000 #G
c = 20 #C: 20?
N = int(c*G/L)

class Read(object):
    """A read is a segment of DNA of length L that has been read.

    It contains information about the base pairs that were read (seq) as well
    as the k-mers that were constructed from the read (kmers).
    """

    #Dictionary of all reads by read base sequence
    reads = {}

    def __init__(self, bases, kmers):
        self.bases = bases
        self.kmers = kmers #Array of associated K-mers
        Read.reads[bases] = self

    def add_read(read_bases):
        """Given a read READ_BASES as a string, construct a read and associated
        K-mers and add them to the database.
        """
        #If already in reads, do nothing.
        if read_bases in Read.reads:
            return
        
        #Otherwise, make the read        
        rkmers = set() #K-mers from this read

        #Construct K-mers
        prev_k = None
        for start_i in range(L-K):
            k_bases = read_bases[start_i:start_i+K]
            prev_k = Kmer.add_kmer(k_bases, prev_k)
            rkmers.add(prev_k)
                
        #Add read to database
        read = Read(read_bases, list(rkmers))

    def __str__(self):
        return self.bases

class Kmer(object):
    """A K-mer is a segment of DNA of length K that has been pulled from a read.

    It contains information about the base pairs that were read (seq), the read
    that it was constructed from (read), and the K-mers which connect to it in
    the K-mer graph (in_k and out_k).
    """

    #Dictionary of all kmers by base sequence
    kmers = {}

    def __init__(self, bases, prev_k):
        self.bases = bases
        self.reads = []
        self.in_k = []
        self.out_k = []

        #If there was a previous K-mer, link the two together
        self.link(prev_k)

        #Add to database
        Kmer.kmers[bases] = self
        
    def add_kmer(bases, prev_k):
        """Add a K-mer with bases BASES and linked to previous K-mer PREV_K to
        the database, or add the link to an existing K-mer if necessary. Return
        the created/existing K-mer.
        """
        
        #Make new K-mer if one does not exist
        if bases not in Kmer.kmers:
            k = Kmer(bases, prev_k)
        #Otherwise, augment existing K-mer
        else:
            k = Kmer.kmers[bases]
            k.link(prev_k)
                
        return k

    def link(self, prev_k):
        #Don't link to nothing, or existing links
        if not prev_k or prev_k in self.in_k:
            return
        
        self.in_k.append(prev_k)
        prev_k.out_k.append(self)

    def __str__(self):
        return self.bases

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
    for _ in range(N):
        #Add a read
        read_i = random.randrange(G-L)
        Read.add_read(genome[read_i:read_i+L])
