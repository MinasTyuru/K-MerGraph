import random
import numpy
from sys import getsizeof

#What happens to reads when K-mers are condensed?

#Constants
L = 100
K = 50
G = 1000000
c = 20
N = int(c*G/L)

class Read(object):
    """A read is a segment of DNA of length L that has been read.

    It contains information about the base pairs that were read (seq) as well
    as the k-mers that were constructed from the read (kmers).
    """

    #Dictionary of all reads by read base sequence
    reads = {}

    def __init__(self, bases):
        self.bases = bases
        self.kmers = None #Array of associated K-mers
        Read.reads[bases] = self

    def add_read(read_bases):
        """Given a read READ_BASES as a string, construct a read and associated
        K-mers and add them to the database.
        """
        #If already in reads, do nothing.
        if read_bases in Read.reads:
            return
        
        #Otherwise, make the read
        read = Read(read_bases)
        rkmers = set() #K-mers from this read

        #Construct K-mers
        prev_k = None
        for start_i in range(L-K):
            k_bases = read_bases[start_i:start_i+K]
            prev_k = Kmer.add_kmer(k_bases, read, prev_k)
            rkmers.add(prev_k)

        read.kmers = list(rkmers)

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

    def __init__(self, bases, reads, prev_k=None):
        self.bases = bases
        self.reads = reads
        self.in_k = []
        self.out_k = []

        #If there was a previous K-mer, link the two together
        self.link(prev_k)

        #Add to database
        Kmer.kmers[bases] = self
        
    def add_kmer(bases, read, prev_k):
        """Add a K-mer with bases BASES, extracted from read READ, and linked
        to previous K-mer PREV_K to the database, or add the link to an
        existing K-mer if necessary. Return the created/existing K-mer.
        """
        
        #Make new K-mer if one does not exist
        if bases not in Kmer.kmers:
            k = Kmer(bases, [read], prev_k)
        #Otherwise, augment existing K-mer
        else:
            k = Kmer.kmers[bases]
            k.reads.append(read)
            k.link(prev_k)
                
        return k

    def link(self, prev_k):
        #Don't link to nothing, or existing links
        if not prev_k or prev_k in self.in_k:
            return
        
        self.in_k.append(prev_k)
        prev_k.out_k.append(self)

    def condense():
        """Condense edges until no more unambiguous edges remain.
        """
        condensed = 0
        while Kmer.condense_one():
            condensed += 1
            if condensed % 10000 == 0:
                print(condensed, " condensed.")

    def condense_one():
        """Look for an unambiguous edge and condense it.
        """
        #Find K-mers with only one outgoing edge
        for src in (k for b, k in Kmer.kmers.items() if len(k.out_k) == 1):
            #If destination node also has only one incoming edge, collapse
            dest = src.out_k[0]
            if len(dest.in_k) == 1:                
                new_bases = Kmer.overlap(src, dest)
                collapsed = Kmer(new_bases, [])

                #Collapse edges and reads
                collapsed.in_k = list(src.in_k)
                collapsed.out_k = list(dest.out_k)
                for k in src.in_k:
                    k.out_k.remove(src)
                    k.out_k.append(collapsed)
                for k in dest.out_k:
                    k.in_k.remove(dest)
                    k.in_k.append(collapsed)

                #Remove src and dest references in parent reads
                for r in src.reads:
                    r.kmers.remove(src)
                for r in dest.reads:
                    r.kmers.remove(dest)
                
                reads = set(src.reads)
                reads.update(dest.reads)
                #Only include reads that might bridge collapsed sequence
                collapsed.reads = [r for r in reads if new_bases in r.bases]
                for r in collapsed.reads:
                    r.kmers.append(collapsed)

                #Update database
                Kmer.kmers[new_bases] = collapsed
                Kmer.kmers.pop(src.bases)
                Kmer.kmers.pop(dest.bases)
                return collapsed
            
        #Return nothing if could not find condensed
        return None
                
    def overlap(km1, km2):
        """Given two K-mers which have overlapping base pair sequences, return
        a sequence combining them.

        >>> k1 = Kmer('ABCDEFG')
        >>> k2 = Kmer('EFGH')
        >>> Kmer.overlap(k1, k2)
        'ABCDEFGH'
        """
        b1, b2 = km1.bases, km2.bases
        #Search for overlapping section, biggest first
        for ov in reversed(range(min(len(b1), len(b2)))):
            if b1[-ov:] == b2[:ov]:
                 return b1 + b2[ov:]

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
    print("Generating genome...")
    genome = randDNA(G)
    print("Generating reads...")
    
    #Generate reads
    for _ in range(N):
        #Add a read
        read_i = random.randrange(G-L)
        Read.add_read(genome[read_i:read_i+L])
        
    print("Generated " + str(len(Read.reads)) + " unique reads.")
    print("Generated " + str(len(Kmer.kmers)) + " unique K-mers.")
    print("Condensing graph...")

    Kmer.condense()
        
    print("Now " + str(len(Read.reads)) + " unique reads.")
    print("Now " + str(len(Kmer.kmers)) + " unique K-mers.")
