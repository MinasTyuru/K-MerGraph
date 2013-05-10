import random
import numpy
import doctest
import time

#What happens to reads when K-mers are condensed?


class Read(object):
    """A read is a segment of DNA of length L that has been read.

    It contains information about the base pairs that were read (seq) as well
    as the k-mers that were constructed from the read (kmers).

    >>> Read.reads = {} #Clear databases
    >>> Kmer.kmers = {}
    >>> Read.L, Read.K = 8, 4 #Simulate CAGAATCCACAACAGAATTA
    >>> Read.add_read('CAGAATCC')
    >>> print(Read.reads['CAGAATCC'])
    CAGAATCC: [AATC, AGAA, ATCC, CAGA, GAAT]
    >>> print(Kmer.kmers['AATC'])
    [GAAT] => AATC => [ATCC]
    >>> print(Kmer.kmers['CAGA'])
    [] => CAGA => [AGAA]
    >>> Read.add_read('CAGAATTA') #Try a repeat
    >>> print(Read.reads['CAGAATTA'])
    CAGAATTA: [AATT, AGAA, ATTA, CAGA, GAAT]
    >>> print(Kmer.kmers['GAAT'])
    [AGAA] => GAAT => [AATC, AATT]
    """

    #Dictionary of all reads by read base sequence
    reads = {}
    
    #Constants
    L = 100
    K = 50
    G = 1000000
    c = 20
    N = int(c*G/L)

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
        for start_i in range(Read.L-Read.K+1):
            k_bases = read_bases[start_i:start_i+Read.K]
            prev_k = Kmer.add_kmer(k_bases, read, prev_k)
            rkmers.add(prev_k)

        read.kmers = list(rkmers)

    def __str__(self):
        out_str = self.bases + ": ["
        out_str += ", ".join(sorted(k.bases for k in self.kmers))
        out_str += "]"
        return out_str

class Kmer(object):
    """A K-mer is a segment of DNA of length K that has been pulled from a read.

    It contains information about the base pairs that were read (seq), the read
    that it was constructed from (read), and the K-mers which connect to it in
    the K-mer graph (in_k and out_k).

    >>> Read.reads = {} #Clear databases
    >>> Kmer.kmers = {}
    >>> genome = 'CAGAATCCACAACAGAATTACAGAATCC' #Genome should be cyclic
    >>> Read.L, Read.K = 8, 4 #Simulate full read
    >>> for read_i in range(len(genome)-Read.L+1):
    ...     Read.add_read(genome[read_i:read_i+Read.L])
    >>> print(Read.reads['CAGAATTA'])
    CAGAATTA: [AATT, AGAA, ATTA, CAGA, GAAT]
    >>> print(Kmer.kmers['GAAT'])
    [AGAA] => GAAT => [AATC, AATT]
    >>> Kmer.condense_all()
    >>> sorted([b for b, k in Kmer.kmers.items()])
    ['AATCCACAACA', 'AATTACA', 'ACAGAAT']
    >>> print(Kmer.kmers['ACAGAAT'])
    [AATCCACAACA, AATTACA] => ACAGAAT => [AATCCACAACA, AATTACA]
    """

    """
    Now we condense the graph.
    >>> k = Kmer.condense(Kmer.kmers['AATC'], Kmer.kmers['ATCC'])
    >>> k = Kmer.condense(Kmer.kmers['AATCC'], Kmer.kmers['TCCA'])
    >>> k = Kmer.condense(Kmer.kmers['AATCCA'], Kmer.kmers['CCAC'])
    >>> k = Kmer.condense(Kmer.kmers['AATCCAC'], Kmer.kmers['CACA'])
    >>> k = Kmer.condense(Kmer.kmers['AATCCACA'], Kmer.kmers['ACAA'])
    >>> k = Kmer.condense(Kmer.kmers['AATCCACAA'], Kmer.kmers['CAAC'])
    >>> k = Kmer.condense(Kmer.kmers['AATCCACAAC'], Kmer.kmers['AACA'])
    >>> k = Kmer.condense(Kmer.kmers['AATCCACAACA'], Kmer.kmers['ACAG'])
    >>> k = Kmer.condense(Kmer.kmers['AATCCACAACAG'], Kmer.kmers['CAGA'])
    >>> k = Kmer.condense(Kmer.kmers['AATCCACAACAGA'], Kmer.kmers['AGAA'])
    >>> k = Kmer.condense(Kmer.kmers['AATT'], Kmer.kmers['ATTA'])
    >>> print(Kmer.kmers['GAAT'])
    [AATCCACAACAGAA] => GAAT => [AATCCACAACAGAA, AATT]
    >>> print(Kmer.kmers['AATCCACAACAGAA'])
    [GAAT] => AATCCACAACAGAA => [GAAT]
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

    def condense_all():
        """Condense edges until no more unambiguous edges remain.
        """
        condensed = 0
        while Kmer.condense_one():
            condensed += 1
            if condensed % 100000 == 0:
                print(time.time(), ": ", condensed, " condensed. ",
                      len(Kmer.kmers), " K-mers remain. "
                      len(Read.reads), " reads remain. ")

    def condense_one():
        """Look for an unambiguous edge and condense it.
        """
        #Find K-mers with only one outgoing edge
        for src in (k for b, k in Kmer.kmers.items() if len(k.out_k) == 1):
            #If destination node also has only one incoming edge, collapse
            dest = src.out_k[0]
            if len(dest.in_k) == 1:
                return Kmer.condense(src, dest)
        
        #Return nothing if could not find unambiguous edge
        return None

    def condense(src, dest):
        new_bases = Kmer.overlap(src, dest)
        collapsed = Kmer(new_bases, [])

        #Collapse edges
        collapsed.in_k = list(src.in_k)
        collapsed.out_k = list(dest.out_k)
        for k in [k for k in src.in_k if k is not src and k is not dest]:
            k.out_k.remove(src)
            k.out_k.append(collapsed)
        for k in [k for k in dest.out_k if k is not src and k is not dest]:
            k.in_k.remove(dest)
            k.in_k.append(collapsed)

        #Update reads
        reads = set(src.reads)
        reads.update(dest.reads)
        for r in reads:
            #Remove old references
            r.kmers.remove(src)
            r.kmers.remove(dest)

            #Only include reads that might bridge collapsed sequence
            if new_bases in r.bases:
                collapsed.reads.append(r)
                r.kmers.append(collapsed)
            else:
                Read.reads.pop(r)

        #Update database
        Kmer.kmers[new_bases] = collapsed
        Kmer.kmers.pop(src.bases)
        Kmer.kmers.pop(dest.bases)
        return collapsed
                
    def overlap(km1, km2):
        """Given two K-mers which have overlapping base pair sequences, return
        a sequence combining them.

        >>> k1 = Kmer('ABCDEFG', [])
        >>> k2 = Kmer('EFGH', [])
        >>> Kmer.overlap(k1, k2)
        'ABCDEFGH'
        """
        b1, b2 = km1.bases, km2.bases
        #Search for overlapping section, biggest first
        for ov in reversed(range(min(len(b1), len(b2)))):
            if b1[-ov:] == b2[:ov]:
                 return b1 + b2[ov:]

    def __str__(self):
        instr = "["+", ".join(sorted(k.bases for k in self.in_k))+"]"
        outstr = "["+", ".join(sorted(k.bases for k in self.out_k))+"]"
        return instr + " => " + self.bases + " => " + outstr

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
    genome = randDNA(Read.G)
    genome += genome[:Read.L]
    print("Generating reads...")
    
    #Generate reads
    for _ in range(Read.N):
        #Add a read
        read_i = random.randrange(Read.G-Read.L+1)
        Read.add_read(genome[read_i:read_i+Read.L])
        
    print("Generated " + str(len(Read.reads)) + " unique reads.")
    print("Generated " + str(len(Kmer.kmers)) + " unique K-mers.")
    print("Condensing graph...")

    Kmer.condense_all()
        
    print("Now " + str(len(Read.reads)) + " unique reads.")
    print("Now " + str(len(Kmer.kmers)) + " unique K-mers.")
