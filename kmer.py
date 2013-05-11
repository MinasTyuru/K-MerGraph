import random
import numpy
from doctest import testmod
import time

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
    K = 30
    G = 1000000
    c = 20
    N = int(c*G/L)

    def __init__(self, bases):
        assert bases not in Read.reads
        
        self.bases = bases
        self.kmers = [] #Array of associated K-mers
        Read.reads[bases] = self

    def link(self, kmer):
        """Constructs a link between SELF and KMER.
        """
        assert kmer not in self.kmers
        assert self not in kmer.reads
        self.kmers.append(kmer)
        kmer.reads.append(self)

    def unlink(self, kmer):
        """Destroys the link between SELF and KMER.
        """
        assert kmer in self.kmers
        assert self in kmer.reads
        self.kmers.remove(kmer)
        kmer.reads.remove(self)

        #Clean up self if no longer have references
        if len(self.kmers) == 0:
            Read.reads.pop(self.bases)

    def add_read(read_bases):
        """Given a read READ_BASES as a string, construct a read and associated
        K-mers and add them to the database.
        """
        #If already in reads, do nothing.
        if read_bases in Read.reads:
            return
        
        #Otherwise, make the read
        read = Read(read_bases)

        #Construct K-mers
        prev_k = None
        for start_i in range(Read.L-Read.K+1):
            k_bases = read_bases[start_i:start_i+Read.K]
            prev_k = Kmer.add_kmer(k_bases, read, prev_k)

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

    #Dictionary of all kmers by base sequence
    kmers = {}

    def __init__(self, bases):
        self.bases = bases
        self.reads = []
        self.in_k = []
        self.out_k = []
        
        #Add to database
        assert bases not in Kmer.kmers
        Kmer.kmers[bases] = self
        
    def add_kmer(bases, read, prev_k):
        """Add a K-mer with bases BASES, extracted from read READ, and linked
        to previous K-mer PREV_K to the database, or add the link to an
        existing K-mer if necessary. Return the created/existing K-mer.
        """
        #Make new K-mer/use existing one
        if bases not in Kmer.kmers:
            k = Kmer(bases)
        else:
            k = Kmer.kmers[bases]

        #Update read and prev_k links
        if prev_k and prev_k not in k.in_k:
            k.link(prev_k)
        if k not in read.kmers:
            read.link(k)
                
        return k

    def link(self, prev_k):
        """Makes a link from PREV_K to SELF.
        """
        assert prev_k not in self.in_k
        assert self not in prev_k.out_k
        self.in_k.append(prev_k)
        prev_k.out_k.append(self)

    def unlink(self, prev_k):
        """Destroys the link from PREV_K to SELF.
        """
        assert prev_k in self.in_k
        assert self in prev_k.out_k
        self.in_k.remove(prev_k)
        prev_k.out_k.remove(self)

    def replace_in(self, old_in, new_in):
        """Replaces link from OLD_IN to SELF with a link from NEW_IN to SELF.
        """
        assert old_in in self.in_k
        assert new_in not in self.in_k
        self.unlink(old_in)
        self.link(new_in)

    def replace_out(self, old_out, new_out):
        """Replaces link from SELF to OLD_OUT with a link from SELF to NEW_OUT.
        """
        assert old_out in self.out_k
        assert new_out not in self.out_k
        old_out.unlink(self)
        new_out.link(self)

    def destroy(self):
        """Remove SELF from dictionary of K-mers.
        """
        assert self.bases in Kmer.kmers
        assert len(self.in_k) == 0, "Remaining links: %s" % self
        assert len(self.out_k) == 0, "Remaining links: %s" % self
        assert len(self.reads) == 0, "Remaining reads: %s" % self.reads[0]

        Kmer.kmers.pop(self.bases)

    def condense_all():
        """Condense edges until no more unambiguous edges remain.
        """
        condensed = 0
        while Kmer.condense_one():
            condensed += 1
            if condensed % 100000 == 0:
                print(int(time.time()) % 10000, ":",
                      condensed, "condensed,",
                      len(Kmer.kmers), "K-mers,",
                      len(Read.reads), "reads")
        print(int(time.time()) % 10000)

    def condense_one():
        """Look for an unambiguous edge and condense it.
        """
        #Find K-mers with only one outgoing edge
        for src in (k for b, k in Kmer.kmers.items() if len(k.out_k) == 1):
            #If destination node also has only one incoming edge, collapse
            dest = src.out_k[0]
            if len(dest.in_k) == 1 and src is not dest:
                return Kmer.condense(src, dest)
        
        #Return nothing if could not find unambiguous edge
        return None

    def condense(src, dest):
        """Given two K-mers with an unambiguous edge between them (SRC has only
        this outgoing edge, and DEST has only this incoming edge), condenses the
        two K-mers into one K-mer, inheriting edges and read parentage.
        """
        new_bases = Kmer.overlap(src, dest)
        collapsed = Kmer(new_bases)
        dest.unlink(src)

        #Update edges
        for k in list(src.in_k):
            k.replace_out(src, collapsed)
        for k in list(dest.out_k):
            k.replace_in(dest, collapsed)

        #Make new set of reads
        reads = set(src.reads)
        reads.intersection_update(dest.reads)

        #Only include reads that might bridge collapsed sequence
        for r in reads:
            if collapsed.bases in r.bases[1:-1]:
                r.link(collapsed)
        
        #Remove old read references
        for r in list(src.reads):
            r.unlink(src)
        for r in list(dest.reads):
            r.unlink(dest)

        #Clean up database
        src.destroy()
        dest.destroy()
        return collapsed
                
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
        read_i = random.randrange(Read.G)
        Read.add_read(genome[read_i:read_i+Read.L])
        
    print("Generated " + str(len(Read.reads)) + " unique reads.")
    print("Generated " + str(len(Kmer.kmers)) + " unique K-mers.")
    print("Condensing graph...")

    Kmer.condense_all()
        
    print("Now " + str(len(Read.reads)) + " unique reads.")
    print("Now " + str(len(Kmer.kmers)) + " unique K-mers.")

    return genome
