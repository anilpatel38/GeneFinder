"""
Gene Finder

@author: Anil Patel

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

from load import load_seq
dna = load_seq("./data/X73525.fa")


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))
    # YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == "A":
        return("T")
    elif nucleotide == "T":
        return("A")
    elif nucleotide == "C":
        return("G")
    elif nucleotide == "G":
        return("C")


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("ATGCGAATGTAGCATCAAA")
    'TGAACGCGG'
    """
    index = 1
    a = ""
    while index-1 < len(dna):
        b = get_complement(dna[-1*index])
        a = a+b
        index = index + 1
    return(a)


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATdGTGAA")
    'ATdGTGAA'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    length = len(dna)
    b = "no"
    for i in range(0, length, 3):
        a = dna[i:(i+3)]
        if a == "TGA" or a == "TAG" or a == "TAA":
            b = "yes"
            return(dna[0:(i)])
        else:
            b = "no"
    if b == "no":
        return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    dnas = []
    i = 0
    while True:
        i = find_start(dna)
        if i == -1:
            break
        else:
            a = rest_of_ORF(dna[i:])
            dnas.append(a)
            temp = i + len(a)
            dna = dna[temp:]

    return dnas


def find_start(dna):
    for i in range(0, len(dna), 3):
        if dna[i:i+3] == "ATG":
            return i
    return -1


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    dnas = []
    dnas += find_all_ORFs_oneframe(dna)
    dnas += find_all_ORFs_oneframe(dna[1:])
    dnas += find_all_ORFs_oneframe(dna[2:])

    return dnas


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    dnas = find_all_ORFs(dna)
    reverse = get_reverse_complement(dna)
    dnas += find_all_ORFs(reverse)
    return dnas


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("AAAAAAAAA")
    'ATGCTACATTCGCAT'
    """
    ORFs = find_all_ORFs_both_strands(dna)
    if len(ORFs) == 0:
        return 0
    elif len(ORFs) == 1:
        return ORFs[0]
    else:
        for i in range(0, len(ORFs)-1):
            if len(ORFs[i]) < len(ORFs[i+1]):
                a = ORFs[i+1]
            else:
                a = ORFs[i]
    return a


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        >>> longest_ORF_noncoding("ATGCGAATGTAGCATCAAA", 10)
        'any'
        """
    lengths = []
    for i in range(0, num_trials):
        dna1 = shuffle_string(dna)
        if longest_ORF(dna1) == 0:
            lengths.append(0)
        else:
            lengths.append(len(longest_ORF(dna1)))
    lengths.sort()
    return lengths[-1]


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    cods = []
    AAs = []
    for i in range(0, len(dna), 3):
        a = dna[i:i+3]
        cods.append(a)
    for i in range(0, len(cods)):
        if len(cods[i]) == 3:
            AAs.append(aa_table[cods[i]])
    AAs = ''.join(AAs)
    return AAs


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    good_orfs = []
    aminos = []
    threshold = longest_ORF_noncoding(dna, 1500)
    all_orfs = find_all_ORFs_both_strands(dna)
    for orf in all_orfs:
        if len(orf) > int(threshold):
            good_orfs.append(orf)
    for i in range(0, len(good_orfs)):
        aminos.append(coding_strand_to_AA(good_orfs[i]))

    print(aminos)
    return aminos


gene_finder(dna)
"""
if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA, globals(),
                                   verbose=True)
"""
