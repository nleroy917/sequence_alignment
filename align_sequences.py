#!/usr/bin/env python

"""
    usage:
        align_sequences [options] seq1.fa seq2.fa
    where the options are:
        -h,--help : print usage and quit
        -m,--match: score of a match in the alignment [2]
        -x,--mismatch: penalty for a mismatch in the alignment [1]
        -g,--gapopen: penalty for opening a new gap [4]
        -e,--gapextend: penalty for extending a gap [1]
"""

from sys import argv, stderr
from getopt import getopt, GetoptError
import numpy as np

# a simple function to read the name and sequence from a file
# The file is expected to have just one contig/sequence. This function
# checks the assumption and complains if it is not the case.
def read_single_contig_fasta(filename):
    names = []
    sequences = []
    with open(filename, 'r') as f:
        line = f.readline()
        assert line.startswith(">")
        names.append(line.strip().split("\t"))
        sequence = ""
        for line in f:
            if line.startswith(">"):
                sequences.append(sequence)
                names.append(line.strip().split("\t"))
                sequence = ""
            else:
                for x in line.strip():
                    if x not in ["A", "C", "G", "T"]:
                        print("Unknown nucleotide {}".format(x), file=stderr)
                        exit(3)
                sequence += line.strip()

    sequences.append(sequence)
    assert len(names) == 1
    assert len(sequences) == 1
    return names[0], sequences[0]

def _score_sequences(seq1: str, seq2: str, match: int, mismatch: int, gapopen: int, gapextend: int):
    """Create scoring matrix for two sequences"""
    LEN1 = len(seq1)
    LEN2 = len(seq2)

    # init scoreing matrix
    H = np.zeros(
        (LEN1+1, LEN2+1), dtype=int
    )

    # gap matrix
    D = np.zeros(
        (LEN1+1, LEN2+1), dtype=int
    )

    # gap matrix
    I = np.zeros(
        (LEN1+1, LEN2+1), dtype=int
    )

    # score sequences
    for i in range(1,H.shape[0]):
        for j in range(1,H.shape[1]):
            # update gap matrix
            D[i][j] = max(
                D[i-1][j] - gapextend,
                H[i-1][j] - gapextend - gapopen
            )
            # update gap matrix
            I[i][j] = max(
                I[i][j-1] - gapextend,
                H[i][j-1] - gapextend - gapopen
            )
            # update scoring matrix
            H[i][j] += max(
                0,
                H[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else -mismatch),
                D[i][j],
                I[i][j]
            )
    
    return H

def _max_score(S: np.ndarray):
    """Find the last occurance of the maximum score"""
    _val =np.argmax(S)
    idx = np.array(np.unravel_index(((S==_val)[::-1,::-1]).argmax(), S.shape))
    return S.shape - idx - 1

def _traceback(S: np.ndarray):
    """Create a traceback path for a given score matrix"""
    # get index of maximum score
    _max_indx = _max_score(S)
    i, j = _max_indx

    # init path
    PATH = []

    while S[i][j] != 0:
        PATH.append(
            max(
                (S[i-1][j-1], "DIAG"), # move diagonaly -> first, so we move here by default
                (S[i-1][j], "UP"), # move UP
                (S[i][j-1], "LEFT"), # move LEFT
                key=lambda s: s[0])[1] # extract the score, return the P
        )
        _last = PATH[-1]
        if _last == "DIAG":
            i -= 1
            j -= 1
        if _last == "UP":
            i -= 1
        if _last == "LEFT":
            j -= 1

    return PATH

def smith_waterman(seq1: str, seq2: str, match: int, mismatch: int, gapopen: int, gapextend: int):
    max_score = 0
    alnseq1 = ""
    alnseq2 = ""

    H = _score_sequences(seq1, seq2, match, mismatch, gapopen, gapextend)
    PATH = _traceback(H)

    _max_indx = _max_score(H)
    i, j = _max_indx

    for P in PATH:
        if P == "DIAG":
            alnseq1 += seq1[i-1]
            alnseq2 += seq2[j-1]
            i-=1
            j-=1
        elif P == "LEFT":
            alnseq1 += "-"
            alnseq2 += seq2[j-1]
            j-=1
        elif P == "UP":
            alnseq1 += seq1[i-1]
            alnseq2 += "-"
            i-=1
        else:
            raise ValueError(f"Unknown P: {P}")

        if i == 0 or j == 0:
            break

    max_score = np.max(H)

    return max_score, alnseq1[::-1], alnseq2[::-1]
    

def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    name1, seq1 = read_single_contig_fasta(filename1)
    name2, seq2 = read_single_contig_fasta(filename2)

    # this function takes as input two nucleotide sequences along with
    # scores for an alignment match, mismatch, opening a new gap, and 
    # extending an existing gap. This should return the maximum alignment
    # score as well as the alignment. For examples see the testdriver script
    max_score, alnseq1, alnseq2 = smith_waterman(seq1, seq2, 
                                  match, mismatch, gapopen, gapextend)
    
    print("Maximum alignment score: {}".format(max_score))
    print("Sequence1 : {}".format(alnseq1))
    print("Sequence2 : {}".format(alnseq2))

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:],
                     "hm:x:g:e:",
                     ["help", "match=", "mismatch=", "gapopen=", "gapextend="])
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1) 

    match = 2
    mismatch = 1
    gapopen = 4
    gapextend = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        elif o in ("-m", "--match"):
            match = float(a)
        elif o in ("-x", "--mismatch"):
            mismatch = float(a)
        elif o in ("-g", "--gapopen"):
            gapopen = float(a)
        elif o in ("-e", "--gapextend"):
            gapextend = float(a)
        else:
            assert False, "unhandled option"

    if len(args) != 2:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0], args[1], match, mismatch, gapopen, gapextend)
