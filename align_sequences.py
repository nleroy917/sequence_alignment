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
from itertools import groupby

# set flags
STOP, DIAG, UP, LEFT= range(4)

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

    # init scoring matrix
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

    # traceback matrix
    T = np.zeros(
        (LEN1+1, LEN2+1), dtype=int
    )

    # score sequences
    for i in range(1,H.shape[0]):
        for j in range(1,H.shape[1]):
        
            # update deletion gap matrix
            D[i][j] = max(
                0,
                D[i-1][j] - gapextend,
                H[i-1][j] - gapopen
            )

            # update insertion matrix
            I[i][j] = max(
                0,
                I[i][j-1] - gapextend,
                H[i][j-1] - gapopen
            )

            # update scoring matrix
            MATCH = H[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else -mismatch)
            DEL = D[i][j]
            INSR = I[i][j]

            H[i][j] += max(
                0,
                MATCH,
                DEL,
                INSR,
            )

            # update traceback matrix
            if H[i][j] == 0:
                T[i][j] = STOP
            elif H[i][j] == MATCH:
                T[i][j] = DIAG
            elif H[i][j] == DEL:
                T[i][j] = UP
            elif H[i][j] == INSR:
                T[i][j] = LEFT
            else:
                T[i][j] = -1

    return H, T

def _max_score_indx(S: np.ndarray):
    """Find the *last* occurance of the maximum score"""
    rows, cols = np.where(S == S.max())
    return rows[-1], cols[-1]

def smith_waterman(seq1: str, seq2: str, match: int, mismatch: int, gapopen: int, gapextend: int):
    max_score = 0
    alnseq1 = ""
    alnseq2 = ""

    # create scoring and traceback matrix
    H, T = _score_sequences(seq1, seq2, match, mismatch, gapopen, gapextend)

    # get index of maximum value
    _max_indx = _max_score_indx(H)
    i, j = _max_indx

    # traceback
    while T[i][j] != STOP:
        _dir = T[i][j]
        if _dir == DIAG:
            j-=1
            i-=1
            alnseq1 += seq1[i]
            alnseq2 += seq2[j]
            
        elif _dir == LEFT:
            j-=1
            alnseq1 += "-"
            alnseq2 += seq2[j]
        elif _dir == UP:
            i-=1
            alnseq1 += seq1[i]
            alnseq2 += "-"
        else:
            raise ValueError(f"Unknown direction: {_dir}")

        if i == 0 or j == 0:
            break

    max_score = np.max(H)

    return max_score, alnseq1[::-1], alnseq2[::-1]
    

def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    _, seq1 = read_single_contig_fasta(filename1)
    _, seq2 = read_single_contig_fasta(filename2)

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
