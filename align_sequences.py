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

# set flags
STOP, DIAG, UP, LEFT = range(4)

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
        (LEN1+1, LEN2+1), dtype=float
    )

    # traceback matrix
    T = np.zeros(
        (LEN1+1, LEN2+1), dtype=float
    )

    # score sequences
    for i in range(1,H.shape[0]):
        for j in range(1,H.shape[1]):

            # calculate diagonal score
            if seq1[i-1] == seq2[j-1]:
                MATCHSCORE = H[i-1][j-1] +match
            else:
                MATCHSCORE = H[i-1][j-1] -mismatch
            
            # calculate up score
            if T[i-1][j] == UP:
                num_ups = 0
                k = i-1
                while T[k][j] == UP:
                    num_ups += 1
                    k -= 1
                UPSCORE = H[i-num_ups-1][j]  -gapopen-(num_ups+1)*(gapextend)
            else:
                UPSCORE = H[i-1][j] -gapopen
            
            # calculate left score
            if T[i][j-1] == LEFT:
                num_lefts = 0
                l = j-1
                while T[i][l] == LEFT:
                    num_lefts += 1
                    l -= 1
                LEFTSCORE = H[i][j-num_lefts-1] -gapopen-(num_lefts+1)*(gapextend)
            else:
                LEFTSCORE = H[i][j-1] -gapopen

            # get the max
            H[i][j] += max(
                0,
                MATCHSCORE,
                UPSCORE,
                LEFTSCORE,
            )

            # update traceback matrix
            if H[i][j] == MATCHSCORE:
                T[i][j] = DIAG

            elif H[i][j] == UPSCORE:
                T[i][j] = UP

            elif H[i][j] == LEFTSCORE:
                T[i][j] = LEFT

            else:
                T[i][j] = STOP

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
    while T[i][j]:
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
