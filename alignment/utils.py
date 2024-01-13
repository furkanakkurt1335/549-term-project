import numpy as np
from Bio.Align import substitution_matrices
from Bio import Align, SeqIO
import json, argparse

# write_msf is taken from: https://github.com/mgtools/PEbA/blob/master/scripts/utility.py
def write_msf(seq1: str, seq2: str, id1: str, id2: str, method: str,
               gopen: float, gext: float, path: str, beg: list, end: list):
    """Writes alignment to file in msf format

    :param seq1: first aligned sequence
    :param seq2: second aligned sequence
    :param id1: first sequence id
    :param id2: second sequence id
    :param method: scoring method used
    :param gopen: gap penalty for opening a new gap
    :param gext: gap penalty for extending a gap
    :param path: directory to write alignment to
    :param list: beginning positions of seqs in alignment
    :param list: end positions of seqs in alignment
    """

    fpath = f'{path}/{id1}-{id2}.msf'

    # Add space every 10 characters
    seq1 = [seq1[i:i+10] for i in range(0, len(seq1), 10)]
    seq1 = ' '.join(seq1)
    seq2 = [seq2[i:i+10] for i in range(0, len(seq2), 10)]
    seq2 = ' '.join(seq2)

    # Split sequences every 50 characters
    seq1_split = [seq1[i:i+55] for i in range(0, len(seq1), 55)]
    seq2_split = [seq2[i:i+55] for i in range(0, len(seq2), 55)]

    # Add extra spaces to either id if they are not the same length
    if len(id1) != len(id2):
        if len(id1) > len(id2):
            id2 = id2 + ' ' * (len(id1) - len(id2))
        else:
            id1 = id1 + ' ' * (len(id2) - len(id1))

    # Put alignment in string format
    length = len(seq1)
    alignment = 'PileUp\n\n\n\n'
    alignment += f'   MSF: {length}  Type: P  Method: {method}  Gopen: {gopen}  Gext: {gext}\n\n'
    alignment += f' Name: {id1} oo  Len:  {length}  Start/End:  {beg[0]},{end[0]}\n'
    alignment += f' Name: {id2} oo  Len:  {length}  Start/End:  {beg[1]},{end[1]}\n\n//\n\n\n\n'
    for i in range(len(seq1_split)): # pylint: disable=C0200
        alignment += f'{id1}      {seq1_split[i]}\n'
        alignment += f'{id2}      {seq2_split[i]}\n\n'

    # If no path is determined then print to console, otherwise write to file
    if path == 'n':
        print(alignment)
    else:
        with open(fpath, 'w', encoding='utf8') as file:
            file.write(alignment)

def get_scaled_similarity_matrix(unscaled_matrix_path = "similarity_matrix.json",factor = 10.0):
    with open(unscaled_matrix_path) as f:
        our_m = json.load(f)
    aa_l = ''.join(our_m["aa's"])
    m = our_m["matrix"]

    blosum62_m = substitution_matrices.load("BLOSUM62") # Load the BLOSUM62 matrix
    embedding_similarity_matrix = blosum62_m.copy()

    for i in range(len(aa_l)):
        for j in range(len(aa_l)):
            try:
                # Multiply with 10 to bring cosine sims to similar scale with blosum (as done in PEbA paper)
                embedding_similarity_matrix[aa_l[i]][aa_l[j]] = m[i][j]*factor
            except:
                pass

    return embedding_similarity_matrix


def align_sequences_biopython(seq1, seq2, penalty_params, matrix, mode = 'global'):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    if penalty_params["penalty_type"] == "linear":
        aligner.open_gap_score = penalty_params["gap_penalty"]
        aligner.extend_gap_score = penalty_params["gap_penalty"]
    elif penalty_params["penalty_type"] == "affine":
        aligner.open_gap_score = penalty_params["gap_penalty"]
        aligner.extend_gap_score = penalty_params["gap_extension_penalty"]

    aligner.mode = mode
    return aligner.align(seq1, seq2)
