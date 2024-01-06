from Bio.Align import substitution_matrices
from Bio import Align, SeqIO
import json, argparse, os

parser = argparse.ArgumentParser(description='Align two sequences')
parser.add_argument('-m', '--matrix', help='our matrix to use for alignment', required=True)
args = vars(parser.parse_args())

# Load the matrix
with open(args["matrix"]) as f:
    our_m = json.load(f)
aa_l = ''.join(our_m["aa's"])
m = our_m["matrix"]

blosum62_m = substitution_matrices.load("BLOSUM62") # Load the BLOSUM62 matrix
pam250_m = substitution_matrices.load("PAM250") # Load the PAM250 matrix
our_matrix = blosum62_m.copy()

for i in range(len(aa_l)):
    for j in range(len(aa_l)):
        try:
            # Multiply with 10 to bring cosine sims to similar scale with blosum (as done in PEbA paper)
            our_matrix[aa_l[i]][aa_l[j]] = m[i][j]*10
        except:
            pass


aligner = Align.PairwiseAligner()
def align_sequences_biopython(seq1, seq2, penalty_params, matrix, mode = 'global'):
    aligner.substitution_matrix = matrix
    if penalty_params["penalty_type"] == "linear":
        aligner.open_gap_score = penalty_params["gap_penalty"]
        aligner.extend_gap_score = penalty_params["gap_penalty"]
    elif penalty_params["penalty_type"] == "affine":
        aligner.open_gap_score = penalty_params["gap_penalty"]
        aligner.extend_gap_score = penalty_params["gap_extension_penalty"]

    aligner.mode = mode
    return aligner.align(seq1, seq2)

#seq1 = "PANAMABANANA"
#seq2 = "PANKLBANA"
fa_dir = 'data/sequences/RV11/BB11001'
seq1_file = os.path.join(fa_dir, '1aab_.fa')
seq1 = SeqIO.read(seq1_file, "fasta")
seq2_file = os.path.join(fa_dir, '1j46_A.fa')
seq2 = SeqIO.read(seq2_file, "fasta")

print("Linear penalty")
alignment = align_sequences_biopython(seq1, seq2, {"penalty_type": "linear", "gap_penalty": -8}, blosum62_m)
print("Score: ", alignment.score)
print("BLOSUM62")
print(alignment[0][0])
print(alignment[0][1])
print('----------------')
alignment = align_sequences_biopython(seq1, seq2, {"penalty_type": "linear", "gap_penalty": -8}, pam250_m)
print("Score: ", alignment.score)
print("PAM250")
print(alignment[0][0])
print(alignment[0][1])
print('----------------')
alignment = align_sequences_biopython(seq1, seq2, {"penalty_type": "linear", "gap_penalty": -8}, our_matrix)
print("Score: ", alignment.score)
print("Our matrix")
print(alignment[0][0])
print(alignment[0][1])

print("\nAffine Penalty")
alignment = align_sequences_biopython(seq1, seq2, {"penalty_type": "affine", "gap_penalty": -10, "gap_extension_penalty": -1}, blosum62_m)
print("Score: ", alignment.score)
print("BLOSUM62")
print(alignment[0][0])
print(alignment[0][1])
print('----------------')
alignment = align_sequences_biopython(seq1, seq2, {"penalty_type": "affine", "gap_penalty": -10, "gap_extension_penalty": -1}, pam250_m)
print("Score: ", alignment.score)
print("PAM250")
print(alignment[0][0])
print(alignment[0][1])
print('----------------')
alignment = align_sequences_biopython(seq1, seq2, {"penalty_type": "affine", "gap_penalty": -10,"gap_extension_penalty": -1}, our_matrix)
print("Score: ", alignment.score)
print("Our matrix")
print(alignment[0][0])
print(alignment[0][1])

# TODO: Compute Sum of pairs scores according to the reference alignments
