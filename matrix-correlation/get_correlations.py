from Bio.Align import substitution_matrices
import json, argparse, random
from scipy.stats import pearsonr

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

# get pearson correlation
blosum62_l = []
our_l = []
pam250_l = []
random_l = []
for i in range(len(aa_l)):
    for j in range(len(aa_l)):
        try:
            blosum62_m[aa_l[i], aa_l[j]]
            pam250_m[aa_l[i], aa_l[j]]
            blosum62_l.append(blosum62_m[aa_l[i], aa_l[j]])
            pam250_l.append(pam250_m[aa_l[i], aa_l[j]])
            our_l.append(m[i][j])
            random_l.append(random.random())
        except:
            pass
print("Pearson correlation between our matrix and BLOSUM62:", '%.3f' % pearsonr(blosum62_l, our_l)[0])
print("Pearson correlation between our matrix and PAM250:", '%.3f' % pearsonr(pam250_l, our_l)[0])
print("Pearson correlation between BLOSUM62 and PAM250:", '%.3f' % pearsonr(blosum62_l, pam250_l)[0])
print("Pearson correlation between BLOSUM62 and random:", '%.3f' % pearsonr(blosum62_l, random_l)[0])
# Pearson correlation between our matrix and BLOSUM62: 0.862
# Pearson correlation between our matrix and PAM250: 0.690
# Pearson correlation between BLOSUM62 and PAM250: 0.840
# Pearson correlation between BLOSUM62 and random: -0.006