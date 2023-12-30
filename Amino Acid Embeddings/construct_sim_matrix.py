import json, argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-aa', '--amino_acids', type=str, required=True, help='Path to amino acids file')
    parser.add_argument('-s', '--similarities', type=str, required=True, help='Path to similarities file')
    args = parser.parse_args()

    with open(args.similarities, 'r') as f:
        similarities = json.load(f)
    sim_key_l = list(similarities.keys())

    with open(args.amino_acids, 'r') as f:
        aa_d = json.load(f)
    aa_key_l = list(aa_d.keys())
    for aa in aa_key_l:
        if aa not in sim_key_l:
            del aa_d[aa]

    matrix = np.zeros((len(aa_d), len(aa_d)))
    aa_letters = list(aa_d.keys())
    for i, aa1 in enumerate(aa_letters):
        for j, aa2 in enumerate(aa_letters[i+1:]):
            matrix[i][j+i+1] = similarities[aa1][aa2]
            matrix[j+i+1][i] = similarities[aa1][aa2]
    with open('similarity_matrix.json', 'w') as f:
        json.dump({'aa\'s': aa_letters, 'matrix': matrix.tolist()}, f, indent=4)
    
if __name__ == '__main__':
    main()