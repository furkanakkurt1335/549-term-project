import json, argparse
from sklearn.metrics.pairwise import cosine_similarity

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--embeddings', type=str, required=True, help='Path to embeddings file')
    parser.add_argument('-aa', '--amino_acids', type=str, required=True, help='Path to amino acids file')
    args = parser.parse_args()

    with open(args.embeddings, 'r') as f:
        embeddings = json.load(f)

    with open(args.amino_acids, 'r') as f:
        aa_d = json.load(f)

    amino_acids = list(aa_d.keys())

    filtered_amino_acids = []
    for aa in amino_acids:
        if aa in embeddings:
            filtered_amino_acids.append(aa)
        else:
            print(aa)

    amino_acids = filtered_amino_acids

    similarity_d = {}
    for i, aa1 in enumerate(amino_acids):
        if aa1 not in similarity_d:
            similarity_d[aa1] = {}
        for j, aa2 in enumerate(amino_acids[i:]):
            if aa2 not in similarity_d:
                similarity_d[aa2] = {}
            sim_t = cosine_similarity([embeddings[aa1]], [embeddings[aa2]])[0][0]
            similarity_d[aa1][aa2] = sim_t
            similarity_d[aa2][aa1] = sim_t

    with open('similarities.json', 'w') as f:
        json.dump(similarity_d, f, indent=4)

if __name__ == '__main__':
    main()