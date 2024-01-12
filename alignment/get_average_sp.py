import argparse, json

def get_args():
    parser = argparse.ArgumentParser(description='Align embeddings')
    parser.add_argument('-i', '--input', help='Path to input file', required=True)
    parser.add_argument('-s', '--set', help='Selected set', required=True)
    return parser.parse_args()

def main():
    args = get_args()
    input_path = args.input
    with open(input_path, 'r') as f:
        sps = json.load(f)
    
    sum_t, count_t = 0, 0
    set_t = args.set
    for box_t in sps[set_t]:
        for seq_t in sps[set_t][box_t]:
            for seq2_t in sps[set_t][box_t][seq_t]:
                sum_t += sps[set_t][box_t][seq_t][seq2_t]
                count_t += 1
    print(sum_t / count_t)

if __name__ == '__main__':
    main()