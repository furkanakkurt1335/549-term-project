# Example run can be made with the command:
#   `python3 alignment/align_batch_embeddings.py -sj alignment/all_sequences.json -sd alignment/data/sequences -a alignment/data/alignments --set RV11 -scr alignment/peba/peba.py`

import os, json, argparse, subprocess, sys, re
from proteinbert import load_pretrained_model
from proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description='Align embeddings')
    parser.add_argument('-sj', '--seq-json', help='Path to JSON file of all sequences', required=True)
    parser.add_argument('-sd', '--seq-dir', help='Path to sequences directory', required=True)
    parser.add_argument('-a', '--alignments', help='Path to alignments directory', required=True)
    parser.add_argument('--set', help='Set to align', required=True)
    parser.add_argument('-scr', '--script', help='Path to PEbA alignment script', required=True)
    return parser.parse_args()

def get_local_representation(model, input_encoder, sequence, seq_len):
    encoded_x = input_encoder.encode_X([sequence], seq_len)
    local_representations, _ = model.predict(encoded_x, batch_size = 1, verbose = 0)
    mask_to_exclude_special_tokens =  encoded_x[0][0] < 23
    local_rep = local_representations[0][mask_to_exclude_special_tokens]
    local_rep = local_rep[:, -154:-26]
    return local_rep

def main():
    THIS_DIR = os.path.dirname(os.path.abspath(__file__))
    args = get_args()
    seq_path = args.seq_json
    with open(seq_path, 'r') as f:
        sequences = json.load(f)
    seq_dir = args.seq_dir
    selected_set = args.set
    set_dir = os.path.join(seq_dir, selected_set)
    box_l = os.listdir(set_dir)

    aln_dir = args.alignments
    ref_dir = os.path.join(aln_dir, 'refs')

    python_path = sys.executable
    
    script_path = args.script

    seq_len = 3072
    pretrained_model_generator, input_encoder = load_pretrained_model()
    model = get_model_with_hidden_layers_as_outputs(pretrained_model_generator.create_model(seq_len))

    repr1_path, repr2_path = os.path.join(THIS_DIR, 'repr1.txt'), os.path.join(THIS_DIR, 'repr2.txt')
    filtered_sequences = [seq for seq in sequences if seq.startswith(args.set)]
    print('Aligning', len(filtered_sequences), 'sequences')
    alignment_save_path = os.path.join(THIS_DIR, 'alignment.msf')
    sp_pattern = re.compile(r'SP: ([0-9.]+)')
    sp_path = os.path.join(THIS_DIR, 'sps.json')
    if not os.path.exists(sp_path):
        json.dump({}, open(sp_path, 'w'))
    sps = json.load(open(sp_path, 'r'))
    if selected_set not in sps:
        sps[selected_set] = {}
    for box_t in box_l:
        filtered_sequences = [seq for seq in sequences if seq.startswith('{}-{}'.format(selected_set, box_t))]
        if box_t not in sps[selected_set]:
            sps[selected_set][box_t] = {}
        for i, sequence1_name in enumerate(filtered_sequences):
            if i == len(filtered_sequences) - 1:
                continue
            elif sequence1_name not in sps[selected_set][box_t]:
                sps[selected_set][box_t][sequence1_name] = {}
            if set(sps[selected_set][box_t][sequence1_name].keys()) == set(filtered_sequences[i+1:]):
                continue
            print('Aligning', sequence1_name)
            seq1_split = sequence1_name.split('-')
            seq1_id = seq1_split[-1]
            seq1_t = sequences[sequence1_name]
            seq1_replaced = sequence1_name.replace('-', '/')
            file1_t = seq1_replaced + '.fa'
            path1_t = os.path.join(seq_dir, file1_t)
            repr1 = get_local_representation(model, input_encoder, seq1_t, seq_len)
            np.savetxt(repr1_path, repr1)
            while 1:
                current_repr1 = np.loadtxt(repr1_path)
                if np.array_equal(current_repr1, repr1):
                    break
            for sequence2_name in filtered_sequences[i+1:]:
                if sequence2_name in sps[selected_set][box_t][sequence1_name]:
                    continue
                print('with', sequence2_name)
                seq2_split = sequence2_name.split('-')
                seq2_id = seq2_split[-1]
                seq2_t = sequences[sequence2_name]
                seq2_replaced = sequence2_name.replace('-', '/')
                file2_t = seq2_replaced + '.fa'
                path2_t = os.path.join(seq_dir, file2_t)
                repr2 = get_local_representation(model, input_encoder, seq2_t, seq_len)
                np.savetxt(repr2_path, repr2)
                while 1:
                    current_repr2 = np.loadtxt(repr2_path)
                    if np.array_equal(current_repr2, repr2):
                        break
                cmd_l = [python_path, script_path, '-f1', path1_t, '-f2', path2_t, '-e1', repr1_path, '-e2', repr2_path]
                output = subprocess.check_output(cmd_l, stderr=subprocess.STDOUT)
                output = output.decode('utf-8')
                with open(alignment_save_path, 'w') as f:
                    f.write(output)
                while 1:
                    with open(alignment_save_path, 'r') as f:
                        current_output = f.read()
                    if current_output == output:
                        break
                ref_alignment_path1 = os.path.join(ref_dir, selected_set, box_t, '{}-{}.msf'.format(seq1_id, seq2_id))
                ref_alignment_path2 = os.path.join(ref_dir, selected_set, box_t, '{}-{}.msf'.format(seq2_id, seq1_id))
                if os.path.exists(ref_alignment_path1):
                    ref_alignment_path = ref_alignment_path1
                elif os.path.exists(ref_alignment_path2):
                    ref_alignment_path = ref_alignment_path2
                score_cmd_l = [python_path, os.path.join(THIS_DIR, 'evaluation_scripts/compute_score.py'), '-align1', alignment_save_path, '-align2', ref_alignment_path, '-score', 'sp']
                score_output = subprocess.check_output(score_cmd_l, stderr=subprocess.STDOUT)
                score_output = score_output.decode('utf-8')
                sp_match = sp_pattern.search(score_output)
                if sp_match:
                    sp = float(sp_match.group(1))
                    sps[selected_set][box_t][sequence1_name][sequence2_name] = sp
                    with open(sp_path, 'w') as f:
                        json.dump(sps, f, indent=4)
                else:
                    print('Error: no SP score found', sequence1_name, sequence2_name)
            print('Finished aligning', sequence1_name)
            print('Remaining:', len(filtered_sequences) - i - 1, 'sequences')
        print('Finished aligning box', box_t)
        print('Remaining:', len(box_l) - box_l.index(box_t) - 1, 'boxes')
    print('Finished aligning set', selected_set)

if __name__ == '__main__':
    main()
