# taken from notebook by Melikşah in this dir ; updated by Furkan
import json, os
from proteinbert import load_pretrained_model
from proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs

seq_len = 3072
pretrained_model_generator, input_encoder = load_pretrained_model()
model = get_model_with_hidden_layers_as_outputs(pretrained_model_generator.create_model(seq_len))

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
seq_path = os.path.join(THIS_DIR, '../alignment/all_sequences.json')

with open(seq_path, 'r') as f:
    sequences = json.load(f)

for sequence in sequences:
    seq_t = sequences[sequence]
    encoded_x = input_encoder.encode_X([seq_t], seq_len)
    local_representations, _ = model.predict(encoded_x, batch_size = 1, verbose = 0)

    # SPECIAL TOKENS: 22: UNK, 23: BOS, 24: EOS, 25: PAD
    mask_to_exclude_special_tokens =  encoded_x[0][0] < 23 # exclude 23: BOS, 24: EOS, 25: PAD
    local_rep = local_representations[0][mask_to_exclude_special_tokens]
