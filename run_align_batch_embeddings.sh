source .venv/bin/activate
python3 alignment/align_batch_embeddings.py -sj alignment/all_sequences.json -sd alignment/data/sequences -a alignment/data/alignments --set RV913 -scr alignment/peba/peba.py -l 5
python3 alignment/align_batch_embeddings.py -sj alignment/all_sequences.json -sd alignment/data/sequences -a alignment/data/alignments --set RV911 -scr alignment/peba/peba.py -l 5
python3 alignment/align_batch_embeddings.py -sj alignment/all_sequences.json -sd alignment/data/sequences -a alignment/data/alignments --set RV912 -scr alignment/peba/peba.py -l 5