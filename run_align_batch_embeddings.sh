source .venv/bin/activate
python3 alignment/align_batch_embeddings.py -sj alignment/all_sequences.json -sd alignment/data/sequences -a alignment/data/alignments --set RV913 -scr alignment/peba/peba.py