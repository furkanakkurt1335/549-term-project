import glob, os, json

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
glob_path = os.path.join(THIS_DIR, 'data/sequences/*/*/*.fa')
files = glob.glob(glob_path, recursive=True)

sequences = {}
for file in files:
    with open(file, 'r') as f:
        content = f.readlines()
    
    name = '-'.join(file.split('/')[-3:]).replace('.fa', '')

    sequence = ''.join([line.strip() for line in content[1:]])

    sequences[name] = sequence

with open(os.path.join(THIS_DIR, 'all_sequences.json'), 'w') as f:
    json.dump(sequences, f)
