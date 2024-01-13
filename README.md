# 549-term-project
Two main scripts that are used for obtaining pairwise alignments are write_alignments.py and write_alignments_peba.py in the [alignment folder](https://github.com/furkanakkurt1335/549-term-project/tree/main/alignment). The first one is used to obtain PWAs with BLOSUM and custom amino-acid similarity matrix with the PairwiseAligner of Biopython, whereas the second one use the implementation of [PEbA](https://github.com/mgtools/PEbA) with the embeddings we extracted from ProteinBERT. After obtaining the alignments as .msf files, we use the [evaluation scripts](https://github.com/furkanakkurt1335/549-term-project/tree/main/alignment/evaluation_scripts) to obtain evaluation metric scores on .log files, which are then parsed to obtain the average scores.

The resulting alignments and the log files can be found in the [data folder](https://github.com/furkanakkurt1335/549-term-project/tree/main/alignment/data).
<!-- Fill for repo presentation -->
Base Reference Papers: 
1. [ProteinBERT: a universal deep-learning model of protein sequence and function](https://academic.oup.com/bioinformatics/article/38/8/2102/6502274)
2. [Protein Embedding based Alignment
](https://www.authorea.com/users/623259/articles/646069-protein-embedding-based-alignment)
3. [BAliBASE 3.0: latest developments of the multiple sequence alignment benchmark](https://pubmed.ncbi.nlm.nih.gov/16044462/)
