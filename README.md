# Transcriptologs
Transcriptologs is a python software for the detection of orthology relationships based on transcript sequences of two different species.
Transcriptologs takes as input files the results of two tBLASTx all-versus-all output files (the output format of the BLAST+ software must be in the extended alignment view, obtained with the outfmt "0" option). As an example, if you want to predict the orthologs between human and mouse, you should launch two tBLASTx: one considering all the human transcript sequences as queries and all the mouse transcript sequences as subjects, and the other one considering all the mouse transcript sequences as queries and all the human transcript sequences as subjects.
Transcriptologs gives as output a file of orthology relationships between all the genes of the two considered species, predicted using an in-house techique of alignment extension and using a Bidirectional Best Hit approach.
To use the software, python version 3.3 or above is required. 

Example of the command-line use: transcriptologs.py -i1 input_1_file -i2 input_2_file -o outputfile
