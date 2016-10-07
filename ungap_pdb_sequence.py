# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:45:07 2016

@author: sshank
"""

import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


taed_id = 'taed_descendent'
pdb_id = 'pdb'
alignment = AlignIO.read('fasta/taed_aligned.fasta', 'fasta')
for seqrecord in alignment:
    if seqrecord.name == taed_id:
        taed_descendent = seqrecord.seq
    elif seqrecord.name == pdb_id:
        pdb = seqrecord.seq

# pdb is in first row
sequences = np.array([list(pdb), list(taed_descendent)])
indices = sequences[0,:] != '-'
sequences = sequences[:, indices]
convert_numpy_sequence = lambda sequence: ''.join(sequence)
pdb_string = convert_numpy_sequence(sequences[0,:])
taed_descendent_string = convert_numpy_sequence(sequences[1,:])
pdb = SeqRecord(Seq(pdb_string, IUPAC.protein), id=pdb_id)
taed_descendent = SeqRecord(Seq(taed_descendent_string, IUPAC.protein), id=taed_id)
alignment = MultipleSeqAlignment([pdb, taed_descendent])
AlignIO.write(alignment, 'fasta/taed_aligned_pdb_ungapped.fasta', 'fasta')
