# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:30:11 2016

@author: sshank
"""

import re
import os
import subprocess
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC


#DIRECTORY_ROOT = '/var/www/TAED/'
DIRECTORY_ROOT = '/home/sshank/Documents/TAED/'
DATABASE_ROOT = DIRECTORY_ROOT+'DATABASE99/'
TAEDPV_ROOT = DIRECTORY_ROOT+'taed-pv/'


def parse_rst_file(rst_file, ancestral_index, descendent_index):
    directory = rst_file['directory']
    file_number = rst_file['file_number']
    paml_subtree = rst_file['paml_subtree']
    input_filename = str(file_number)+'_'+str(paml_subtree)+'.RST'
    input_path = '%s%d/%s' % (DATABASE_ROOT, directory, input_filename)
    ancestral_string = ''
    descendent_string = ''
    p = re.compile('\s+\d+\s+\d+\s+[ATCG-]+.+')
    with open(input_path, 'r') as input_file:
        for line in input_file:
            if p.match(line):
                split_line = line.split()
                # This requires careful inspection of RST files before and after :
                index = split_line.index(':')
                first_half = split_line[:index]
                second_half = split_line[index+1:]
                first_codons = first_half[2::2]
                second_codons = second_half[::5]
                all_codons = first_codons+second_codons
                ancestral_codon = all_codons[ancestral_index-1]
                ancestral_string += ancestral_codon
                descendent_codon = all_codons[descendent_index-1]
                descendent_string += descendent_codon
            if line[:3] == 'Sum':
                break
    ancestor = Seq(ancestral_string)
    descendent = Seq(descendent_string)
    return ancestor, descendent


def process_rst_for_alignment(famMapID, rst_file, ancestral_index, descendent_index, pdb_id):
    _, descendent = parse_rst_file(rst_file, ancestral_index, descendent_index)
    pdb_filename = '%sfasta/%s.fasta' % (TAEDPV_ROOT, pdb_id)
    pdb_record = SeqIO.read(pdb_filename, 'fasta')
    descendent_record = SeqRecord(descendent.ungap('-').translate(), id=famMapID)
    output_path = '%sfasta/%s_unaligned.fasta' % (TAEDPV_ROOT, famMapID)
    SeqIO.write([descendent_record, pdb_record], output_path, 'fasta')


def align(famMapID):
    input_path = '%sfasta/%s_unaligned.fasta' % (TAEDPV_ROOT, famMapID)
    output_path = '%sfasta/%s_gaps.fasta' % (TAEDPV_ROOT, famMapID)
    alignment_command = 'mafft %s > %s' % (input_path, output_path)
    p = subprocess.Popen(alignment_command, shell=True)
    os.waitpid(p.pid, 0)


def prepare_for_visualizer(rst_file, ancestral_index, descendent_index, famMapID, pdb_id):
    # File paths
    gapped_input_path = '%sfasta/%s_gaps.fasta' % (TAEDPV_ROOT, famMapID)
    fasta_output_path = '%sfasta/%s.fasta' % (TAEDPV_ROOT, famMapID)
    info_output_path = '%sinfo/%s.txt' % (TAEDPV_ROOT, famMapID)
    ancestor, descendent = parse_rst_file(rst_file, ancestral_index, descendent_index)
    
    # Get information about changes
    descendent_annotations = []
    descendent_changes = []
    sequence_length = len(ancestor)
    for k in range(0, sequence_length, 3):
        descendent_codon = str(descendent[k:k+3])
        ancestral_codon = str(ancestor[k:k+3])
        if descendent_codon != '---':
            descendent_amino_acid = Seq(descendent_codon).translate()
            if descendent_codon == ancestral_codon or ancestral_codon == '---':
                # No change or missing information
                descendent_annotations.append(0)
                descendent_changes.append('-')
            else:
                ancestral_amino_acid = Seq(ancestral_codon).translate()
                if descendent_amino_acid == ancestral_amino_acid:
                    # Synonymous change
                    descendent_annotations.append(1)
                    change = ancestral_codon + '->' + descendent_codon
                    descendent_changes.append(change)
                else:
                    # Nonsynonymous change
                    descendent_annotations.append(2)
                    change = str(ancestral_amino_acid) + '->' + str(descendent_amino_acid)
                    descendent_changes.append(change)

    # Annotate pdb/descendent alignment
    pdb_annotations = []
    pdb_changes = []
    pdb_indices = []
    alignment = AlignIO.read(gapped_input_path, 'fasta')
    assert alignment[1].id == pdb_id
    d_index = 0
    p_index = 0
    
    for k in range(alignment.get_alignment_length()):
        descendent_amino_acid, pdb_amino_acid = alignment[:, k]
        if pdb_amino_acid != '-' and descendent_amino_acid != '-':
            # There is a chance that something happened... append and increment both
            pdb_annotations.append(descendent_annotations[d_index])
            pdb_changes.append(descendent_changes[d_index])
            p_index += 1
            d_index += 1
            if pdb_annotations[-1] == 2:
                # Order ensures indices are 1-based, for end user
                pdb_indices.append(p_index)
        else:
            if pdb_amino_acid != '-':
                pdb_annotations.append(0)
                pdb_changes.append('-')
                p_index += 1
            if descendent_amino_acid != '-':
                d_index += 1

    # Ungap sequence
    for seqrecord in alignment:
        if seqrecord.name == pdb_id:
            taed_descendent = seqrecord.seq
        else:
            pdb = seqrecord.seq
    
    # pdb is in first row
    sequences = np.array([list(pdb), list(taed_descendent)])
    indices = sequences[0,:] != '-'
    sequences = sequences[:, indices]
    convert_numpy_sequence = lambda sequence: ''.join(sequence)
    pdb_string = convert_numpy_sequence(sequences[0,:])
    taed_descendent_string = convert_numpy_sequence(sequences[1,:])
    pdb = SeqRecord(Seq(pdb_string, IUPAC.protein), id=pdb_id)
    taed_descendent = SeqRecord(Seq(taed_descendent_string, IUPAC.protein), id=famMapID)
    alignment = MultipleSeqAlignment([pdb, taed_descendent])
    AlignIO.write(alignment, fasta_output_path, 'fasta')
    
    # Write output
    with open(info_output_path, 'w') as info_file:
        info_file.write(','.join([str(i) for i in pdb_annotations]))
        info_file.write('\n')
        info_file.write(','.join([str(i) for i in pdb_changes]))
        info_file.write('\n')
        info_file.write(','.join([str(i) for i in pdb_indices]))
        
        
    

if __name__ == '__main__':
    rst_file={'directory':111,'file_number':10005,'paml_subtree':1}
    ancestral_index, descendent_index = 32, 7
    parse_rst_file(rst_file, ancestral_index, descendent_index)
    famMapID = '10005_3778'
    pdb_id = '4GA7_A'
    process_rst_for_alignment(famMapID,
                              rst_file,
                              ancestral_index,
                              descendent_index,
                              pdb_id)
    align(famMapID)
    prepare_for_visualizer(rst_file, ancestral_index, descendent_index, famMapID, pdb_id)