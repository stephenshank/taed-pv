# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:30:11 2016

@author: sshank
"""

import re
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


DIRECTORY_ROOT = '/var/www/TAED/'
TAEDPV_ROOT = DIRECTORY_ROOT+'taed-pv/'


def parse_rst_file(directory, file_number, paml_subtree):
    input_filename = str(file_number)+'_'+str(paml_subtree)+'.RST'
    input_path = DIRECTORY_ROOT+str(directory)+'/'+input_filename
    output_filename = input_filename+'.parsed'
    output_path = DIRECTORY_ROOT+str(directory)+'/'+output_filename
    p = re.compile('\s+\d+\s+\d+\s+[ATCG-]+.+')
    with open(input_path, 'r') as input_file, open(output_path, 'w') as output_file:
        for line in input_file:
            if p.match(line):
                output_file.write(line)
            if line[:3] == 'Sum':
                break


#def extract_sequences(parsed_rst_filename, ancestral_index, descendent_index):
#    ancestral_string = ''
#    descendent_string = ''
#    with open(parsed_rst_filename, 'r') as input_file:
#        for line in input_file:
#            split_line = line.split()
#            # This requires careful inspection of RST files before and after :
#            index = split_line.index(':')
#            first_half = split_line[:index]
#            second_half = split_line[index+1:]
#            first_codons = first_half[2::2]
#            second_codons = second_half[::5]
#            all_codons = first_codons+second_codons
#            
#            # Append codons only if they are not blank
#            ancestral_codon = all_codons[ancestral_index]
#            if ancestral_codon != '---':
#                ancestral_string += ancestral_codon
#            descendent_codon = all_codons[descendent_index]
#            if descendent_codon != '---':
#                descendent_string += descendent_codon
#    
#    ancestor = Seq(ancestral_string).translate()
#    descendent = Seq(descendent_string).translate()
#    return ancestor, descendent
#
#
#def process_rst_for_alignment(parsed_rst_filename, nodes, pdb_filename):
#    ancestral_index = nodes[0]
#    descendent_index = nodes[1]
#    ancestor, descendent = extract_sequences(parsed_rst_filename,
#                                             ancestral_index,
#                                             descendent_index)
#    pdb = SeqIO.read(pdb_filename, fasta)
#    descendent_record = SeqRecord(descendent)
#    pdb_record = SeqIO.SeqRecord(Seq(pdb_sequence), id='pdb')
#    SeqIO.write([taed_descendent, pdb_record], output_filename, 'fasta')
