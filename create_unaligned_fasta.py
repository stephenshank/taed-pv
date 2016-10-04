# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:27:48 2016

@author: sshank
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser
from argparse import ArgumentParser


parser = ArgumentParser()
rst_help = 'Path to parsed RST file (created with parse_rst.py).'
parser.add_argument('-r', '--rst', metavar='RST', help=rst_help, dest='rst')
pdb_help = 'Path to PDB file.'
parser.add_argument('-p', '--pdb', metavar='PDB', help=pdb_help, dest='pdb')
output_help = 'Path to output fasta file (not aligned).'
parser.add_argument('-o', '--output', metavar='OUTPUT', help=output_help, dest='output')
args = parser.parse_args()
rst_filename = args.rst
pdb_filename = args.pdb
output_filename = args.output

descendent_sequence = ''
with open(rst_filename, 'r') as file:
    for line in file:
        split = line.split()
        # This is fine tuned to the current RST file... MUST CHANGE!!
        descendent_codon = split[16]
        if descendent_codon != '---':
            descendent_sequence += descendent_codon
taed_descendent = SeqRecord(Seq(descendent_sequence).translate(), id='taed_descendent')

pdb_parser = PDBParser()
structure = pdb_parser.get_structure('example', pdb_filename)
model = structure[0]
chainA, chainB = model.get_chains()
sequenceA = []
for residue in chainA.get_residues():
    sequenceA.append(residue.get_resname())
wrong_residues =['HOH',' NA', 'MLI']
desired_residues = [residue for residue in sequenceA if residue not in wrong_residues]
pdb_verbose_sequence = ','.join(desired_residues)

conversion_dict = {'PHE':'F', 'LEU':'L', 'ILE':'I', 'MET':'M',
                   'VAL':'V', 'SER':'S', 'PRO':'P', 'THR':'T',
                   'ALA':'A', 'TYR':'Y', 'HIS':'H', 'GLN':'Q',
                   'ASN':'N', 'LYS':'K', 'ASP':'D', 'GLU':'E',
                   'CYS':'C', 'TRP':'W', 'ARG':'R', 'GLY':'G'}

def convert_verbose(sequence_str):
    sequence_list = sequence_str.split(',')
    concise_sequence_list = [conversion_dict[residue] for residue in sequence_list]
    concise_sequence = ''.join(concise_sequence_list)
    return concise_sequence

pdb_sequence = convert_verbose(pdb_verbose_sequence)
pdb_record = SeqIO.SeqRecord(Seq(pdb_sequence), id='pdb')
SeqIO.write([taed_descendent, pdb_record], output_filename, 'fasta')
