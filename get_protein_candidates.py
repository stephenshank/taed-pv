#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:37:56 2016

@author: sshank
"""

import os
import csv

import pymysql


db_credentials = { 
    'host': os.getenv('TAEDIP'),
    'user': os.getenv('TAEDUSER'),
    'password': os.getenv('TAEDPASS'),
    'db': 'TAED2'
}
connection = pymysql.connect(**db_credentials)
cursor = connection.cursor(pymysql.cursors.DictCursor)

all_queries = False
current = True
fammap = False
taedfile = False

# FamMap data
if all_queries or fammap:
    query = '''
        SELECT *
        FROM famMap4;
    '''
    cursor.execute(query)
    
    with open('csv/fammap_data.csv', 'w') as csvfile:
        fieldnames = ['dndsValue',
                      'pamlSubTree',
                      'taedFileNumber',
                      'mappedBranchEnd',
                      'mappedBranches',
                      'famMapID',
                      'pamlNodes',
                      'mappedBranchStart']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for result in cursor.fetchall():
            writer.writerow(result)

# TAEDFile data
if all_queries or taedfile:
    query = '''
        SELECT *
        FROM taedfile3;
    '''
    cursor.execute(query)
    
    with open('csv/taedfile_data.csv', 'w') as csvfile:
        fieldnames = ['interleafed',
                      'reconciledTree',
                      'mapPDBLine',
                      'baseDirectory',
                      'familyName',
                      'pValue',
                      'keggNames',
                      'nhxRooted',
                      'keggPathways',
                      'speciesTreeMapLocation',
                      'Directory',
                      'nhx',
                      'paml1RatioDNDS',
                      'pamlRooted',
                      'alignmentFile',
                      'subtrees',
                      'namedTreeFile',
                      'significant',
                      'root',
                      'pdb',
                      'taedFileNumber',
                      'complete',
                      'geneTree',
                      'positiveRatio',
                      'paml1RatioRooted',
                      'paml',
                      'famMapLine']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for result in cursor.fetchall():
            writer.writerow(result)

connection.close()