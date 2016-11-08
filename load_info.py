# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 08:38:37 2016

@author: sshank
"""

import pymysql
import os

connection = pymysql.connect(host='liberlesBox.bio.temple.edu',
                             user=os.environ.get('TAEDUSER'),
                             password=os.environ.get('TAEDPASS'),
                             db='TAED2',
                             charset='utf8mb4',
                             cursorclass=pymysql.cursors.DictCursor)

with connection.cursor() as cursor:
    # Read a single record
    sql = '''SELECT file.taedFileNumber, file.familyName, fam.dndsValue,
                    fam.mappedBranchEnd, fam.mappedBranchStart, file.pdb,
                    fam.famMapID
             FROM taedfile3 as file
             JOIN famMap4 as fam
             ON file.taedFileNumber=fam.taedFileNumber;'''
    cursor.execute(sql)
    results = cursor.fetchall()
    for result in results:
        try:
            with open('info/%s.txt' % result['famMapID'],'r') as input_file:
                file_text = input_file.read()
            annotations, changes, indices = file_text.split('\n')
            sql = "INSERT INTO proteinViewer VALUES ('%s','%s','%s',%f,'%s','%s','%s','%s','%s');"
            info = (result['famMapID'],
                    result['familyName'],
                    result['pdb'],
                    result['dndsValue'],
                    result['mappedBranchStart'],
                    result['mappedBranchEnd'],
                    annotations,
                    changes,
                    indices)
            query = sql % info
            cursor.execute(query)
            print(query)
            print('Inserted famMapID %s' % result['famMapID'])
        except:
            print('Had a problem with famMapID %s' % result['famMapID'])
connection.commit()
connection.close()
