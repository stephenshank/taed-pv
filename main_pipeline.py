# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:15:54 2016

@author: sshank
"""

import os
import pymysql
from pipeline_tools import *

connection = pymysql.connect(host='127.0.0.1',
                             user=os.environ.get('TAEDUSER'),
                             password=os.environ.get('TAEDPASS'),
                             db='TAED2',
                             charset='utf8mb4',
                             cursorclass=pymysql.cursors.DictCursor,
                             unix_socket='/var/run/mysqld/mysqld.sock')

with connection.cursor() as cursor:
    # Read a single record
    sql = '''select file.taedFileNumber,
                    file.Directory,
                    file.pdb,
                    fam.pamlSubtree,
                    fam.pamlNodes,
                    fam.famMapID
             from taedfile2 as file 
             join famMap3 as fam
             on file.taedFileNumber=fam.taedFileNumber limit 50;'''
    cursor.execute(sql)
    result = cursor.fetchone()
#    while result:
#        result = cursor.fetchone()
    
connection.close()