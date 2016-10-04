wget -O pdbs/4PVP.pdb https://files.rcsb.org/download/4PVP.pdb
wget -O paml/15640_8.RST https://liberles.cst.temple.edu/TAED/DATABASE99/102/15640_8.RST
python parse_rst.py -f paml/15640_8.RST > paml/15640_8.RST.parsed
python -m http.server