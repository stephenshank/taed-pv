# TAED Protein Viewer

Visualizing proteins from TAED, The Adaptive Evolution Database. Work in progress.

To run, enter:

```
git clone https://github.com/stephenshank/taed-pv
cd taed-pv
wget -O pdbs/4PVP.pdb https://files.rcsb.org/download/4PVP.pdb
python -m http.server
```

and visit `http://0.0.0.0:8000/TAEDProteinViewer.html`.