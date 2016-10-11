# TAED Protein Viewer

Visualizing proteins from TAED, [The Adaptive Evolution Database](https://liberles.cst.temple.edu/TAED/). Work in progress.

![Snapshot of TAED Protein Viewer](./images/snapshot.png =500x250 "Snapshot of TAED Protein Viewer")

To run (which involves scraping some data that gets ran through a pipeline), enter:

```
git clone https://github.com/stephenshank/taed-pv
cd taed-pv
source sample_pipeline.sh
python -m http.server
```

and visit `http://0.0.0.0:8000/TAEDProteinViewer.html` in a WebGL enabled browser. This requires `biopython` and `mafft` to be installed (among others).
