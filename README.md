# TAED Protein Viewer

Visualizing proteins from TAED, The Adaptive Evolution Database. Work in progress.

To run (which involves scraping some data that gets ran through a pipeline), enter:

```
git clone https://github.com/stephenshank/taed-pv
cd taed-pv
source run.sh
python -m http.server
```

and visit `http://0.0.0.0:8000/TAEDProteinViewer.html` in a WebGL enabled browser.
