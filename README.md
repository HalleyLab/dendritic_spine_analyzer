# Dendritic spine analyzer

Measure 2D fluorescence-positive dendritic spine candidate area per calibrated dendrite length (um2/um) from CZI images. The notebook includes source-image selection, horizontal views, manual shaft edges and add/erase brushes.

Use Python 3.12 or later. From the repository directory:

```bash
python -m pip install -r requirements.txt
python -m jupyterlab
```

Open `dendrite_spine_analysis.ipynb`, set `CZI_PATH` in the input-parameters cell to your CZI file, and run cells in order up to the selection interface. Start with a region near the first dendritic bifurcation above the soma. The shared notebook has no image-specific trace or saved outputs. Input parameters and short code comments are in English. `PROJECT_DIR` defaults to the directory where the notebook kernel starts; set it explicitly if needed. Notebook caches, optional dependencies and results stay inside that directory.

The optional desktop interface is `python app.py` (or `start.bat` on Windows); use the notebook for manual brush editing.

1. **Select region:** drag a new rectangle, or repeatedly drag its interior/handles to move/resize it. Rectangle mode stays active until you click **2 Trace shaft** or change Edit mode explicitly.
2. **Trace shaft:** place centerline points from the soma-near end outward; cover offset + length. Point placement does not run segmentation.
3. **Analyze / refresh:** run the segmentation once. Threshold edits then refresh the current result on release/Enter.
4. Run the final **Save result** code cell. Editing controls does not save files.

After analysis, choose **Paint erase spines** in Edit and drag over a misclassified branch on Original CZI. **Paint add spines** restores missed positive regions. Set **Brush diameter (um)**; release the mouse to update area. **Undo brush** (or right-click in paint mode) removes the last stroke; **Clear brush edits** restores the automatic result. Pan/zoom with the toolbar, then turn navigation off before painting.
Add only includes threshold-positive pixels inside the measured corridor, outside the shaft and exclusion polygons. It can override automatic connectivity, area and reach filters. Lower the spine threshold for dim signal; use the shaft-edge tools to change the shaft. Strokes survive threshold changes and Analyze refreshes; selecting a new rectangle or clearing the trace resets them.

Defaults match the supplied parameter screenshot: length 15 um, start offset 2 um, spine threshold 14024.49000011921, shaft reference 32744.2421875, maximum spine reach 10 um, minimum candidate area 0.001 um2 and maximum candidate area 0 (unlimited).
Keep `Shaft reference` fixed while lowering `Spine fluorescence threshold` to recover dim heads/necks. Setting `SHAFT_THRESHOLD=None` initializes an independent Otsu reference; 0 explicitly links it to the spine threshold.
Inspect View layers: All threshold-positive, Before area/reach filters, Rejected / detached signal, Accepted candidates. Missing in positive view suggests threshold/preprocessing; orange rejection suggests size/reach filtering; a spine inside green suggests an incorrect shaft edge.
The area table shows accepted and excluded areas separately. Lowering the threshold does not guarantee a larger final area: merging with a long branch can reject the whole component. Relax filters only after checking morphology and connectivity.

Advanced settings include shaft width/margin, Gaussian smoothing and background correction. Open the horizontal view to edit upper/lower shaft edges manually. Maximum candidate area 0 disables the upper area limit; components are not spine counts.
Two horizontal views show the same red fluorescence and detections: filled overlays, and outlines only. Green marks the shaft and magenta marks accepted spine candidates. Both views support manual shaft-edge editing and are saved by the final save cell as `horizontal_filled_VIEW_ONLY.png` and `horizontal_outlines_VIEW_ONLY.png`. Straightened pixels are for viewing; area is measured on the source image.

Save result writes full-size source PNG/RGB TIFF annotations, original raw TIFF, labels, masks, CSV and session JSON under `PROJECT_DIR/notebook_results/`. Yellow=region, green=shaft, magenta=accepted candidates. CZI raw data remains unchanged.
The final mask/area and both horizontal views include brush edits. `source_auto_spine_mask.tif`, `source_manual_added_mask.tif` and `source_manual_removed_mask.tif` retain the automatic result and net corrections. `session.json` records calibrated source-coordinate strokes for replay. `components.csv` describes automatic filtering; `final_components.csv` describes the final corrected mask.
Detailed Chinese manual ImageJ reproduction: `IMAGEJ_MANUAL_CN.md`. Short English reference: `IMAGEJ_GUIDE.md`.
Review saturation and biological segmentation manually; no expert ground-truth annotations were provided.

Core checks: `python -B -m unittest -v test_core`. See `VALIDATION.md` for the validation scope.
