FIGURE 7 VERSION 12

Panel A illustrates the three landmarks on a generic pulse: peak, rising half-height and falling half-height. This is a schematic, not an experimental profile.
Panel B pools all 233 eligible peak and half-height comparisons in one summary graph, grouped into pair-rule, gap-gene and gap/pair-rule pairs.
The purple category is labeled exactly: Unresolved in time, separated in space.
Panels C-F retain the four original examples, now arranged horizontally in a single row. Each case is translated as an intact block without scaling or reflow. All 470 original SVG objects, including text, styling and 271 paths, are unchanged from version 11. The four panel letters remain C-F.

FORMATS
Fig7_ver12.svg: one editable layer with live text and vector paths.
Fig7_ver12.pdf: vector artwork, with zero embedded raster images.
Fig7_ver12.tiff: RGB TIFF, 600 dpi, LZW compression, 9000 x 5334 pixels.
Fig7_ver12.png: preview.
Landscape canvas: 1080 x 640 points, or 15 x 8.89 inches.

REBUILD
Install Python packages lxml, PyMuPDF and Pillow, plus Inkscape on the executable path.
Run: python build_fig7_ver12.py

INPUTS
source/Fig7_ver10.svg supplies the four original case panels.
source/S2_Data_Landmark_Comparisons.csv supplies all 859 candidates and eligibility/category fields for independent count verification.
source/summary_counts_and_percentages.csv is the current S3 Data. The builder plots only the 15 rows with panel=all_landmarks and figure_panel=B. The retained 30 stratified rows are not plotted in panel B.

VALIDATION
The builder checks all pooled counts against both S3 and eligible S2 rows. It verifies preservation of all 470 original SVG objects, including 271 plotted paths. The saved QA report also compares each repositioned case with its version-11 rendering and confirms that all SVG content remains exact. The final landscape figure has been visually inspected for clipping, overlap and readability.

Final submission ver04: the existing S2/S3 CSV filenames are retained as stable internal sources. Their tables are consolidated into the Landmark inventory and Summary source sheets of S1 Data. The figure artwork is unchanged.
