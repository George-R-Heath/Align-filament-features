# Correlate-Filaments
Matlab code used to straighten and align filament/fibril/polymer features such as helical repeats in image files.
Run each of the following codes in succession to detect (F1), analyze (F2) and correlate (F3) features.
- F1_detect_skel: Filament detection is performed using thresholding and skeletonization (does not perform well for high density, branched or overlapping filaments) 
- F2_analyze_skel: Filaments are then contoured smoothly, analyzed and digitally straightened 
- F3_filament_correlate: Using a reference frame/section features features on the filament are cross correlated to find similar features.

This code was developed for analysis of ribonucleoprotein, if using please refer to and cite:

The structure of a native orthobunyavirus ribonucleoprotein reveals a key role for viral RNA in maintaining its helical architecture
Francis R. Hopkins, Beatriz Álvarez-Rodríguez, George R. Heath, Kyriakoulla Panayi, Samantha Hover, Thomas A. Edwards, John N. Barr, Juan Fontana
bioRxiv 2021.10.27.466080; doi: https://doi.org/10.1101/2021.10.27.466080
