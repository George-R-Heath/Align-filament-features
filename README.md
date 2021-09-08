# filament-correlate
Matlab code used to align filament/fibril/polymer features such as helical repeats in image files.
Run each of the following codes in succession to detect (F1), analyze (F2) and correlate (F3) features.
- F1_detect_skel: Filament detection is performed using thresholding and skeletonization (does not perform well for high density, branched or overlapping filaments) 
- F2_analyze_skel: Filaments are then contoured smoothly, analyzed and digitally straightened 
- F3_filament_correlate: Using a reference frame/section features features on the filament are cross correlated to find similar features.
