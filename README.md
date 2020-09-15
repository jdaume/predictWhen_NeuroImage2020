# predictWhen_NeuroImage_2020

Analysis code from Daume et al. (2020) Non-rhythmic temporal prediction involves phase resets of low-frequency delta oscillations. NeuroImage.

--------------------------------------------

<b> Simulation ITPC vs evoked potential (CNV): </b>

simulate_ITPCvsCNV.m 

<br>

<b> ITPC and power analysis:</b> 

wlconv_sensor_wholeTrial.m        (Whole-trial ITPC, power, and cross-spectra analysis on sensor level)

wlconv_sensor_stimOffset.m        (ITPC and power analysis in enlarged window time-locked to disappearance on sensor level)

wlconv_power_source_wholeTrial.m  (Whole-trial power analysis on source level)

wlconv_source_stimOffset.m        (ITPC and power analysis in enlarged window time-locked to disappearance on source level)

<br>

<b> Mixed-model regression analysis: </b>

meg_roi_mixed_model_regression.R


<br>

<b> Plotting and cluster-based permutation statistics: </b>

clusterTest_sensor_itpc_pow_stimOffset.m

clusterTest_sensor_itpc_wholeTrial.m

clusterTest_sensor_pow_wholeTrial.m

clusterTest_source_itpc_stimOffset.m (Includes correlation analysis between ITPC and behavior)

clusterTest_source_pow_wholeTrial.m





