# stSCA_scripts
Analysis scripts used for stSCA manuscript

signal_proc_functions is a folder that contains custom scripts for various signal processing procedures used throughout the scripts 
utah_maps is a folder containing spatial map of the utah array channels. They are specific to the site of data collection (i.e. Columbia and Harvard) 
all_files_metadata.xlsx contains information such as seizure start and end times for the recordings 

Recordings are identified by (1) the patient ID, and (2) the recording ID
Example: c5_s1 is the recording for patient c5, seizure 1. Interictal recordings do not have a number (e.g. c5_interictal)

All analyses assume that the microelectrode array recordings have been standardized into the following files: 
(c5_s1 is used as an example to demonstrate the naming convention)

c5_s1_mua.mat is the data bandpass filtered at 300-3000Hz
c5_s1_llfp.mat is the data bandpass filtered at 2-50Hz
c5_s1_whitened_llfp.mat is the xx_xx_llfp.mat file that has been decorrelated using the script whiten_signals.m

These standardized data files may be shared upon request.
