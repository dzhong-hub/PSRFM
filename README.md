# PSRFM
This PSRFM (Prediction Smooth Reflectance Fusion Model) software is a feature rich C++ program for blending Landsat or Sentinel-2 satellite images (higher spatial resolution but lower temporal frequency) with MODIS images ((lower spatial resolution but higher temporal frequency) to generate an exploitation-ready time series of synthetic images (higher spatial resolution and temporal frequency). The current version was developed and compiled by using the Microsoft Visual Studio 2017.

To run the compiled executable PSRFM_Main.exe, you just need to prepare an input parameter file (e.g. PSRFM_input.txt) to specify the input data files, output folders and your processing options etc. Then type the dos command "PSRFM_Main.exe PSRFM_Input.txt" in a command execution window of Windows 7 or 10.

<b>Important:</b> 
All the input image files must be:
1. reprojected to a same mapping projection system (e.g. UTM), 
2. saved in the <b>Int16</b> binary format and 
3. named with their acquisition date in the format <b>"*_ddd_dd-mmm-yyyy.dat"</b> at the end of their file name, for example: E:\\PSRFMTest\\Input\\L2A_T13UCT_clip<b>_143_23-May-2018.dat</b> (the ddd is the number of days of year).
4. croped to match the image size (number of rows and colums) which is an integer N times the resolution ratio between the fine and coarse images. For example, the fine image resolution is 30m, the corase image resolution is 500m, take the resolution ratio 500/30 as an integer 16, then the numbers of rows and colums of the input fine images must be an integer N times 16, i.e. N*16. The numbers of rows and colums of the input coarse images must be the same for co-registered images or larger for not co-registered images with the option for co-registration enabled. 

The mask files (with value 0 for normal and 1 for cloudy pixel) used to mark cloudy pixels must be saved in the <b>Int8</b> binary format. For the co-registration between the input fine (e.g. the Landsat or Sentinel-2 images) and coarse images (e.g. the MODIS images), the coarse images should cover a bigger area than the fine images. A coverage extension of 1 or 2 km around the boarder of the fine images is recommended. 

<b>More details refer to the example of the input parameter file: PSRFM_Input.txt</b>

If you have questions, please contact Detang Zhong (detang.zhong@canada.ca) and Fuqun Zhou (Fuqun.Zhou@canada.ca) at Canada Centre for Remote Sensing, Canada Centre for Mapping and Earth Observation, Natural Resources Canada.

If the software is helpful for your research and application developments, <b>please cite the following papers:</b>

1.	Zhong, D.; Zhou, F. Improvement of Clustering Methods for Modelling Abrupt Land Surface Changes in Satellite Image Fusions, Remote Sens. 2019, 11, 1759; doi:10.3390/rs11151759.

2.	Zhong, D.; Zhou, F. A Prediction Smooth Method for Blending Landsat and Moderate Resolution Imagine Spectroradiometer Images, Remote Sens. 2018, 10 (9), 1371; doi:10.3390/rs10091371


