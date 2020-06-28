# SO_Mooring_Flux_Array

(1) Download monthly JRA reanalysis data from the link below: 
http://sose.ucsd.edu/DATA/OTHER/JRA_Qnet_197901_201612.mat

(2) Then run ‘SO_JRA_Qnet_DataFiltering.m’ to generate some pre-processed data. 
The subroutine ‘smooth2D_per.m’ is used to generate spatially filtered data, and it will take 1-2 hrs in this program.

(3) Run ‘SO_JRA_Qnet_EOF_Analysis.m’ for EOF analysis. 

(4) Run ‘SO_PlaceMooring.m’ to site moorings based on the criteria of constraining the largest local standard deviation or constraining the largest fraction of total variance. 

Matlab toolboxes m_map and seawater are required. A subroutine 'gft.m' is required and it can be downloaded from: https://github.com/mmazloff/gft

