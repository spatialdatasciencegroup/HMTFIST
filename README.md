
## Compilation Code
```
g++ -std=c++11 -lgdal -I(PATH_TO_GDAL_INCLUDE)  _L(PATH_TO_GDAL_LIB) GeotiffRead.cpp GeotiffWrite.cpp HMTFIST.cpp -o hmtfst --verbose
```
## Usage:

1. List of input data:
  - a. Source direction layer. (tif file)
  - b. Bank layer.(tif file)
  - c. Cost layer.(tif file)
  - d. Pits layer.(tif file)
  - e. Tree canopy layer.(tif file)
  - f. Roughness layer.(tif file)
  - g. Elevation layer.(tif file)
  - h. parameter file. (text file)
  - i. Stream Data (csv file)

2. Create Configuration File (See: Setup.xlsx->"ConfigurationSetup" worksheet for more details) Example: config.txt

3. Run

	HMTFIST.exe config.txt

## Reference

Jiang, Zhe, and Arpan Man Sainju. "Hidden Markov Contour Tree: A Spatial Structured Model for Hydrological Applications." Proceedings of the 25th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining. ACM, 2019.
