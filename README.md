
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

4. Output files:
  1. Tif file: Image file with flood map. 0 implies dry and 1 implies flood. -1 implies no data.
  2. ProfileTable.xlsx   (See Output_File_Details.xlsx-> ProfileTable)
  3. ProfileTable_preInterpolation.xlsx (See Output_File_Details.xlsx-> ProfileTable_preInterpolation)

## Reference
1. Jiang, Zhe, Miao Xie, and Arpan Man Sainju. "Geographical hidden markov tree." IEEE Transactions on Knowledge and Data Engineering (2019).

2. Jiang, Zhe, and Arpan Man Sainju. "Hidden Markov Contour Tree: A Spatial Structured Model for Hydrological Applications." Proceedings of the 25th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining. ACM, 2019.
