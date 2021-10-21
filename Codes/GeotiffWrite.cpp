#include <iostream>
#include <string>
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <gdalwarper.h>
#include <stdlib.h>

using namespace std;
typedef std::string String; 
 
class GeotiffWrite { 
 
  private: 
    const char* filename;        // name of Geotiff
    GDALDataset *geotiffDataset; // Geotiff GDAL dataset object. 
	GDALDataset* srcDataset; // Geotiff GDAL dataset object. 
    GDALDriver *driverTiff;
    int dimensions[3];           // X,Y, and Z dimensions. 
    int NROWS,NCOLS,NLEVELS;     // dimensions of data in Geotiff. 
 
  public: 
	  //GeotiffWrite(const char* dst, const char* src) {
		 // const char* pszFormat = "GTiff";
		 // GDALDriver* poDriver;
		 // char** papszMetadata;
		 // poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		 // if (poDriver == NULL)
			//  exit(1);
		 // papszMetadata = poDriver->GetMetadata();
		 // if (CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, FALSE))
			//  printf("Driver %s supports Create() method.\n", pszFormat);
		 // if (CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATECOPY, FALSE))
			//  printf("Driver %s supports CreateCopy() method.\n", pszFormat);
		 // GDALDataset* poSrcDS =
			//  (GDALDataset*)GDALOpen(src, GA_ReadOnly);
		 // GDALDataset* poDstDS;
		 // poDstDS = poDriver->CreateCopy(dst, poSrcDS, FALSE,
			//  NULL, NULL, NULL);
		 // /* Once we're done, close properly the dataset */
		 // //if (poDstDS != NULL)
			// // GDALClose((GDALDatasetH)poDstDS);
		 // GDALClose((GDALDatasetH)poSrcDS);

	  //}
	  GeotiffWrite(const char* dst, const char* src) {
		  srcDataset = (GDALDataset*)GDALOpen(src, GA_ReadOnly);
		  // set the dimensions of the Geotiff 
		  NROWS = GDALGetRasterYSize(srcDataset);
		  NCOLS = GDALGetRasterXSize(srcDataset);
		  NLEVELS = GDALGetRasterCount(srcDataset);
		  cout << NROWS << " " << NCOLS << " " << NLEVELS << endl;
		  driverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");

		  // set pointer to Geotiff dataset as class member.  
		  geotiffDataset = driverTiff->Create(dst, NCOLS, NROWS, NLEVELS, GDT_Float32, NULL);
		  double adfGeoTransform[6]; 
		  srcDataset->GetGeoTransform(adfGeoTransform);
		  geotiffDataset->SetGeoTransform(adfGeoTransform);
		  geotiffDataset->SetProjection(srcDataset->GetProjectionRef());

	  }

    GeotiffWrite( const char* tiffname, int rows, int cols, int levels,double* geotransform, const OGRSpatialReference* spatialreference ) {
      filename = tiffname ; 
      NROWS = rows;
      NCOLS = cols;
      NLEVELS = levels;

      driverTiff = GetGDALDriverManager()->GetDriverByName("GTiff");
      // set pointer to Geotiff dataset as class member.  
      geotiffDataset = driverTiff->Create(filename, NCOLS, NROWS, NLEVELS, GDT_Float32, NULL);
	  geotiffDataset->SetGeoTransform(geotransform);
	  geotiffDataset->SetSpatialRef(spatialreference);
    }
 
    ~GeotiffWrite() {
      // close the Geotiff dataset, free memory for array.  
      GDALClose(geotiffDataset);
      // GDALDestroyDriverManager();
    }
 
    const char *GetFileName() { 
      return filename; 
    }
 
    int *GetDimensions() {
      /* 
       * int *GetDimensions(): 
       * 
       *  This function returns a pointer to an array of 3 integers 
       *  holding the dimensions of the Geotiff. The array holds the 
       *  dimensions in the following order:
       *   (1) number of columns (x size)
       *   (2) number of rows (y size)
       *   (3) number of bands (number of bands, z dimension)
       */
      dimensions[0] = NROWS; 
      dimensions[1] = NCOLS;
      dimensions[2] = NLEVELS; 
      return dimensions;  
    }

    void write(float** data) {
        float *rowBuff = (float*) CPLMalloc(sizeof(float) * NCOLS);
        for (int row = 0; row < NROWS; row++)
        {
            for (int col = 0; col < NCOLS; col++)
            {
                rowBuff[col] = data[row][col];
            }
            geotiffDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, row, NCOLS, 1, rowBuff, NCOLS, 1, GDT_Float32, 0, 0);
        }
        CPLFree( rowBuff );
        cout << "Wrote the data to the file: " << filename << endl;
    }
 
};