#include "tools.h"
#include "netcdf.h"

// NOTE: the below is only applied for MODIS LC data
//
/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define error(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

// Define the scope
double lon_min = -180.0; double lon_max = 180.0;
double lat_min = -90.0;  double lat_max =  90.0;
double lon_min; double lon_max;
double lat_min; double lat_max;

// Spatial resolution (input)
double delta = 1.0/240.0;   // 15 sec
int rows, cols;

// Start and end year
int sy = 2000, ey = 2020;

// Data define
DTYPE **xydata;
unsigned short **xydata1;

// Initialization
void Init();

// Process data
void ProcessData();

// Write out data for one day
void WriteOutData(int, int, int, double, double, double, double);

// Free memory
void Free();

int main(int argc, char *argv[])
{
  int year, day;
  int ilat, ilon;
  int hh, vv;
  double lat, lon;
  POINT pp;

  // Initialization for map projection
  sinforint(6371007.181, 0, 0, 0);

  // Initial MODIS products parameters
  int id = atoi(argv[1]);
  InitModisDataType(id);

  Init();

  // Day loop
  for (day = 0; day < NDAYS; day++) // Loop for each 8-day step
  {
    // Year loop
    for (year = sy; year <= ey; year=year+5)
    {
      // Reproject the data using nearest sampling method
      Echos("Reproject the data using nearest sampling method ...")
//#pragma omp parallel for num_threads(2) private(ilat, ilon, lat, lon, pp)
      for (ilat = 0; ilat < rows; ilat++)
      {
        for (ilon = 0; ilon < cols; ilon++)
        {
          lon = lon_min + delta/2 + ilon*delta;
          lat = lat_max - delta/2 - ilat*delta;
          pp = LL2XY(lon, lat);
          xydata[ilat][ilon] = GetXY(&pp, idx, year, day);
        }
      }

      if (trueopen == 1) {
        ProcessData();
        WriteOutData(idx, year, day, lat_max, lon_min, lat_min, lon_max);
      }

      ClearOpened();
      trueopen = 0;
    } // End year loop
  } // End day loop

  // free memory
  Free();
  return 0;
}

// ===================================================
// Process data
// ===================================================
void ProcessData()
{
  int ncid, varid;
  int i, j, i1, j1, retval;
  int kg;

  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
   * to the file.*/
  if ((retval = nc_open("/home/yuanhua/tera02/Beck/Beck_KG_V1_present_0p0083.nc", NC_NOWRITE, &ncid)))
    error(retval);

  /* Get the varid of the data variable, based on its name. */
  if ((retval = nc_inq_varid(ncid, "zonecode", &varid)))
    error(retval);

  /* Read the data. */
  if ((retval = nc_get_var(ncid, varid, &xydata1[0][0])))
    error(retval);

  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    error(retval);

  printf("*** SUCCESS reading file %s!\n", "Beck climate zone code data");

  printf("*** Processing data with Koppen-Geiger zone code...\n");
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
      i1 = (int)(i/2);
      j1 = (int)(j/2);
      kg = xydata1[i1][j1];

      // in case of OCEAN, set land cover type to 0
      if (kg == 0) {xydata[i][j] = 0; continue;}

// yuan, 1/2/2020: set it to ocean 0
// barren (16) may be inconsistant with soil data
// NOTE!!!: need to re-run this progrom 与mkmod.r相匹配
// mkmod.r需要重新run
      // in case of FILLVALUE, set it to barren 16
      //if (xydata[i][j] == 255) {xydata[i][j] = 16;}
      if (xydata[i][j] == 255) {xydata[i][j] = 0;}
    }
  }
}

// ===================================================
// Output one day data in binary format
// ===================================================

void WriteOutData(int idx, int year, int day,
    double lat1, double lon1, double lat2, double lon2)
{
  FILE *po, *pv;
  char foutdata[200];
  int i;

  sprintf(foutdata, "/home/yuanhua/tera02/mksrf/lc_15s/MODIS_IGBP_%s.%04d%03d.global15s",
      modis[idx].name, year, day*modis[idx].day+1);

  printf("*** Write out file %s!\n", foutdata);
  po = fopen(foutdata, "wb");

  for (i = 0; i < rows; i++) {
    fwrite(xydata[i], sizeof(DTYPE), cols, po);
  }

  fclose(po);
}

void Init()
{
  int i, j;

  // High resolution lines/samples
  rows = (int)((lat_max-lat_min)/delta);
  cols = (int)((lon_max-lon_min)/delta);

  // Initial xydata
  xydata = (DTYPE **)malloc(sizeof(DTYPE *)*rows);
  xydata[0] = (DTYPE *)malloc(sizeof(DTYPE)*rows*cols);

  if (!xydata[0] || !xydata) {
    Err("Alloc memory error! stop!\n");
  }

  for (i = 1; i < rows; i++) {
    xydata[i] = xydata[i-1] + cols;
  }

  // Initial xydata1
  xydata1 = (unsigned short **)malloc(sizeof(unsigned short *)*rows/2);
  xydata1[0] = (unsigned short *)malloc(sizeof(unsigned short)*rows*cols/4);

  if (!xydata1[0] || !xydata1) {
    Err("Alloc memory error! stop!\n");
  }

  for (i = 1; i < rows/2; i++) {
    xydata1[i] = xydata1[i-1] + cols/2;
  }

  InitTile();
}

void Free()
{
  int i, j;

  free(xydata[0]);
  free(xydata);

  FreeTile();
}
