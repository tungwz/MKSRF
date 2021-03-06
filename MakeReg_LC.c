#include "tools.h"

// Define the scope
double lon_min; double lon_max;
double lat_min; double lat_max;

// Spatial resolution (input)
double delta = 1.0/240.0;   // 15 sec
int rows, cols;

// Start and end year
int sy = 2000, ey = 2020;

// Data define
DTYPE **xydata;

// Initialization
void Init();

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

  lat_max = atof(argv[2]);
  lon_min = atof(argv[3]);
  lat_min = atof(argv[4]);
  lon_max = atof(argv[5]);
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
// Output one day data in binary format
// ===================================================

void WriteOutData(int idx, int year, int day,
    double lat1, double lon1, double lat2, double lon2)
{
  FILE *po, *pv;
  char foutdata[200];
  int i, j;

  sprintf(foutdata, "/tera02/yuanhua/mksrf/lc_5x5/RG_%d_%d_%d_%d.%s%04d%03d",
      (int)lat1, (int)lon1, (int)lat2, (int)lon2,
      modis[idx].name, year, day*modis[idx].day+1);

  po = fopen(foutdata, "wb");

  fwrite(xydata[0], sizeof(DTYPE), rows*cols, po);

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

  InitTile();
}

void Free()
{
  int i, j;

  free(xydata[0]);
  free(xydata);

  FreeTile();
}
