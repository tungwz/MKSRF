#include "tools.h"

// Define the scope
double lon_min; double lon_max;
double lat_min; double lat_max;

// Spatial resolution (input)
double delta = 1.0/480.0;   // 7.5 sec
int rows, cols;

// Start and end year
int sy = 2000, ey = 2020;

// Data define
int **xydata;
DTYPE **outdata;

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
  for (day = 8; day < 9; day++) // Loop for each 8-day step
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

      WriteOutData(idx, year, day, lat_max, lon_min, lat_min, lon_max);
      ClearOpened();
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
  int i, j, cnt, sum2;
  double sum;

  sprintf(foutdata, "/tera02/yuanhua/mksrf/vcf_5x5/RG_%d_%d_%d_%d.%s%04d%03d",
      (int)lat1, (int)lon1, (int)lat2, (int)lon2,
      modis[idx].name, year, day*modis[idx].day+1);

  po = fopen(foutdata, "wb");

//#pragma omp parallel for private(i, j, sum, cnt)
  for (i = 0; i < rows; i=i+2)
  {
    for (j = 0; j < cols; j=j+2)
    {
      sum  = 0;
      sum2 = 0;
      cnt  = 0;

      // water: 200, fillvalue: 253
      if (xydata[i][j]     <= 100) {
        sum = sum + xydata[i][j]; cnt = cnt + 1;
      } else {sum2 = sum2 + xydata[i][j];}

      if (xydata[i][j+1]   <= 100) {
        sum = sum + xydata[i][j+1]; cnt = cnt + 1;
      } else {sum2 = sum2 + xydata[i][j+1];}

      if (xydata[i+1][j]   <= 100) {
        sum = sum + xydata[i+1][j]; cnt = cnt + 1;
      } else {sum2 = sum2 + xydata[i+1][j];}

      if (xydata[i+1][j+1] <= 100) {
        sum = sum + xydata[i+1][j+1]; cnt = cnt + 1;
      } else {sum2 = sum2 + xydata[i+1][j+1];}

      if (cnt != 0) {               // if data exist
        outdata[i/2][j/2] = (int)(sum/cnt + 0.5);
      } else if (sum2%10 > 0) {     // if fillvalue exsits
        outdata[i/2][j/2] = 253;
      } else {                      // all 200, water case
        outdata[i/2][j/2] = 0;
      }
    }
  }

  fwrite(outdata[0], sizeof(DTYPE), rows*cols/4, po);

  fclose(po);
}

void Init()
{
  int i, j;

  // High resolution lines/samples
  rows = (int)((lat_max-lat_min)/delta);
  cols = (int)((lon_max-lon_min)/delta);

  // Initial xydata
  xydata = (int **)malloc(sizeof(int *)*rows);
  xydata[0] = (int *)malloc(sizeof(int)*rows*cols);

  if (!xydata[0] || !xydata) {
    Err("Alloc memory error! stop!\n");
  }

  for (i = 1; i < rows; i++) {
    xydata[i] = xydata[i-1] + cols;
  }

  // Initial outdata
  outdata = (DTYPE **)malloc(sizeof(DTYPE *)*rows/2);
  outdata[0] = (DTYPE *)malloc(sizeof(DTYPE)*rows*cols/4);

  if (!outdata[0] || !outdata) {
    Err("Alloc memory error! stop!\n");
  }

  for (i = 1; i < rows/2; i++) {
    outdata[i] = outdata[i-1] + cols/2;
  }

  InitTile();
}

void Free()
{
  int i, j;

  free(xydata[0]);
  free(xydata);
  free(outdata[0]);
  free(outdata);

  FreeTile();
}
