#include "modis.h"

void GetFileName(int year, int day, int hh, int vv, 
    int idx, char filename[])
{  
  if (modis[idx].id == 24 || modis[idx].id == 25) {
    sprintf(filename, "%s%s.h%02dv%02d/%s%d%03d.h%02dv%02d.tsf.sat.hdf",
        modis[idx].dir, "MCD15A2H", hh, vv,
        modis[idx].header, year, day*modis[idx].day+1, hh, vv);
  } else {
    sprintf(filename, "%s%s%d%03d/%s%d%03d.h%02dv%02d.hdf", 
        modis[idx].dir, modis[idx].header, year, day*modis[idx].day+1,
        modis[idx].header, year, day*modis[idx].day+1, hh, vv);
  }
}
