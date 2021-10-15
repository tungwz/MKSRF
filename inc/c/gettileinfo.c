#include "modis.h"

int GetTileInfo()
{
  intn            status;
  int32           gdfid, GDid;
  int             err = 0;
  char            filename[200];

  // GET THE LAI FILE NAME
  GetFileName(2000, 57, hh, vv, 1, filename);

  // OPEN HDF FILE
  gdfid = GDopen(filename, DFACC_READ);

  if (gdfid != -1)
  {
    GDid = GDattach(gdfid, modis[1].gdname);

    // GET THE GRID INFO AND PROJECTION INFO
    if (GDid != -1)
    {
      // GET THE GRID LOCATON INFORMATION
      status = GDgridinfo(GDid, &xdim, &ydim, uplft, lowrgt);

      // GET THE PROJECTION INFORMATION
      status = GDprojinfo(GDid, &projcode, &zonecode, &spherecode, projparm);
    } else {
      err = 2;
    }
  } else {
    err = 1;
  }

  // CLOSE HDF FILE
  GDdetach(GDid);
  GDclose(gdfid);

  return err;
}
