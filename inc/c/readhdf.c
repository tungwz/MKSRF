#include "modis.h"

int ReadHdf(char filename[], int idx, DTYPE **data)
{
  intn status;
  int32 gdfid, GDid;

  int err = 0;
  int i, j;
  unsigned short int tmp[XDIM][YDIM][3];

  gdfid = GDopen(filename, DFACC_READ);
  printf("--- %s ---\n", filename);

  if (gdfid != -1)
  {
    GDid = GDattach(gdfid, modis[idx].gdname);

    if (GDid != -1)
    {
      if (modis[idx].pid == 4) {
        status = GDreadfield(GDid, modis[idx].fldname, 
            NULL, NULL, NULL, tmp);
      } else {
        status = GDreadfield(GDid, modis[idx].fldname, 
            NULL, NULL, NULL, data[0]);
      }
      if (status != 0) {
        err = 3;
      }
    } else {
      err = 2;
    }
  } else {
    err = 1;
  }

  GDdetach(GDid);
  GDclose(gdfid);

  if (modis[idx].pid == 4 && (!err)) {
    for (i = 0; i < XDIM; i++) {
      for (j = 0; j < YDIM; j++) {
        data[i][j] = tmp[i][j][modis[idx].dim];
      }
    } 
  }

  return err;
}
