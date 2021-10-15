#include "modis.h"

int ReadLaiQcHdf(char filename[], DTYPE **lai, DTYPE **qc)
{
  intn status;
  int32 gdfid, GDid;

  int err = 0;

  // open hdf file
  gdfid = GDopen(filename, DFACC_READ);

  if (gdfid != -1)
  {
    GDid = GDattach(gdfid, "MOD_Grid_MOD15A2");

    if (GDid != -1)
    {
      status = GDreadfield(GDid, "Lai_1km", 
          NULL, NULL, NULL, lai[0]);
      if (status != 0) {
        err = 3;
      }

      status = GDreadfield(GDid, "FparLai_QC", 
          NULL, NULL, NULL, qc[0]);

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

  return err;
}
