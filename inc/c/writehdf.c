#include "modis.h"

void WriteHdf(char filename[], DTYPE **lai)
{
  intn            status;
  int32           gdfid, GDid;

  gdfid = GDopen(filename, DFACC_CREATE);

  // SET GRID INFORMATION
  GDid = GDcreate(gdfid, "MOD_Grid_MOD15A2", xdim, ydim, uplft, lowrgt);

  // SET PROJECTION INFORMATION
  status = GDdefproj(GDid, GCTP_SNSOID, zonecode, spherecode, projparm);

  // DEFINE COMPRESSION INFO
  compparm[0] = 5;
  status = GDdefcomp(GDid, HDFE_COMP_DEFLATE, compparm);

  // DEFINE LAI
  status = GDdeffield(GDid, "Lai_1km", "YDim,XDim", 
      DFNT_UINT8, HDFE_NOMERGE);

  // WRITE LAI
  status = GDwritefield(GDid, "Lai_1km", NULL, NULL, NULL, lai[0]);

  // CLOSE HDF FILE
  GDdetach(GDid);
  GDclose(gdfid);
}
