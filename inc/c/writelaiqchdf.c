#include "modis.h"

void WriteLaiQcHdf(char filename[], DTYPE **lai, 
    DTYPE **qc)
{
  intn            status;
  int32           gdfid, GDid;

  gdfid = GDopen(filename, DFACC_CREATE);

  // set grid information
  GDid = GDcreate(gdfid, "MOD_Grid_MOD15A2", xdim, ydim, uplft, lowrgt);

  // set projection information
  status = GDdefproj(GDid, GCTP_SNSOID, zonecode, spherecode, projparm);

  // define compression info
  compparm[0] = 5;
  status = GDdefcomp(GDid, HDFE_COMP_DEFLATE, compparm);

  // define lai and qc field
  status = GDdeffield(GDid, "Lai_1km", "YDim,XDim", 
      DFNT_UINT8, HDFE_NOMERGE);
  status = GDdeffield(GDid, "FparLai_QC", "YDim,XDim", 
      DFNT_UINT8, HDFE_NOMERGE);

  // write lai and qc data
  status = GDwritefield(GDid, "Lai_1km", NULL, NULL, NULL, lai[0]);
  status = GDwritefield(GDid, "FparLai_QC", NULL, NULL, NULL, qc[0]);

  // close hdf file
  GDdetach(GDid);
  GDclose(gdfid);
}
