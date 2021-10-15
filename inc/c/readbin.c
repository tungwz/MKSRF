#include "modis.h"

int ReadBin(char filename[], DTYPE **data)
{
  FILE *pf;
  int err = 0;

  pf = fopen(filename, "rb");

  if (!pf) {
    fread(data[0], sizeof(DTYPE), YDIM*XDIM, pf);
  } else {
    err = 1; 
  }

  fclose(pf);

  return err;
}
