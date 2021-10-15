#include "modis.h"

void WriteBin(char filename[], DTYPE **data)
{
  FILE *pf;
  
  pf = fopen(filename, "wb");
  fwrite(data[0], sizeof(DTYPE), YDIM*XDIM, pf);

  fclose(pf);
}
