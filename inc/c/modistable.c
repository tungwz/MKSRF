#include "modis.h"

// ===================================================
// charactor replacement
// ===================================================
char* CharReplace(char *src,char oldChar,char newChar){
   char *head=src;
   while(*src!='\0'){
     if(*src==oldChar) *src=newChar;
     src++;
   }
   return head;
}

// ===================================================
// Initialize MODIS products data type
// ===================================================

void InitModisDataType(int id)
{
  ReadModisTable();
  idx = GetModisDataIndex(id);
  SetResPar(idx);
}

// ===================================================
// Read MODIS products table
// ===================================================

void ReadModisTable()
{
  int i;
  FILE *pf;
  char buf[1000];
  char fmodistable[] = "/home/yuanhua/modis/tools/data/modistable.dat";
  
  if ((pf=fopen(fmodistable, "r")) == NULL) {
    printf("Can't find file %s, stop!\n", fmodistable);
    exit(1);
  }

  // Read the table header
  fscanf(pf, "%[^\n]s", buf);

  i = 0;
  while (!feof(pf)) {
    fscanf(pf, "%d%s%d%s%s%s%d%s%d%d%d%d%s%d%d%d%lf%s%[^\n]s",
        &modis[i].pid, modis[i].pname,
        &modis[i].id, modis[i].name,
        modis[i].gdname, modis[i].fldname, &modis[i].dim, modis[i].header, 
        &modis[i].xdim, &modis[i].ydim, &modis[i].day, &modis[i].ndays,
        modis[i].dtype, &modis[i].novalue, &modis[i].fillvalue, &modis[i].valid, &modis[i].scale,
        modis[i].dir,
        buf);
    if (modis[i].pid == 6) {
      CharReplace(modis[i].fldname, '_', ' ');
    }
    i++;
  }

  fclose(pf);
}

// ===================================================
// Get the table index according to the product data id
// ===================================================

int GetModisDataIndex(int id) 
{
  int i;

  for (i = 0; i < 50; i++) {
    if (modis[i].id == id) {
      return i;
    }
  }
  
  return -1;
}

// ===================================================
// Set MODIS data resolution parameters
// ===================================================

void SetResPar(int idx)
{
  XDIM = modis[idx].xdim;
  YDIM = modis[idx].ydim;
  NDAYS = modis[idx].ndays;
}
