#ifndef _MODIS_H_
#define _MODIS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

// Gctpc header files
#include "cproj.h"
#include "proj.h"

// HDF and HDF-EOS header files
#include "hdf.h"
#include "HdfEosDef.h"

#define UINT8
//#define UINT16

#ifdef UINT8
typedef unsigned char DTYPE;
#endif

#ifdef UINT16
typedef unsigned short int DTYPE;
#endif

// ===================================================
// Define some macros
// ===================================================

#define Err(s)      {printf("%s\n", s); exit(1);}
#define Echos(s)    {printf("%s\n", s);}
#define Echoi(s)    {printf("%d\n", s);}
#define NoFile(s)   {printf("\nCan't Read %s, Set The Filled Value!\n", s);}
#define ErrFile(s)  {printf("\nCan't Read %s! Stop!\n", s); exit(1);}

// ===================================================
// MODIS products namelist sturcture
// ===================================================

typedef struct modis_table {
  int pid, id;
  int xdim, ydim, dim, day, ndays;
  int novalue, fillvalue, valid;
  double scale;
  char pname[50];
  char name[50];
  char gdname[50];
  char fldname[50];
  char dtype[50];
  char header[50];
  char dir[200];
} MODIS_TABLE;

// ===================================================
// Global variables
// ===================================================

// Modis products table
MODIS_TABLE modis[50];
int idx;

// Modis resolution parameters
int XDIM, YDIM;
int NDAYS;

// Tile information
char *hhvv;
char h[3];
char v[3];
int hh, vv;

// ===================================================
// Projection information and other
// ===================================================

int32   xdim, ydim, zonecode, projcode, spherecode;
int32   compcode;
intn    compparm[5];
float64 projparm[16], uplft[2], lowrgt[2];

// ===================================================
// Functions
// ===================================================

// Read MODIS products table
void ReadModisTable();

// Get MODIS data index
int GetModisDataIndex(int);

// Set MODIS data resolution parameters
void SetResPar(int);

// Initialize MODIS products data type
void InitModisDataType(int);

// To get the projection information
int GetTileInfo();

// Get file name
void GetFileName(int, int, int, int, int, char []);

// Read data
int ReadLaiQcHdf(char [], DTYPE **, DTYPE **);
int ReadHdf(char [], int, DTYPE **);
int ReadBin(char [], DTYPE **);

// Write out data
void WriteLaiQcHdf(char [], DTYPE **, DTYPE **);
void WriteHdf(char [], DTYPE **);
void WriteBin(char [], DTYPE **);

// Print the processing bar
void PrintProcBar(int, int);

#endif
