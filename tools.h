#ifndef _TOOLS_H_
#define _TOOLS_H_

#include "modis.h"

// Typedef point structure
typedef struct point {
  int hh;
  int vv;
  int x;
  int y;
  int index;
  struct point *ppoint;
} POINT;

#define LEN sizeof (struct point)

// Typedef tile structure
typedef struct tile {
  int opened;
  int newed;
  DTYPE **pdata;
} TILE;

// PFT name of MOD12Q1 land cover type 5
static char pftname[13][30] = {
  "water",
  "evergreen needleleaf tree",
  "evergreen broadleaf tree",
  "deciduous needleleaf tree",
  "deciduous broadleaf tree",
  "shrub",
  "grass",
  "cereal crop",
  "broadleaf crop",
  "urban and built-up",
  "snow and ice",
  "barren or sparely vegetated",
  "unclassified"
};

// Global tiles
TILE tiles[36][18];

int trueopen = 0;

// ===================================================
// Set the value of no file data
// ===================================================

void SetNoValue(DTYPE **data)
{
  int i, j;

  if (sizeof(DTYPE) == 1) {
    memset(data[0], modis[idx].novalue, XDIM*YDIM);
  } else {
    for (i = 0; i < XDIM; i++) {
      for (j = 0; j < YDIM; j++) {
        data[i][j] = modis[idx].novalue;
      }
    }
  }
}

// ===================================================
// open hdf tile file
// ===================================================

void OpenTile(POINT *point, int idx, int year, int day)
{
  char filename[200];
  int hh, vv, i;

  hh = point->hh; vv = point->vv;

  if (!tiles[hh][vv].newed) {

    tiles[hh][vv].pdata =
      (DTYPE **)malloc(sizeof(DTYPE *)*XDIM);
    tiles[hh][vv].pdata[0] =
      (DTYPE *)malloc(sizeof(DTYPE)*XDIM*YDIM);

    for (i = 1; i < YDIM; i++) {
      tiles[hh][vv].pdata[i] = tiles[hh][vv].pdata[i-1] + YDIM;
    }

    if (!tiles[hh][vv].pdata[0] || !tiles[hh][vv].pdata) {
      Err("Alloc memory error! stop!\n");
    }

    tiles[hh][vv].newed = 1;
  }

  tiles[hh][vv].opened = 1;

  GetFileName(year, day, point->hh, point->vv, idx, filename);

  if (ReadHdf(filename, idx, tiles[hh][vv].pdata)) {
    NoFile(filename);
    SetNoValue(tiles[hh][vv].pdata);
  } else { trueopen = 1;}
}

void OpenTileD(int hh, int vv, int idx, int year, int day)
{
  int i;
  char filename[200];

  if (!tiles[hh][vv].newed) {

    tiles[hh][vv].pdata =
      (DTYPE **)malloc(sizeof(DTYPE *)*XDIM);
    tiles[hh][vv].pdata[0] =
      (DTYPE *)malloc(sizeof(DTYPE)*XDIM*YDIM);

    for (i = 1; i < YDIM; i++) {
      tiles[hh][vv].pdata[i] = tiles[hh][vv].pdata[i-1] + YDIM;
    }

    if (!tiles[hh][vv].pdata[0] || !tiles[hh][vv].pdata) {
      Err("Alloc memory error! stop!\n");
    }

    tiles[hh][vv].newed = 1;
  }

  tiles[hh][vv].opened = 1;

  GetFileName(year, day, hh, vv, idx, filename);

  if (ReadHdf(filename, idx, tiles[hh][vv].pdata)) {
    NoFile(filename);
    SetNoValue(tiles[hh][vv].pdata);
  } else { trueopen = 1;}
}

// ===================================================
// Get the one point value in the tile
// ===================================================

DTYPE SetXY(POINT *point, int idx, int year, int day)
{
  int hh, vv, x, y;
  DTYPE **p;
  int value;

  hh = point->hh; vv = point->vv;
  x = point->x-1; y = point->y-1;

  p = tiles[hh][vv].pdata;
  value = p[y][x];

  return value;
}

DTYPE GetXY(POINT *point, int idx, int year, int day)
{
  if (!(tiles[point->hh][point->vv]).opened) {
    OpenTile(point, idx, year, day);
  }
  return SetXY(point, idx, year, day);
}

// ===================================================
// Lat/lon to tile hh/vv
// ===================================================

POINT LL2XY(double lon, double lat)
{
  int hh, vv, dx, dy;
  double x1, x2, y1, y2, x, y, temp;
  POINT pp;

  hh = (int)(18+lon*cos(lat*D2R)/10);
  vv = (int)(9-lat/10);

  sinfor((hh-18)*10*D2R, 0*D2R, &x1, &temp);
  sinfor((hh-17)*10*D2R, 0*D2R, &x2, &temp);
  sinfor(0*D2R, (9-vv)*10*D2R, &temp, &y1);
  sinfor(0*D2R, (8-vv)*10*D2R, &temp, &y2);
  sinfor(lon*D2R, lat*D2R, &x, &y);

  dx = (int)((x-x1)*modis[idx].xdim/(x2-x1)+1);
  dy = (int)((y-y1)*modis[idx].ydim/(y2-y1)+1);

  pp.x = dx; pp.y = dy;
  pp.hh = hh; pp.vv = vv;

  return pp;
}

// ===================================================
// Initial tile open/new flag
// ===================================================

void InitTile()
{
  int i, j;

  // Loop for hh/vv
  for (i = 0; i < 36; i++) {
    for (j = 0; j < 18; j++) {
      tiles[i][j].opened = 0;
      tiles[i][j].newed = 0;
      tiles[i][j].pdata = NULL;
    }
  }
}

// ===================================================
// Clear the opend flags
// ===================================================

void ClearOpened()
{
  int i, j;

  // Loop for hh/vv
  for (i = 0; i < 36; i++) {
    for (j = 0; j < 18; j++) {
      tiles[i][j].opened = 0;
    }
  }
}

// ===================================================
// Free memory
// ===================================================

void FreeTile()
{
  int i, j;

  // Loop for hh/vv
  for (i = 0; i < 36; i ++) {
    for (j = 0; j < 18; j++) {
      if (tiles[i][j].newed) {
        free(tiles[i][j].pdata[0]);
        free(tiles[i][j].pdata);
      }
    }
  }

  Echos("Free tiles memory");
}

#endif
