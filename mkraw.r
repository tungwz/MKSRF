
# ----------------------------------------
# make the vegetation rawdata of land surface
# ----------------------------------------

library("ncdf4")

# set year
year = "2005" # default value, year 2005
args = commandArgs(T)
if (!is.na(args[1])) { year = args[1] }

# define directory
# ----------------------------------------
ROOT_DIR = "/home/yuanhua/tera02/mksrf/"

RAW_DIR  = paste(ROOT_DIR, "raw_5x5/", sep="")

LC_DIR   = paste(ROOT_DIR, "lc_5x5/", sep="")
LC_SUF   = paste(".LC", year, "001",  sep="")

PCTT_DIR = paste(ROOT_DIR, "vcf_5x5/", sep="")
PCTH_DIR = paste(ROOT_DIR, "vcf_5x5/", sep="")
PCTB_DIR = paste(ROOT_DIR, "vcf_5x5/", sep="")
PCTT_SUF = paste(".PCTT", year, "065", sep="")
PCTH_SUF = paste(".PCTH", year, "065", sep="")
PCTB_SUF = paste(".PCTB", year, "065", sep="")

HTOP_DIR = paste(ROOT_DIR, "htop_5x5/", sep="")
HTOP_SUF = ".TREETOP"

BT_DIR   = paste(ROOT_DIR, "avhrr_5x5/", sep="")
NT_DIR   = paste(ROOT_DIR, "avhrr_5x5/", sep="")
ET_DIR   = paste(ROOT_DIR, "avhrr_5x5/", sep="")
DT_DIR   = paste(ROOT_DIR, "avhrr_5x5/", sep="")
BT_SUF   = ".B"
NT_SUF   = ".N"
ET_SUF   = ".E"
DT_SUF   = ".D"

KG_DIR   = paste(ROOT_DIR, "kg_5x5/", sep="")
KG_SUF   = ".KG.nc"

LAI_DIR  = paste(ROOT_DIR, "lai_5x5/", sep="")
LAI_SUF  = paste(".RLAI", year, sep="")

PREC_DIR = paste(ROOT_DIR, "wc_5x5/", sep="")
TAVG_DIR = paste(ROOT_DIR, "wc_5x5/", sep="")
TMAX_DIR = paste(ROOT_DIR, "wc_5x5/", sep="")
TMIN_DIR = paste(ROOT_DIR, "wc_5x5/", sep="")
PREC_SUF = ".PREC"
TAVG_SUF = ".TAVG"
TMAX_SUF = ".TMAX"
TMIN_SUF = ".TMIN"

# get regions paras from input file
reg=read.csv(file="reg_5x5", sep=' ', header=F)

# define resolution
xydim    = 1200
xydim1   = 600

days     = seq(1, 46, 1)
mons     = seq(1, 12, 1)
laidata  = array(0, c(xydim,  xydim,  length(days)))
precdata = array(0, c(xydim1, xydim1, length(mons)))
tavgdata = array(0, c(xydim1, xydim1, length(mons)))
tmaxdata = array(0, c(xydim1, xydim1, length(mons)))
tmindata = array(0, c(xydim1, xydim1, length(mons)))

# process the regions one by one
for (i in 1:dim(reg)[1]) {

  cat("\n")
  regname = paste(reg[i,1],reg[i,2],reg[i,3],reg[i,4],sep="_")
  print(paste("Start to precess region: ", regname, sep=""))

  # define dimensions
  dll   = (reg[i,4]-reg[i,2])/xydim
  dll1  = (reg[i,4]-reg[i,2])/xydim1
  lons  = reg[i,2] + c(1:xydim) *dll  - dll/2
  lats  = reg[i,1] - c(1:xydim) *dll  + dll/2
  lons1 = reg[i,2] + c(1:xydim1)*dll1 - dll1/2
  lats1 = reg[i,1] - c(1:xydim1)*dll1 + dll1/2

  londim  <- ncdim_def("lon",  "degrees_east",  as.single(lons),  longname="Longitude"    )
  latdim  <- ncdim_def("lat",  "degrees_north", as.single(lats),  longname="Latitude"     )
  londim1 <- ncdim_def("lon1", "degrees_east",  as.single(lons1), longname="Longitude"    )
  latdim1 <- ncdim_def("lat1", "degrees_north", as.single(lats1), longname="Latitude"     )
  daydim  <- ncdim_def("day",  "day",           as.integer(days), longname="8-day of year")
  mondim  <- ncdim_def("mon",  "month",         as.integer(mons), longname="month of year")

  # define variables
  # --------------------------------------------------

  # land cover data
  fillvalue <- 255
  dlname <- "MODIS Land Cover Type (LC_Type1) data product, MCD12Q1 V006"
  lc <- ncvar_def("LC", "",
      list(londim, latdim),
      fillvalue, dlname, prec="short", compression=6)

  # vcf data
  fillvalue <- 253
  dlname <- "Percent tree cover, MODIS Vegetation Continuous Fields product, MOD43B V006"
  pctt <- ncvar_def("PCTT", "%",
      list(londim, latdim),
      fillvalue, dlname, prec="short", compression=6)

  dlname <- "Percent non-tree cover, MODIS Vegetation Continuous Fields product, MOD43B V006"
  pcth <- ncvar_def("PCTH", "%",
      list(londim, latdim),
      fillvalue, dlname, prec="short", compression=6)

  dlname <- "Percent non-vegetated cover, MODIS Vegetation Continuous Fields product, MOD43B V006"
  pctb <- ncvar_def("PCTB", "%",
      list(londim, latdim),
      fillvalue, dlname, prec="short", compression=6)

  # tree height
  fillvalue <- 255
  dlname <- "Global forest canopy height"
  htop <- ncvar_def("HTOP", "m",
      list(londim1, latdim1),
      fillvalue, dlname, prec="short", compression=6)

  # avhrr data
  fillvalue <- 255
  dlname <- "Broadleaf%, AVHRR Tree Cover Continuous Fields"
  bt <- ncvar_def("BT", "%",
      list(londim1, latdim1),
      fillvalue, dlname, prec="short", compression=6)

  dlname <- "Needleleaf%, AVHRR Tree Cover Continuous Fields"
  nt <- ncvar_def("NT", "%",
      list(londim1, latdim1),
      fillvalue, dlname, prec="short", compression=6)

  dlname <- "Evergreen%, AVHRR Tree Cover Continuous Fields"
  et <- ncvar_def("ET", "%",
      list(londim1, latdim1),
      fillvalue, dlname, prec="short", compression=6)

  dlname <- "Deciduous%, AVHRR Tree Cover Continuous Fields"
  dt <- ncvar_def("DT", "%",
      list(londim1, latdim1),
      fillvalue, dlname, prec="short", compression=6)

  # Köppen-Geiger data
  fillvalue <- 256
  dlname <- "Present Köppen-Geiger climate classification maps at 1-km resolution"
  kg <- ncvar_def("KG", "",
      list(londim1, latdim1),
      fillvalue, dlname, prec="short", compression=6)

  # lai data
  fillvalue <- 255
  dlname <- "Reprocessed MODIS Version 6.1 leaf area index, MCD15A2H V061"
  lai <- ncvar_def("LAI", "m2/m2",
      list(londim, latdim, daydim),
      fillvalue, dlname, prec="short", compression=6)

  # Worldclim 2: New 1-km spatial resolution climate surfaces for global land areas
  fillvalue <- -32768
  dlname <- "Precipitation, Worldclim 2"
  prec <- ncvar_def("PREC", "mm",
      list(londim1, latdim1, mondim),
      fillvalue, dlname, prec="short", compression=6)

  # Worldclim 2: New 1-km spatial resolution climate surfaces for global land areas
  fillvalue <- -3.4e+38
  dlname <- "Monthly mean temperature, Worldclim 2"
  tavg <- ncvar_def("TAVG", "Degree",
      list(londim1, latdim1, mondim),
      fillvalue, dlname, prec="float", compression=6)

  # Worldclim 2: New 1-km spatial resolution climate surfaces for global land areas
  fillvalue <- -3.4e+38
  dlname <- "Monthly maximum temperature, Worldclim 2"
  tmax <- ncvar_def("TMAX", "Degree",
      list(londim1, latdim1, mondim),
      fillvalue, dlname, prec="float", compression=6)

  # Worldclim 2: New 1-km spatial resolution climate surfaces for global land areas
  fillvalue <- -3.4e+38
  dlname <- "Monthly minimum temperature, Worldclim 2"
  tmin <- ncvar_def("TMIN", "Degree",
      list(londim1, latdim1, mondim),
      fillvalue, dlname, prec="float", compression=6)

  # create netCDF file
  # --------------------------------------------------
  filename = paste(RAW_DIR, "RG_", regname, ".RAW", year, ".nc", sep="")
  cmd = paste("rm -f ", filename, sep="")
  system(cmd)
  ncout <- nc_create(filename,
      list(lc, pctt, pcth, pctb, htop, bt, nt, et, dt,
          kg, lai, prec, tavg, tmax, tmin),
      force_v4=TRUE, verbose=FALSE)


  # read data
  # --------------------------------------------------

  # land cover
  filename = paste(LC_DIR, 'RG_', regname, LC_SUF, sep="")
  lcdata   = readBin(filename, integer(), xydim*xydim, size=1, signed=F)
  lcdata   = matrix(lcdata, xydim, xydim)

  # vcf data
  filename = paste(PCTT_DIR, 'RG_', regname, PCTT_SUF, sep="")
  pcttdata = readBin(filename, integer(), xydim*xydim, size=1, signed=F)
  pcttdata = matrix(pcttdata, xydim, xydim)

  filename = paste(PCTH_DIR, 'RG_', regname, PCTH_SUF, sep="")
  pcthdata = readBin(filename, integer(), xydim*xydim, size=1, signed=F)
  pcthdata = matrix(pcthdata, xydim, xydim)

  filename = paste(PCTB_DIR, 'RG_', regname, PCTB_SUF, sep="")
  pctbdata = readBin(filename, integer(), xydim*xydim, size=1, signed=F)
  pctbdata = matrix(pctbdata, xydim, xydim)

  # tree top data
  filename = paste(HTOP_DIR, 'RG_', regname, HTOP_SUF, sep="")
  htopdata = readBin(filename, integer(), xydim1*xydim1, size=1, signed=F)
  htopdata = matrix(htopdata, xydim1, xydim1)

  # avhrr data
  filename = paste(BT_DIR, 'RG_', regname, BT_SUF, sep="")
  btdata   = readBin(filename, integer(), xydim1*xydim1, size=1, signed=F)
  btdata   = matrix(btdata, xydim1, xydim1)

  filename = paste(NT_DIR, 'RG_', regname, NT_SUF, sep="")
  ntdata   = readBin(filename, integer(), xydim1*xydim1, size=1, signed=F)
  ntdata   = matrix(ntdata, xydim1, xydim1)

  filename = paste(ET_DIR, 'RG_', regname, ET_SUF, sep="")
  etdata   = readBin(filename, integer(), xydim1*xydim1, size=1, signed=F)
  etdata   = matrix(etdata, xydim1, xydim1)

  filename = paste(DT_DIR, 'RG_', regname, DT_SUF, sep="")
  dtdata   = readBin(filename, integer(), xydim1*xydim1, size=1, signed=F)
  dtdata   = matrix(dtdata, xydim1, xydim1)

  # kg data
  filename = paste(KG_DIR, 'RG_', regname, KG_SUF, sep="")
  ncfile   = nc_open(filename)
  kgdata   = ncvar_get(ncfile, "zonecode")
  nc_close(ncfile)

  # lai data
  for (iday in 1:length(days)) {
    filename = paste(LAI_DIR, 'RG_', regname, LAI_SUF, sprintf("%0.3d", (iday-1)*8+1), sep="")
    tmpdata  = readBin(filename, integer(), xydim*xydim, size=1, signed=F)
    tmpdata  = matrix(tmpdata, xydim, xydim)
    laidata[,,iday] = tmpdata
  }

  # wc data
  for (imon in 1:length(mons)) {
    filename = paste(PREC_DIR, 'RG_', regname, PREC_SUF, sprintf("%0.2d", imon), ".nc", sep="")
    ncfile   = nc_open(filename)
    tmpdata  = ncvar_get(ncfile, "prec")
    precdata[,,imon] = tmpdata
    nc_close(ncfile)
  }

  for (imon in 1:length(mons)) {
    filename = paste(TAVG_DIR, 'RG_', regname, TAVG_SUF, sprintf("%0.2d", imon), ".nc", sep="")
    ncfile   = nc_open(filename)
    tmpdata  = ncvar_get(ncfile, "tavg")
    tavgdata[,,imon] = tmpdata
    nc_close(ncfile)
  }

  for (imon in 1:length(mons)) {
    filename = paste(TMAX_DIR, 'RG_', regname, TMAX_SUF, sprintf("%0.2d", imon), ".nc", sep="")
    ncfile   = nc_open(filename)
    tmpdata  = ncvar_get(ncfile, "tmax")
    tmaxdata[,,imon] = tmpdata
    nc_close(ncfile)
  }

  for (imon in 1:length(mons)) {
    filename = paste(TMIN_DIR, 'RG_', regname, TMIN_SUF, sprintf("%0.2d", imon), ".nc", sep="")
    ncfile   = nc_open(filename)
    tmpdata  = ncvar_get(ncfile, "tmin")
    tmindata[,,imon] = tmpdata
    nc_close(ncfile)
  }

  # put variables
  # --------------------------------------------------
  ncvar_put(ncout, lc,   lcdata  )
  ncvar_put(ncout, pctt, pcttdata)
  ncvar_put(ncout, pcth, pcthdata)
  ncvar_put(ncout, pctb, pctbdata)
  ncvar_put(ncout, htop, htopdata)
  ncvar_put(ncout, bt,   btdata  )
  ncvar_put(ncout, nt,   ntdata  )
  ncvar_put(ncout, et,   etdata  )
  ncvar_put(ncout, dt,   dtdata  )
  ncvar_put(ncout, kg,   kgdata  )
  ncvar_put(ncout, lai,  laidata )
  ncvar_put(ncout, prec, precdata)
  ncvar_put(ncout, tavg, tavgdata)
  ncvar_put(ncout, tmax, tmaxdata)
  ncvar_put(ncout, tmin, tmindata)

  # add global attributes
  # --------------------------------------------------
  ncatt_put(ncout, 0, "Title", "Land surface model input data: vegetation raw data")
  ncatt_put(ncout, 0, "Authors", "Yuan et al.")
  ncatt_put(ncout, 0, "Address", "School of Atmospheric Sciences, Sun Yat-sen University, Guangzhou, China")
  ncatt_put(ncout, 0, "Email", "yuanh25@mail.sysu.edu.cn")

  # close file
  nc_close(ncout)
}
