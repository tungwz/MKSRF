
# read the global region file
reg=read.csv(file="reg_5x5", sep='_', header=F)
regname=read.csv(file="reg_5x5", header=F)

for (year in seq(2000, 2020, 5)) {
  for (day in seq(1, 361, 8)) {

    # read binary data of broadleaf tree
    str_day = sprintf("%03d", day)
    filename = paste("/tera04/yuanhua/modis/global_lai_15s/", "global_lai_15s_", year, "_", str_day, sep="")
    print(filename)
    raw = readBin(filename,integer(), 86400*43200, size=1, signed=F)
    mraw = matrix(raw, 86400, 43200)

    # process the regions one by one
    for (i in 1:dim(reg)[1]) {
      slon = (reg[i,2]+180)*240+1
      slat = (90-reg[i,1])*240+1
      filename=paste("RG_", regname[i,1], '.RLAI', year, str_day, sep="")
      print(filename)
      writeBin(c(mraw[slon:(slon+1199),slat:(slat+1199)]), filename, size=1)
    }
  }
}
