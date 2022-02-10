# ----------------------------------------
# which 5x5deg is land only by using World climate data
# ----------------------------------------

ROOT_DIR  = "/home/yuanhua/tera02/mksrf/"
KG_DIR   = paste(ROOT_DIR, "kg_5x5/", sep="")

# get regions paras from input file
reg=read.csv(file="reg_5x5.full", sep='_', header=F)
regname=read.csv(file="reg_5x5.full", header=F)

# process the regions one by one
for (ireg in 1:dim(reg)[1]) {

  #browser()
  print(ireg)
  filename = paste(KG_DIR, 'RG_', regname[ireg,1], ".KG.nc", sep="")
  fkg      = nc_open(filename)
  zc       = ncvar_get(fkg, "zonecode")

  if (sum(zc!=0) == 0) {
    print(regname[ireg,1])
    cmd = paste("echo '", regname[ireg,1],"' >> noland.txt", sep="")
    system(cmd)
  }

  nc_close(fkg)
}
