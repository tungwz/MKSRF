
# read the global region file
reg=read.csv(file="reg_5x5", sep=' ', header=F)

# read binary data of broadleaf tree
raw=readBin("Forest_Height.bin",integer(), 43200*21600, size=1, signed=F)
mraw = matrix(raw, 43200, 21600)

# process the regions one by one
for (i in 1:dim(reg)[1]) {
  regname = paste(reg[i,1],reg[i,2],reg[i,3],reg[i,4],sep="_")
  slon = (reg[i,2]+180)*120+1
  slat = (90-reg[i,1])*120+1
  filename=paste("RG_", regname, '.TREETOP', sep="")
  writeBin(c(mraw[slon:(slon+599),slat:(slat+599)]), filename, size=1)
}
