
# ======================================================
# make land surface input dataset from ESA CCI data
# 
# History:
#   2019/06: Hua Yuan, initial version
# ======================================================

library(ncdf4)

# tables help to conver 8-day LAI values -> monthly data
# ------------------------------------------------------

idx = matrix(0, 12, 5)
wgt = matrix(0, 12, 5)

# index of 46 8-day's data
idx[1,]  = c(1, 2, 3, 4, NA)
idx[2,]  = c(4, 5, 6, 7, 8)
idx[3,]  = c(8,  9, 10,  11,  12)
idx[4,]  = c(12, 13,  14,  15,  NA)
idx[5,]  = c(16, 17,  18,  19,  NA) 
idx[6,]  = c(19, 20,  21,  22,  23)
idx[7,]  = c(23, 24,  25,  26,  27)
idx[8,]  = c(27, 28,  29,  30,  31)
idx[9,]  = c(31, 32,  33,  34,  35)
idx[10,] = c(35, 36,  37,  38,  NA)
idx[11,] = c(39, 40,  41,  42,  NA)
idx[12,] = c(42, 43,  44,  45,  46)

# weights of 8-day's data
wgt[1,]  = c(8, 8, 8, 7, 0)
wgt[2,]  = c(1, 8, 8, 8, 3)
wgt[3,]  = c(5, 8, 8, 8, 2)
wgt[4,]  = c(6, 8, 8, 8, 0)
wgt[5,]  = c(8, 8, 8, 7, 0)
wgt[6,]  = c(1, 8, 8, 8, 5)
wgt[7,]  = c(3, 8, 8, 8, 4)
wgt[8,]  = c(4, 8, 8, 8, 3)
wgt[9,]  = c(5, 8, 8, 8, 1)
wgt[10,] = c(7, 8, 8, 8, 0)
wgt[11,] = c(8, 8, 8, 6, 0)
wgt[12,] = c(2, 8, 8, 8, 5)

# define resolution
xydim = 1200
xydim1 = 600

days  = seq(1, 46, 1)
mons  = seq(1, 12, 1)
dom   = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# PFT names
pftnum = seq(0, 15, 1)
pftname = c(
"not_vegetated                           ",
"needleleaf_evergreen_temperate_tree     ",
"needleleaf_evergreen_boreal_tree        ",
"needleleaf_deciduous_boreal_tree        ",
"broadleaf_evergreen_tropical_tree       ",
"broadleaf_evergreen_temperate_tree      ",
"broadleaf_deciduous_tropical_tree       ",
"broadleaf_deciduous_temperate_tree      ",
"broadleaf_deciduous_boreal_tree         ",
"broadleaf_evergreen_temperate_shrub     ",
"broadleaf_deciduous_temperate_shrub     ",
"broadleaf_deciduous_boreal_shrub        ",
"c3_arctic_grass                         ",
"c3_non-arctic_grass                     ",
"c4_grass                                ",
"c3_crop                                 ")

# values in the table below are from Lawrence et al., 2007
laimax = c(0, 5, 5, 5, 7, 7, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4)
tbase  = c(0, 0, 0, 2, 0, 0, 0, 5, 5, 0, 5, 5, 5, 5, 5, 5)
sgdd   = c(0, 0, 0, 100, 0, 0, 0, 200, 200, 0, 100, 100, 100, 100, 100, 100)
minfr  = c(0, .7, .7, 0, .8, .8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
saimin = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .1)
sairtn = c(0, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, 0.)

# kopper climate zones
tropical = c(1:4,6)
temperate_warm = c(5,7,8,11,14)
temperate_cool = c(9,10,12,13,15,16)
temperate = c(temperate_warm, temperate_cool)
boreal_warm = c(17,21,25)
boreal_cool = c(18:20,22:24,26:28)
boreal = c(boreal_warm, boreal_cool)
polar = c(29,30)

phi    = array(1, c(length(pftnum), length(mons)))
laiini = array(0, c(length(pftnum), length(mons)))
phimin = rep(0.5, length(mons))

# define output data
lclai  = array(0, c(xydim, xydim, length(mons)))
lcsai  = array(0, c(xydim, xydim, length(mons)))
pftlai = array(0, c(xydim, xydim, length(pftnum), length(mons)))
pftsai = array(0, c(xydim, xydim, length(pftnum), length(mons)))
pcrop    = array(0, c(xydim, xydim))
purban   = array(0, c(xydim, xydim))
pwetland = array(0, c(xydim, xydim))
pice     = array(0, c(xydim, xydim))
pwater   = array(0, c(xydim, xydim))
pocean   = array(0, c(xydim, xydim))
ppft     = array(0, c(xydim, xydim, length(pftnum)))
htop500  = array(0, c(xydim, xydim))

ROOT_DIR  = "/home/yuanhua/hard/mksrf/"
RAW_DIR   = paste(ROOT_DIR, "raw_5x5/", sep="")
CCI_DIR   = paste(ROOT_DIR, "cci_5x5/", sep="")
SRF_DIR   = paste(ROOT_DIR, "srf_5x5/", sep="")

# get regions paras from input file
reg=read.csv(file=REGFILE, sep='_', header=F)
regname=read.csv(file=REGFILE, header=F)

# process the regions one by one
for (ireg in 1:dim(reg)[1]) {
  
  filename = paste(RAW_DIR, 'RG_', regname[ireg,1], ".RAW2005.nc", sep="")
  fraw     = nc_open(filename)
  filename = paste(CCI_DIR, 'RG_', regname[ireg,1], ".CCI2005.nc", sep="")
  fcci     = nc_open(filename)

  # get raw data
  lcdata   = ncvar_get(fcci, "majority_class_1")
  kgdata   = ncvar_get(fraw, "KG")
  laidata  = ncvar_get(fraw, "LAI")
  precdata = ncvar_get(fraw, "PREC")
  tavgdata = ncvar_get(fraw, "TAVG")
  tmaxdata = ncvar_get(fraw, "TMAX")
  tmindata = ncvar_get(fraw, "TMIN")
  htopdata = ncvar_get(fraw, "HTOP")

  tbedata  = ncvar_get(fcci, "Tree_Broadleaf_Evergreen")
  tbddata  = ncvar_get(fcci, "Tree_Broadleaf_Deciduous")
  tnedata  = ncvar_get(fcci, "Tree_Needleleaf_Evergreen")
  tnddata  = ncvar_get(fcci, "Tree_Needleleaf_Deciduous")
  sbedata  = ncvar_get(fcci, "Shrub_Broadleaf_Evergreen")
  sbddata  = ncvar_get(fcci, "Shrub_Broadleaf_Deciduous")
  snedata  = ncvar_get(fcci, "Shrub_Needleleaf_Evergreen")
  ngdata   = ncvar_get(fcci, "Natural_Grass")
  mgdata   = ncvar_get(fcci, "Managed_Grass")
  bsdata   = ncvar_get(fcci, "Bare_Soil")
  wtdata   = ncvar_get(fcci, "Water")
  sidata   = ncvar_get(fcci, "Snow_Ice")
  uadata   = ncvar_get(fcci, "Urban_areas")

  nc_close(fraw)
  nc_close(fcci)

  tbedata[is.na(tbedata)] = 0.
  tbddata[is.na(tbddata)] = 0.
  tnedata[is.na(tnedata)] = 0.
  tnddata[is.na(tnddata)] = 0.
  sbedata[is.na(sbedata)] = 0.
  sbddata[is.na(sbddata)] = 0.
  snedata[is.na(snedata)] = 0.
  ngdata[is.na(ngdata)] = 0.
  mgdata[is.na(mgdata)] = 0.
  bsdata[is.na(bsdata)] = 0.
  wtdata[is.na(wtdata)] = 0.
  sidata[is.na(sidata)] = 0.
  uadata[is.na(uadata)] = 0.

  cat("\n")
  print(paste("Start to precess region: ", 
          regname[ireg,1], sep=""))
  
  # loop for each small 500m grid
  for (i in 1:xydim) {
    cat(i);cat(".")
    for (j in 1:xydim) {
    #for (j in 1:5) {
      #for (i in 8:8) {
      #for (j in 714:714) {

      # get data
      # ------------------------------
      lc = lcdata[j, i]
      lai = laidata[j, i, ] * 0.1

      # NOTE: 1Km, 600 x 600
      j1   = floor((j+1)/2)
      i1   = floor((i+1)/2)
      kg   = kgdata[j1, i1]
      prec = precdata[j1, i1, ]
      tavg = tavgdata[j1, i1, ]
      tmax = tmaxdata[j1, i1, ]
      tmin = tmindata[j1, i1, ]
      htop500[j,i] = htopdata[j1,i1]

      tbe  = tbedata[j,i]*100
      tbd  = tbddata[j,i]*100
      tne  = tnedata[j,i]*100
      tnd  = tnddata[j,i]*100
      sbe  = sbedata[j,i]*100
      sbd  = sbddata[j,i]*100
      sne  = snedata[j,i]*100
      ng   = ngdata[j,i]*100
      mg   = mgdata[j,i]*100
      bs   = bsdata[j,i]*100
      wt   = wtdata[j,i]*100
      si   = sidata[j,i]*100
      ua   = uadata[j,i]*100

      # initialization
      lclai[j,i,]   = 0.
      lcsai[j,i,]   = 0.
      pftlai[j,i,,] = 0.
      pftsai[j,i,,] = 0.
      pcrop[j,i]    = 0.
      purban[j,i]   = 0.
      pwetland[j,i] = 0.
      pice[j,i]     = 0.
      pwater[j,i]   = 0.
      pocean[j,i]   = 0.
      ppft[j,i,]    = 0.

      # set land cover type
      lc = trunc(lc/10)

      # set no-data to barren 
      if (lc == 0) {lc = 20}
      lcdata[j,i] = lc

      if (kg == 0) {   # ocean case
        lcdata[j,i] = 0
        pocean[j,i] = 100
        next
      }

      # set PFT... fraction
      # ------------------------------
      pcrop[j,i]  = mg
      purban[j,i] = ua*0.95 # exclude the water body
      pice[j,i]   = si
      pwater[j,i] = wt

      # vegetation fraction

      if (tbe+tbd+tne+tnd+sbe+sbd+sne+ng+mg+ua*0.2==0) {
        ppft[j,i,1] = 100.
        next 
      }

      urban_ne = NA
      urban_bd = NA
      urban_ng = NA

      #browser()
      # needleleaf evergreen, temperate
      if (kg %in% c(tropical, temperate, boreal_warm)) { 
        ppft[j,i,2] = tne 
        urban_ne = 2
      }

      # needleleaf evergreen, boreal
      if (kg %in% c(boreal_cool, polar)) {
        ppft[j,i,3] = tne
        urban_ne = 3
      }

      # needleleaf deciduous, boreal
      # bug of Poulter et al., 2011
      # no needleleaf deciduous trees in temperate
      #if (kg %in% c(temperate_cool, boreal_cool)) {
      if (kg %in% c(boreal_cool, polar)) {
        ppft[j,i,4] = tnd
      }

      # boreal evergreen, tropical
      if (kg %in% tropical) {
        ppft[j,i,5] = tbe
      }

      # boreal evergreen, temperate
      if (kg %in% c(temperate, boreal, polar)) {
        ppft[j,i,6] = tbe
      }

      # boreal deciduous, tropical
      if (kg %in% tropical) {
        ppft[j,i,7] = tbd + tnd
        urban_bd = 7
      }

      # boreal deciduous, temperate
      if (kg %in% c(temperate, boreal_warm)) {
        ppft[j,i,8] = tbd + tnd
        urban_bd = 8
      }

      # boreal deciduous, boreal
      if (kg %in% c(boreal_cool, polar)) {
        ppft[j,i,9] = tbd
        urban_bd = 9
      }

      # broadleaf evergreen shrub, temperate
      # broadleaf deciduous shrub, temperate
      if (kg %in% c(tropical,temperate)) {
        ppft[j,i,10] = sbe + sne
        ppft[j,i,11] = sbd
      }

      # broadleaf deciduous shurb, boreal
      if (kg %in% c(boreal, polar)) {
        ppft[j,i,12] = sbe + sne + sbd
      }

      # C3 artic grass
      if (kg %in% polar) {
        ppft[j,i,13] = ng
        urban_ng = 13
      } else {
        if (kg %in% c(tropical, temperate_warm)) {
          ppft[j,i,15] = ng
          urban_ng = 15
        } else {
          ppft[j,i,14] = ng
          urban_ng = 14
        }
      }

      ppft[j,i,16] = mg  # crop
      ppft[j,i,1]  = bs

      pctpft = ppft[j,i,]

      # adjust for urban area
      pctpft[1] = pctpft[1] + ua*0.75
      pwater[j,i] = pwater[j,i] + ua*0.05

      if (!is.na(urban_ne)) {
        pctpft[urban_ne] = pctpft[urban_ne] + ua*0.025
      } 

      if (!is.na(urban_bd)) {
        pctpft[urban_bd] = pctpft[urban_bd] + ua*0.025
      }

      if (!is.na(urban_ng)) {
        pctpft[urban_ng] = pctpft[urban_ng] + ua*0.15
      } else {
        print("Error!")
      }

      # calcualte LAI
      # ------------------------------

      # convert 8-day lai to monthly lai
      # loop for each month
      laitot = c(1:12)

      for ( imonth in c(1:12) ) {
        laitot[imonth] = 
        sum(lai[idx[imonth,]]*wgt[imonth,], na.rm=T)/
        sum(wgt[imonth,])
      }

      #browser()
      # deal with urban area
      # area may with vegetation
      #if (100-wt-si-ua*0.05 > 0) {
      #  laitot  = laitot/(100-wt-si-ua*0.05)*100
      #}

      # calcuate GDD and phi
      phi[,]  = 1.

      if (is.na(sum(tavg))) { # when T and P are missing
        print("NA temperature! Error!")
        browser()
      } else {
        itmin = which.min(tavg)
        mcnt  = 0
        dd2   = c(1:12)*0.
        dd5   = c(1:12)*0.
        gdd2  = c(1:12)*0.
        gdd5  = c(1:12)*0.
        rgdd2  = c(1:12)*0.
        rgdd5  = c(1:12)*0.
      }

      while ( mcnt < 12 ) {

        imonth = (itmin+mcnt-1)%%12 + 1

        # calculate gdd according to tmin, tmax, tavg, 2, 5degree

        # calcualte day index of tavg 
        if (tmax[imonth]-tmin[imonth] > 0) {
          davg = (tmax[imonth]-tavg[imonth])/(tmax[imonth]-tmin[imonth])*dom[imonth]
        } else {
          dvag = 0
        }

        # gdd2 in different cases
        if (2. < tmin[imonth]) {
          dd2[imonth] = (tavg[imonth]-2.) * dom[imonth]
        } else if (2. >= tmax[imonth]) {
          dd2[imonth] = 0.
        } else if (2.>=tmin[imonth] && 2.<tavg[imonth] && tavg[imonth]!=tmin[imonth]) {
          dd2[imonth] = (tavg[imonth]-2.)/2.*
          (tavg[imonth]-2.)/(tavg[imonth]-tmin[imonth])*davg +
          (tmax[imonth]+tavg[imonth]-4.)/2.*(dom[imonth]-davg)
        } else if (2.>=tavg[imonth] && 2<tmax[imonth] && tmax[imonth]!=tavg[imonth]) {
          dd2[imonth] = (tmax[imonth]-2.)/2.*
          (tmax[imonth]-2.)/(tmax[imonth]-tavg[imonth])*(dom[imonth]-davg)
        } else {
          dd2[imonth] = max(0, (tavg[imonth]-2.)*dom[imonth])
        }

        # gdd5 in different cases
        if (5. < tmin[imonth]) {
          dd5[imonth] = (tavg[imonth]-5.) * dom[imonth]
        } else if (5. >= tmax[imonth]) {
          dd5[imonth] = 0.
        } else if (5.>=tmin[imonth] && 5.<tavg[imonth] && tavg[imonth]!=tmin[imonth]) {
          dd5[imonth] = (tavg[imonth]-5.)/2.*
          (tavg[imonth]-5.)/(tavg[imonth]-tmin[imonth])*davg +
          (tmax[imonth]+tavg[imonth]-10)/2.*(dom[imonth]-davg)
        } else if (5.>=tavg[imonth] && 5<tmax[imonth] && tmax[imonth]!=tavg[imonth]) {
          dd5[imonth] = (tmax[imonth]-5.)/2.*
          (tmax[imonth]-5.)/(tmax[imonth]-tavg[imonth])*(dom[imonth]-davg) 
        } else {
          dd5[imonth] = max(0, (tavg[imonth]-5.)*dom[imonth])
        }

        mcnt = mcnt + 1
      }

      # split the tmin month gdd into 2 parts.
      # gdd, rgdd

      if (dd2[itmin] == 0) {
        gdd2[itmin] = 0.
        rgdd2[itmin] = 0.
      } else {
        nextmonth = itmin%%12 + 1 
        prevmonth = (10+itmin)%%12 + 1
        sum = max(0, (tmin[nextmonth]-2.)) + 
        max(0, (tmin[prevmonth]-2.))
        if (sum > 0) {
          gdd2[itmin]  = dd2[itmin]*max(0, tmin[nextmonth]-2.)/sum
          rgdd2[itmin] = dd2[itmin]*max(0, tmin[prevmonth]-2.)/sum
        } else {
          if (tavg[nextmonth]+tavg[prevmonth]-2*tavg[itmin] > 0) {
            gdd2[itmin]  = dd2[itmin]*(tavg[nextmonth]-tavg[itmin])/
            (tavg[nextmonth]+tavg[prevmonth]-2*tavg[itmin])
            rgdd2[itmin] = dd2[itmin]*(tavg[prevmonth]-tavg[itmin])/
            (tavg[nextmonth]+tavg[prevmonth]-2*tavg[itmin])
          } else {
            gdd2[itmin]  = dd2[itmin] * 0.5
            rgdd2[itmin] = dd2[itmin] * 0.5
          }
        }
      }

      if (dd5[itmin] == 0) {
        gdd5[itmin] = 0.
        rgdd5[itmin] = 0.
      } else {
        nextmonth = itmin%%12 + 1 
        prevmonth = (10+itmin)%%12 + 1
        sum = max(0, (tmin[nextmonth]-5.)) + 
        max(0, (tmin[prevmonth]-5.))
        if (sum > 0) {
          gdd5[itmin]  = dd5[itmin]*max(0, tmin[nextmonth]-5.)/sum
          rgdd5[itmin] = dd5[itmin]*max(0, tmin[prevmonth]-5.)/sum
        } else {
          if (tavg[nextmonth]+tavg[prevmonth]-2*tavg[itmin] > 0) {
            gdd5[itmin]  = dd5[itmin]*(tavg[nextmonth]-tavg[itmin])/
            (tavg[nextmonth]+tavg[prevmonth]-2*tavg[itmin])
            rgdd5[itmin] = dd5[itmin]*(tavg[prevmonth]-tavg[itmin])/
            (tavg[nextmonth]+tavg[prevmonth]-2*tavg[itmin])
          } else {
            gdd5[itmin]  = dd5[itmin] * 0.5
            rgdd5[itmin] = dd5[itmin] * 0.5
          }
        }
      }

      imonth_prev = itmin
      mcnt = 1
      while ( mcnt < 12 ) {

        imonth = (itmin+mcnt-1)%%12 + 1

        gdd2[imonth] = gdd2[imonth_prev] + dd2[imonth]
        gdd5[imonth] = gdd5[imonth_prev] + dd5[imonth]

        imonth_prev = imonth
        mcnt = mcnt + 1
      }

      imonth_prev = itmin
      mcnt = 1
      while ( mcnt < 12 ) {

        imonth = (12+itmin-mcnt-1)%%12 + 1

        rgdd2[imonth] = rgdd2[imonth_prev] + dd2[imonth]
        rgdd5[imonth] = rgdd5[imonth_prev] + dd5[imonth]

        imonth_prev = imonth
        mcnt = mcnt + 1
      }

      #browser()
      gdd2[itmin]  = dd2[itmin]
      rgdd2[itmin] = dd2[itmin]
      gdd5[itmin]  = dd5[itmin]
      rgdd5[itmin] = dd5[itmin]

      gdd2 = pmin(gdd2, rgdd2)
      gdd5 = pmin(gdd5, rgdd5)

      # calcualte phi using gdd/sgdd
      phi[4,]  = pmax(phimin, pmin(rep(1,12), gdd2/sgdd[4]))
      phi[8,]  = pmax(phimin, pmin(rep(1,12), gdd5/sgdd[8]))
      phi[9,]  = pmax(phimin, pmin(rep(1,12), gdd5/sgdd[9]))
      phi[11,] = pmax(phimin, pmin(rep(1,12), gdd5/sgdd[11]))
      phi[12,] = pmax(phimin, pmin(rep(1,12), gdd5/sgdd[12]))
      phi[13,] = pmax(phimin, pmin(rep(1,12), gdd5/sgdd[13]))
      phi[14,] = pmax(phimin, pmin(rep(1,12), gdd5/sgdd[14]))
      phi[15,] = pmax(phimin, pmin(rep(1,12), gdd5/sgdd[15]))
      phi[16,] = pmax(phimin, pmin(rep(1,12), gdd5/sgdd[16]))

      #browser()
      # distribution LAI
      # sum(lai*PFT%) = MODIS_LAI

      # calcualte initial LAI
      laiini[,] = 0.
      #pctpft = pctpft/100.
      sumpctpft = sum(pctpft)
      pctpft = pctpft/sumpctpft

      for (imonth in 1:12) {
        sumwgt = sum(phi[,imonth]*laimax*pctpft)
        if (sumwgt > 0.) {
          laiini[,imonth] = 
          phi[,imonth]*laimax/sumwgt*laitot[imonth] 
        } else {
          laiini[,imonth] = 0.
        }
      }

      # adjust LAI according to minimum evergreen fraction
      # ---------------------------
      #laiinimax = apply(laiini, 1, max) # check
      laimaxloc = which.max(laitot)
      laiinimax = laiini[,laimaxloc]

      for (imonth in 1:12) {
        # evergreen tree phenology (Zeng et al., 2002)
        # minimum fraction of maximum initial PFT LAI
        # laiini[,imonth] = pmax(laiini[,imonth], laiinimax*minfr)

        # max value limited
        ind = which(laiini[,imonth] > laiinimax, arr.ind=T)
        laiini[ind,imonth] = laiinimax[ind]

        # calculate for non-evergreen PFT
        index  = c(1,2,4,5,9)+1
        sumevg = 
        sum(laiini[index,imonth]*pctpft[index])
        laitot_nonevg = max(laitot[imonth]-sumevg, 0.)

        index  = c(3,6,7,8,10:15)+1
        sumnon = 
        sum(phi[index,imonth]*laimax[index]*pctpft[index])
        if (sumnon > 0.) {
          laiini[index,imonth] = phi[index,imonth] * 
          laimax[index] / sumnon * laitot_nonevg
        } else {
          sumnon =
          sum(laimax[index]*pctpft[index])
          if (sumnon > 0.) {  # don't not consider phenology
            laiini[index,imonth] =
            laimax[index] / sumnon * laitot_nonevg
          } else {  # no percentage cover
            laiini[index,imonth] = 0.
          }
        }
      }

      #browser()
      # max LAI value limited
      if (sum(laiini>10) > 0) {
        #  browser()
        laiini[laiini>10.] = 10
      }

      # calculate SAI
      # ------------------------------
      saimin_ = saimin*apply(laiini, 1, max)
      saimin_[2:16] = saimin_[2:16]/laimax[2:16]
      saiini  = matrix(0, length(pftnum), length(mons))
      saiini_ = saiini
      saiini_[,] = saimin

      # PFT LAI diff
      laidiff = laiini[,c(12, 1:11)]-laiini[,]
      laidiff[laidiff<0] = 0

      imonth = 1
      while (imonth < 13) {
        saires = sairtn*saiini_[,(imonth+10)%%12+1]
        saiini[,imonth] = pmax(saires+laidiff[,imonth]*0.5, saimin_)
        imonth = imonth + 1
        if (imonth == 13) {
          if (abs(sum(saiini-saiini_)) > 1.e-5) {
            imonth = 1
            saiini_ = saiini
          }
        }
      }

      laiini[pctpft<1e-5] = 0.
      saiini[pctpft<1e-5] = 0.

      if (sum(laiini>20) > 0) {
        browser()
      }

      # max SAI value limited
      if (sum(saiini>3) > 0) {
        #browser()
        saiini[saiini>3.] = 3.
      }

      # output data
      pftlai[j,i,,] = laiini
      pftsai[j,i,,] = saiini
      ppft[j,i,]    = pctpft*sumpctpft

      if (abs(sumpctpft+pice[j,i]+pwater[j,i]-100) > 1e-3) {
        print("Sum of area is not equle to 100%! STOP!") 
        browser()
      }
      
      if (abs(pcrop[j,i]-ppft[j,i,16]) > 1e-3) {
        print("Crop area is not conserved! STOP!") 
        browser()
      }

      sum = sum(pctpft)
      if (sum > 0) {
        lclai[j,i,]   = apply(pctpft*laiini, 2, sum)
        lcsai[j,i,]   = apply(pctpft*saiini, 2, sum)
      } else {
        lclai[j,i,]   = laitot
        lcsai[j,i,]   = apply(pctpft*saiini, 2, sum)
      }
      if(max(abs(lclai[j,i,]-laitot))>1) {
        print("LAI not conserved, set it to laitot")
        print(laitot)
        print(lclai[j,i,])
        cat(lc, wt, pctpft)
        lclai[j,i,] = laitot
        if (abs(sum(pctpft)-1.) > 1e-5) browser()
      }
    }
  }

  #browser()
  dll   = (reg[ireg,4]-reg[ireg,2])/xydim
  lons  = reg[ireg,2] + c(1:xydim)*dll - dll/2
  lats  = reg[ireg,1] - c(1:xydim)*dll + dll/2

  londim  <- ncdim_def("lon", "degrees_east",  as.single(lons), longname="Longitude")
  latdim  <- ncdim_def("lat", "degrees_north", as.single(lats), longname="Latitude")
  pftdim  <- ncdim_def("pft", "PFT index", as.integer(pftnum), longname="Index of PFT")
  mondim  <- ncdim_def("mon", "month", as.integer(mons), longname="Month of year")

  # define variables
  # --------------------------------------------------

  # land cover data
  fillvalue <- 255
  dlname <- "ESA CCI Land Cover Type data"
  esa_cci<- ncvar_def("LC", "-", 
      list(londim, latdim), 
      fillvalue, dlname, prec="short", compression=6)

  fillvalue <- -999.
  dlname <- "Monthly landcover LAI values"
  monthly_lc_lai <- ncvar_def("MONTHLY_LC_LAI", "m^2/m^2", 
      list(londim, latdim, mondim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Monthly landcover SAI values"
  monthly_lc_sai <- ncvar_def("MONTHLY_LC_SAI", "m^2/m^2", 
      list(londim, latdim, mondim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Monthly PFT LAI values"
  monthly_lai <- ncvar_def("MONTHLY_LAI", "m^2/m^2", 
      list(londim, latdim, pftdim, mondim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Monthly PFT SAI values"
  monthly_sai <- ncvar_def("MONTHLY_SAI", "m^2/m^2", 
      list(londim, latdim, pftdim, mondim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Percent crop cover"
  pct_crop <- ncvar_def("PCT_CROP", "%", 
      list(londim, latdim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Percent urban cover"
  pct_urban <- ncvar_def("PCT_URBAN", "%", 
      list(londim, latdim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Percent wetland cover"
  pct_wetland <- ncvar_def("PCT_WETLAND", "%", 
      list(londim, latdim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Percent glacier/ice cover"
  pct_glacier <- ncvar_def("PCT_GLACIER", "%", 
      list(londim, latdim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Percent water body cover"
  pct_water <- ncvar_def("PCT_WATER", "%", 
      list(londim, latdim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Percent ocean cover"
  pct_ocean <- ncvar_def("PCT_OCEAN", "%", 
      list(londim, latdim), 
      fillvalue, dlname, prec="float", compression=6)

  dlname <- "Percent PFT cover"
  pct_pft <- ncvar_def("PCT_PFT", "%", 
      list(londim, latdim, pftdim), 
      fillvalue, dlname, prec="float", compression=6)

  # tree height
  fillvalue <- 255
  dlname <- "Global forest canopy height"
  htop <- ncvar_def("HTOP", "m", 
      list(londim, latdim), 
      fillvalue, dlname, prec="short", compression=6)

  # create netCDF file
  # --------------------------------------------------
  filename = paste(SRF_DIR, "RG_", regname[ireg,1], ".ESA2005.nc", sep="")
  print(filename)
  cmd = paste("rm -f ", filename, sep="")
  system(cmd)
  ncout <- nc_create(filename, 
      list(esa_cci, monthly_lc_lai, monthly_lc_sai,
          monthly_lai, monthly_sai, pct_crop, pct_urban,
          pct_wetland, pct_glacier, pct_water, pct_ocean, 
          pct_pft, htop), 
          force_v4=TRUE, verbose=F)

  # put variables
  # --------------------------------------------------
  ncvar_put(ncout, esa_cci,        lcdata)
  ncvar_put(ncout, monthly_lc_lai, lclai)
  ncvar_put(ncout, monthly_lc_sai, lcsai)
  ncvar_put(ncout, monthly_lai,    pftlai)
  ncvar_put(ncout, monthly_sai,    pftsai)
  ncvar_put(ncout, pct_crop,       pcrop)
  ncvar_put(ncout, pct_urban,      purban)
  ncvar_put(ncout, pct_wetland,    pwetland)
  ncvar_put(ncout, pct_glacier,    pice)
  ncvar_put(ncout, pct_water,      pwater)
  ncvar_put(ncout, pct_ocean,      pocean)
  ncvar_put(ncout, pct_pft,        ppft)
  ncvar_put(ncout, htop,           htop500)

  # add global attributes
  # --------------------------------------------------
  ncatt_put(ncout, 0, "Title", "Land surface model input vagetation data")
  ncatt_put(ncout, 0, "Authors", "Yuan et al.")
  ncatt_put(ncout, 0, "Address", "School of Atmospheric Sciences, Sun Yat-sen University, Guangzhou, China")
  ncatt_put(ncout, 0, "Email", "yuanh25@mail.sysu.edu.cn")

  # close file
  nc_close(ncout)
}
