
# ======================================================
# Make high-resolution land surface input dataset
# 
# History:
#   2019/06: Hua Yuan, initial version
# ======================================================

library(ncdf4)

# Table to convert 8-day LAI -> monthly data
# Method: weighted average
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

days  = seq(1, 46, 1)
mons  = seq(1, 12, 1)
dom   = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# PFT names
pftnum = seq(0, 15, 1)
pftname = c(
"not_vegetated                           ", # 0
"needleleaf_evergreen_temperate_tree     ", # 1
"needleleaf_evergreen_boreal_tree        ", # 2 
"needleleaf_deciduous_boreal_tree        ", # 3
"broadleaf_evergreen_tropical_tree       ", # 4
"broadleaf_evergreen_temperate_tree      ", # 5
"broadleaf_deciduous_tropical_tree       ", # 6 
"broadleaf_deciduous_temperate_tree      ", # 7 
"broadleaf_deciduous_boreal_tree         ", # 8
"broadleaf_evergreen_temperate_shrub     ", # 9
"broadleaf_deciduous_temperate_shrub     ", # 10
"broadleaf_deciduous_boreal_shrub        ", # 11
"c3_arctic_grass                         ", # 12
"c3_non-arctic_grass                     ", # 13
"c4_grass                                ", # 14
"c3_crop                                 ") # 15

# values in the table below are from Lawrence et al., 2007
# with modifications according to Sitch et al., 2003
laimax = c(0, 5, 5, 5, 7, 7, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4)
tbase  = c(0, 0, 0, 2, 0, 0, 0, 5, 5, 0, 5, 5, 5, 5, 5, 5)
sgdd   = c(0, 0, 0, 100, 0, 0, 0, 200, 200, 0, 100, 100, 100, 100, 100, 100)
minfr  = c(0, .7, .7, 0, .8, .8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
saimin = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .1)
sairtn = c(0, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, 0.)

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
  lcdata   = ncvar_get(fraw, "LC")
  pcttdata = ncvar_get(fraw, "PCTT")
  pcthdata = ncvar_get(fraw, "PCTH")
  pctbdata = ncvar_get(fraw, "PCTB")
  btdata   = ncvar_get(fraw, "BT")
  ntdata   = ncvar_get(fraw, "NT")
  etdata   = ncvar_get(fraw, "ET")
  dtdata   = ncvar_get(fraw, "DT")
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
  uadata[is.na(uadata)] = 0.

  cat("\n")
  print(paste("Start to precess region: ", 
          regname[ireg,1], sep=""))

  # loop for each small 500m grid
  for (i in 1:xydim) {
    cat(i);cat(".")
    for (j in 1:xydim) {
      #for (i in 125:125) {
      #for (j in 30:30) {
      
      # get data
      # ------------------------------
      lc   = lcdata[j, i]
      pctt = pcttdata[j, i]
      pcth = pcthdata[j, i]
      pctb = pctbdata[j, i]
      lai  = laidata[j, i, ] * 0.1

      # NOTE: 1km, 600 x 600
      j1   = floor((j+1)/2)
      i1   = floor((i+1)/2)
      bt   = btdata[j1, i1]
      nt   = ntdata[j1, i1]
      et   = etdata[j1, i1]
      dt   = dtdata[j1, i1]
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

      # ------------------------------------------------------------
      # set crop/urban/water/glacier/wetland(CoLM) pecent
      # basic steps/rules: 
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # step 1: adjust due to water%, total soil reduced  -> (1-water%)
      #
      # step 2: adjust due to ice% (additonal), and crop%
      #   ice% from bare%
      #   crop% from %grass
      #
      # step 3: adjust due to urban% and wetland%
      #   urban and wetland have the same vegeta comp as natrural ones
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #
      # Further update: 
      # the values above need to be updated by higher revolution
      # data. 
      # ------------------------------------------------------------
      
# yuan, 1/1/2020: do not calculate ocean points
      # not land, supposed to be ocean
      if (kg == 0) {   
        lcdata[j,i] = 0
        pocean[j,i] = 100. 
        next
      }

# yuan, 1/2/2020: set NA data to ocean
# inconsistant with soil data
# right now only process RG_25_150_25_155
# need to deal with other regions
# !!!NOTE: 需要重新run这个程序，使与land cover type保持一致
      # set NA data to barren
      if (is.na(lc)) {  
        #lc = 16
        #lcdata[j,i] = 16
        lcdata[j,i] = 0
        pocean[j,i] = 100. 
        next
      }
      
# yuan, 1/1/2020: move it to the first
      # not land, supposed to be ocean
      #if (kg == 0) {   
      #  lcdata[j,i] = 0
      #  pocean[j,i] = 100. 
      #  next
      #}

      if (lc == 11) {           # wetland ==> all
        pwetland[j, i] = 100.
      } else if (lc == 12) {    # crop ==> grass
        pcrop[j, i] = 100.
      } else if (lc == 13) {    # urban ==> all
        purban[j, i] = 100.
      } else if (lc == 14) {    # crop ==> grass
        pcrop[j, i] = 50.
      } else if (lc == 15) {    # ice ==> bare? yes
        pice[j, i] = 100.
        next
      } else if (lc == 17) {    # water ==> no
        pwater[j, i] = 100.
        next
      } 

      # set vegetation fraction
      # ------------------------------

      # calcualte tree and herb fractional cover
      tree_cover = 0.
      herb_cover = 0.

      sum = pctt+pcth+pctb
      if (is.na(sum) || sum==0) {
        tree_cover = tbe + tbd + tne + tnd + ua*0.05
        herb_cover = sbe + sbd + sne + ng + mg +
                     ua*0.15
      } else {
        tree_cover = pctt/0.8
        herb_cover = pcth/0.8

        sum = tree_cover+herb_cover
        if (sum > 100.) {
          tree_cover = tree_cover * 100./sum
          herb_cover = herb_cover * 100./sum
        }
      }

      #browser()
      if (tree_cover+herb_cover == 0) {
        ppft[j,i,1] = 100.
        pcrop[j,i]  = 0.
        next
      }

      # bare soil fractional cover
      ppft[j,i,1] = max(0., 100. - tree_cover - herb_cover)

      # split tree
      if (tree_cover > 0.) {
        
        # tropical 
        if (kg <=4 || kg==6) {

          # assumed no needleleaf tree here
          bt = bt + nt                   
          nt = 0.

          if (et+dt > 0) {
            bet = tree_cover*et/(et+dt)
            bdt = tree_cover*dt/(et+dt)
          } else {
            bet = tree_cover*0.5
            bdt = tree_cover*0.5
          }

          # dominated by evergreen broadleaf
          # Min fradctional 80% of tree cover
          # type, the same below
          if (lc == 2) {
            bet = max(tree_cover*0.8, bet)
            bdt = tree_cover-bet
          }

          # dominated by deciduous broadleaf
          if (lc == 4) {
            bdt = max(tree_cover*0.8, bdt)
            bet = tree_cover-bdt
          }

          ppft[j,i,4+1] = bet
          ppft[j,i,6+1] = bdt

        } else if (kg <=16 ) { # temperate

          # assume no deceduous needleleaf forest
          if (dt<=bt && nt<=et) {
            bdt = dt
            net = nt
            bet = bt-dt
          } else {
            #print("Deceduous needleleaf in temperate, remove it!")
            bdt = bt + (dt-bt)/2.
            net = et + (nt-et)/2.
            bet = 0.
          }

          sum = sum(bdt + net + bet)
          if (sum > 0.) {
            bdt = tree_cover * bdt/sum
            net = tree_cover * net/sum
            bet = tree_cover * bet/sum
          } else {
            bdt = tree_cover * 1/2.
            net = tree_cover * 1/2.
            bet = 0.
          }

          # dominated by evergreen needleleaf
          if (lc == 1) {
            net = max(tree_cover*0.8, net)
            sum = bdt + bet
            if (sum > 0.) {
              bdt = (tree_cover-net)*bdt/sum
              bet = (tree_cover-net)*bet/sum
            }
          }

          # dominated by evergreen broadleaf
          if (lc == 2) {
            bet = max(tree_cover*0.8, bet)
            sum = bdt + net
            if (sum > 0.) {
              bdt = (tree_cover-bet)*bdt/sum
              net = (tree_cover-bet)*net/sum
            }
          }

          # dominated by deciduous broadleaf
          if (lc == 4) {
            bdt = max(tree_cover*0.8, bdt)
            sum = bet + net
            if (sum > 0.) {
              bet = (tree_cover-bdt)*bet/sum
              net = (tree_cover-bdt)*net/sum
            }
          }

          ppft[j,i,1+1] = net
          ppft[j,i,5+1] = bet
          ppft[j,i,7+1] = bdt

        } else if (kg>=17 && kg<=30) { # boreal and polar

          # assume no evergreen broadleaf forest
          if (bt<=dt && et<=nt) {
            bdt = bt
            net = et
            ndt = dt - bt
          } else {
            #print("Evergreen broadleaf in bereal, remove it!")
            bdt = dt + (bt-dt)/2.
            net = nt + (et-nt)/2.
            ndt = 0.
          }

          sum = sum(bdt + net + ndt)
          if (sum > 0.) {
            bdt = tree_cover * bdt/sum
            net = tree_cover * net/sum
            ndt = tree_cover * ndt/sum
          } else {
            bdt = tree_cover * 1/2.
            net = tree_cover * 1/2.
            ndt = 0.
          }

          # dominated by evergreen needleleaf
          # min fradctional 80% of tree cover
          if (lc == 1) {
            net = max(tree_cover*0.8, net)
            sum = bdt + ndt
            if (sum > 0.) {
              bdt = (tree_cover-net)*bdt/sum
              ndt = (tree_cover-net)*ndt/sum
            }
          }

          # dominated by deciduous needelleaf
          if (lc == 3) {
            ndt = max(tree_cover*0.8, ndt)
            sum = bdt + net
            if (sum > 0.) {
              bdt = (tree_cover-ndt)*bdt/sum
              net = (tree_cover-ndt)*net/sum
            }
          }

          # dominated by deciduous broadleaf
          if (lc == 4) {
            bdt = max(tree_cover*0.8, bdt)
            sum = ndt + net
            if (sum > 0.) {
              ndt = (tree_cover-bdt)*ndt/sum
              net = (tree_cover-bdt)*net/sum
            }
          }

          ppft[j,i,2+1] = net
          ppft[j,i,3+1] = ndt
          ppft[j,i,8+1] = bdt
        }

      } else {
        ppft[j,i,c(1:8)+1] = 0.
      }

      #browser()
      # split herbaceous using cci data
      if (herb_cover > 0.) {

        if (kg < 17) { # tropical or temperate

          sbe = sbe + sne
          ng  = ng + mg
          sum = sbe + sbd + ng

          if (sum > 0) {
            sbe = herb_cover * sbe/sum
            sbd = herb_cover * sbd/sum
            ng  = herb_cover * ng /sum
          } else {
            sbe = herb_cover * 1/4.
            sbd = herb_cover * 1/4.
            ng  = herb_cover * 1/2.
          }

          ppft[j,i,9+1]  = sbe
          ppft[j,i,10+1] = sbd
          ppft[j,i,13+1] = ng

        } else if (kg < 29) { # boreal

          sbd = sbd + sbe + sne
          ng  = ng + mg
          sum = sbd + ng

          if (sum > 0) {
            sbd = herb_cover * sbd/sum
            ng  = herb_cover * ng /sum
          } else {
            sbd = herb_cover * 1/2.
            ng  = herb_cover * 1/2.
          }

          ppft[j,i,11+1] = sbd
          ppft[j,i,13+1] = ng

        } else { # polar
          ppft[j,i,12+1] = herb_cover
        }

      } else {
        ppft[j,i,c(9:15)+1] = 0. 
      }

      pctpft = ppft[j,i,]

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
        # ------------------------------

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

      # distribution LAI
      # sum(lai*PFT%) = MODIS_LAI

      # calcualte initial LAI
      laiini[,] = 0.
      pctpft = ppft[j,i,]/100.

      for (imonth in 1:12) {
        sumwgt = sum(phi[,imonth]*laimax*pctpft)
        if (sumwgt > 0.) {
          laiini[,imonth] = phi[,imonth]*laimax/sumwgt*laitot[imonth] 
        } else {
          laiini[,imonth] = 0.
        }
      }

      # adjust LAI 
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
          if (sumnon > 0.) { # don't not consider phenology
            laiini[index,imonth] =
            laimax[index] / sumnon * laitot_nonevg
          } else { # no percentage cover
            laiini[index,imonth] = 0.
          }
        }
      }

      # max LAI value limited
      if (sum(laiini>10) > 0) {
        #browser()
        laiini[laiini>10.] = 10
      }

      # split C3/C4 nature grass (Still et al., 2003)
      # ------------------------------

      # C4 case
      if (sum(tavg>22.)==12 && sum(prec>25.)==12) {
        ppft[j,i,14+1] = sum(ppft[j,i,c(12:13)+1])
        ppft[j,i,c(12:13)+1] = 0.
      } else if (herb_cover>75 && sum(laiini[13+1,]>0.)) { # mixed C3/C4 case 
        frac_c4 = 
        sum(laiini[13+1, tavg>22. & prec>25.]) / sum(laiini[13+1,])

        ppft[j,i,14+1] = ppft[j,i,13+1] * frac_c4
        ppft[j,i,13+1] = ppft[j,i,13+1] * (1.-frac_c4)
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
          if (abs(sum(saiini-saiini_)) > 1.e-6) {
            imonth = 1
            saiini_ = saiini
          }
        }
      }

      # adjust herb to cropland
      # ----------------------------

      # crop exist
      if (pcrop[j,i] > 0) {
        # if herb frac < crop, set it all to crop
        if (herb_cover <= pcrop[j, i]) {
          ppft[j,i,16] = herb_cover
          ppft[j,i,10:15] = 0.
          laiini[10:15] = 0.
          saiini[10:15] = 0.
          pcrop[j,i] = herb_cover
        } else { 
          # if herb > crop, remove the crop
          # fraction from herb
          ppft[j,i,16] = pcrop[j,i]
          ppft[j,i,10:15] = ppft[j,i,10:15]*
          (herb_cover-pcrop[j,i])/herb_cover
        }
      }

      #browser()
      # use adjusted PCT for crop and C3/C4
      pctpft = ppft[j,i,]/100.
      laiini[pctpft<1e-6] = 0.
      saiini[pctpft<1e-6] = 0.

      if (sum(laiini>20) > 0) {
        browser()
      }

      # max SAI value limited
      if (sum(saiini>3) > 0) {
        #browser()
        saiini[saiini>3.] = 3.
      }
      
      if (abs(sum(pctpft)-1.) > 1e-5) {
        print("Sum of area is not equle to 100%! STOP!") 
        browser()
      }
      
      if (abs(pcrop[j,i]-ppft[j,i,16]) > 1e-3) {
        print("Crop area is not conserved! STOP!") 
        browser()
      }

      # output data
      pftlai[j,i,,] = laiini
      pftsai[j,i,,] = saiini
      lclai[j,i,]   = apply(pctpft*laiini, 2, sum)
      lcsai[j,i,]   = apply(pctpft*saiini, 2, sum)
      if(max(abs(lclai[j,i,]-laitot))>1) {
        print("LAI not conserved, set it to laitot")
        print(laitot)
        print(lclai[j,i,])
        cat(lc, pctpft)
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
  dlname <- "MODIS Land Cover Type (LC_Type1) data product, MCD12Q1 V006"
  modis_igbp <- ncvar_def("LC", "-", 
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
  filename = paste(SRF_DIR, "RG_", regname[ireg,1], ".MOD2005.nc", sep="")
  print(filename)
  cmd = paste("rm -f ", filename, sep="")
  system(cmd)
  ncout <- nc_create(filename, 
      list(modis_igbp, monthly_lc_lai, monthly_lc_sai,
          monthly_lai, monthly_sai, pct_crop, pct_urban,
          pct_wetland, pct_glacier, pct_water, pct_ocean, 
          pct_pft, htop), 
          force_v4=TRUE, verbose=F)

  # put variables
  # --------------------------------------------------
  ncvar_put(ncout, modis_igbp,     lcdata)
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
