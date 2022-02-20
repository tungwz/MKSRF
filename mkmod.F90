PROGRAM mkmod

 use netcdf

 IMPLICIT NONE

 INTEGER , PARAMETER :: r8 = selected_real_kind(12)
 INTEGER , PARAMETER :: xydim = 1200, pfts = 16, mon = 12
 INTEGER , PARAMETER :: day = 46

 CHARACTER (len=*), PARAMETER :: Title   = "Land surface model input vagetation data"
 CHARACTER (len=*), PARAMETER :: Authors = "Yuan et al."
 CHARACTER (len=*), PARAMETER :: Address = "School of Atmospheric Sciences, Sun Yat-sen University, Guangzhou, China"
 CHARACTER (len=*), PARAMETER :: Email   = "yuanh25@mail.sysu.edu.cn"

 INTEGER , DIMENSION(12,5) :: idx, wgt
 INTEGER , DIMENSION(46)   :: days
 INTEGER , DIMENSION(12)   :: mons, dom
 INTEGER , DIMENSION(16)   :: pftnum

 REAL(r8), DIMENSION(12)   :: phimin  !=0.5
 REAL(r8), DIMENSION(16)   :: laimax, tbcase, sgdd, &
                              minfr, saimin, sairtn, saimin1, saires, &
                              pctpft, laiinimax
 REAL(r8), DIMENSION(16,12):: phi, laiini, laiini1, saiini2

 REAL(r8), DIMENSION(1200) :: lats, lons

 REAL(r8), DIMENSION(1200,1200)       :: pcrop, purban, pwetland, &
                                         pice, pwater, pocean, tbedata, &
                                         tbddata, tnedata, tnddata, sbedata, sbddata, ngdata, &
                                         mgdata, uadata, snedata
 REAL(r8), DIMENSION(1200,1200,12)    :: lclai, lcsai
 REAL(r8), DIMENSION(1200,1200,16)    :: ppft
 REAL(r8), DIMENSION(1200,1200,16,12) :: pftlai, pftsai

 CHARACTER (len=255) :: ROOT_DIR = "/home/yuanhua/hard/mksrf/"
 CHARACTER (len=255) :: RAW_DIR  = "raw_5x5/"
 CHARACTER (len=255) :: CCI_DIR  = "cci_5x5/"
 CHARACTER (len=255) :: SRF_DIR  = "srf_5x5/"
 CHARACTER (len=255) :: OUT_DIR = "/hard/dongwz/test/"

 CHARACTER (len=255) :: REGFILE  = "reg_5x5"
 CHARACTER (len=4)   :: year = "2005"
 CHARACTER (len=255) :: filename

 CHARACTER (len=255)  , DIMENSION(16) :: pftname
 CHARACTER (len=4)    , DIMENSION(4)  :: creg

 ! input vars

 INTEGER        , DIMENSION(4)            :: reg
 INTEGER(kind=2), DIMENSION(1200,1200)    :: lcdata, htop500
 REAL(r8)       , DIMENSION(1200,1200)    :: pcttdata, pcthdata, pctbdata
 REAL(r8)       , DIMENSION(600 , 600)    :: btdata, ntdata, etdata, &
                                             kgdata, dtdata, htopdata
 REAL(r8)       , DIMENSION(1200,1200,46) :: laidata
 REAL(r8)       , DIMENSION(600,600,12)   :: precdata, tavgdata, tmaxdata, &
                                             tmindata

 ! input vars id
 INTEGER :: ncid, lcid, pcttid, pcthid, pctbid, btid, ntid, &
            etid, dtid, kgid, laiid, precid, tavgid, tmaxid, tminid, htopid
 INTEGER :: tbeid, tbdid, tneid, tndid, sbeid, sbdid, sneid, &
            ngid, mgid, uaid
 ! output vars id
 INTEGER :: clai_id, crop_id, csai_id, dims, gice_id, htop_id, &
            lat_dimid, lat_vid, lc_id, lon_dimid, lon_vid, mon_dimid, &
            mon_vid, ocea_id, pft_dimid, pft_vid, plai_id, ppft_id, &
            psai_id, slai_id, urbn_id, watr_id, wetl_id

 ! vars
 INTEGER(kind=2) :: shortvalue
 INTEGER         :: loc1(1), loc2(1)
 INTEGER         :: i, j, i1, j1, imonth, itmin, nextmonth, prevmonth, &
                     imonth_prev, laimaxloc, iloop, inx, inx1, ipft, ireg, mcnt

 INTEGER , DIMENSION(5)     :: indx1=(/2,3,5,6,10/)
 INTEGER , DIMENSION(10)    :: indx2=(/4,7,8,9,11,12,13,14,15,16/)
 REAL                       :: fillvalue
 REAL(r8), DIMENSION(46)    :: lai
 REAL(r8), DIMENSION(12)    :: laitot, dd2, dd5, gdd2, gdd5, &
                               rgdd2, rgdd5, tmax, tmin, tavg, prec
 REAL(r8), DIMENSION(16,12) :: saiini, saiini1, laidiff

 REAL(r8) :: lc, pctt, pcth, pctb, bt, nt, et, dt, kg, &
             tbe, tbd, tne, tnd, sbe, sbd, sne, ng, mg, ua
 REAL(r8) :: bdt, bet, dll, laitot_nonevg, laiup, ndt, net
 REAL(r8) :: tree_cover, herb_cover, summ, davg, sumwgt, &
             sumevg, sumnin, frac_c4, sumnon, x1, x2, sum_judg
 INTEGER  :: XY2D(2), XY3D(3), XY4D(4), XY3F(3)

 ! index of 46 8-day's data
 idx(1,:)  = (/ 1,  2,  3,  4,  0/)
 idx(2,:)  = (/ 4,  5,  6,  7,  8/)
 idx(3,:)  = (/ 8,  9, 10, 11, 12/)
 idx(4,:)  = (/12, 13, 14, 15,  0/)
 idx(5,:)  = (/16, 17, 18, 19,  0/)
 idx(6,:)  = (/19, 20, 21, 22, 23/)
 idx(7,:)  = (/23, 24, 25, 26, 27/)
 idx(8,:)  = (/27, 28, 29, 30, 31/)
 idx(9,:)  = (/31, 32, 33, 34, 35/)
 idx(10,:) = (/35, 36, 37, 38,  0/)
 idx(11,:) = (/39, 40, 41, 42,  0/)
 idx(12,:) = (/42, 43, 44, 45, 46/)

 ! weights of 8-day's data
 wgt(1,:)  = (/8, 8, 8, 7, 0/)
 wgt(2,:)  = (/1, 8, 8, 8, 3/)
 wgt(3,:)  = (/5, 8, 8, 8, 2/)
 wgt(4,:)  = (/6, 8, 8, 8, 0/)
 wgt(5,:)  = (/8, 8, 8, 7, 0/)
 wgt(6,:)  = (/1, 8, 8, 8, 5/)
 wgt(7,:)  = (/3, 8, 8, 8, 4/)
 wgt(8,:)  = (/4, 8, 8, 8, 3/)
 wgt(9,:)  = (/5, 8, 8, 8, 1/)
 wgt(10,:) = (/7, 8, 8, 8, 0/)
 wgt(11,:) = (/8, 8, 8, 6, 0/)
 wgt(12,:) = (/2, 8, 8, 8, 5/)

 DO i=1,46,1
    days(i) = i
 ENDDO

 DO i=1,12,1
    mons(i) = i
 ENDDO

 DO i=1,16,1
    pftnum(i) = i-1
 ENDDO

 dom = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

 !PFT names
 pftname = (/"not_vegetated                           ", &
             "needleleaf_evergreen_temperate_tree     ", &
             "needleleaf_evergreen_boreal_tree        ", &
             "needleleaf_deciduous_boreal_tree        ", &
             "broadleaf_evergreen_tropical_tree       ", &
             "broadleaf_evergreen_temperate_tree      ", &
             "broadleaf_deciduous_tropical_tree       ", &
             "broadleaf_deciduous_temperate_tree      ", &
             "broadleaf_deciduous_boreal_tree         ", &
             "broadleaf_evergreen_temperate_shrub     ", &
             "broadleaf_deciduous_temperate_shrub     ", &
             "broadleaf_deciduous_boreal_shrub        ", &
             "c3_arctic_grass                         ", &
             "c3_non-arctic_grass                     ", &
             "c4_grass                                ", &
             "c3_crop                                 "/)

 ! values in the table below are from Lawrence et al., 2007
 ! with modifications according to Sitch et al., 2003
 laimax = (/0, 5, 5, 5, 7, 7, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4/)
 tbcase = (/0, 0, 0, 2, 0, 0, 0, 5, 5, 0, 5, 5, 5, 5, 5, 5/)
 sgdd   = (/0, 0, 0, 100, 0, 0, 0, 200, 200, 0, 100, 100, 100, 100, 100, 100/)
 minfr  = (/0., 0.7, 0.7, 0., 0.8, 0.8, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
 saimin = (/0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.1/)
 sairtn = (/0., 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, &
            0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0./)

 phimin(:)   = 0.5
 phi   (:,:) = 1.
 laiini(:,:) = 0.

 OPEN(11,FILE=REGFILE)
 OPEN(12,FILE=REGFILE)

 DO WHILE(.TRUE.)

    READ(11,*,END=100) reg
    READ(12,*,END=101) creg
    filename = TRIM(ROOT_DIR)//TRIM(RAW_DIR)//'RG_'//TRIM(creg(1))//'_'//&
               TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.RAW'//TRIM(year)//'.nc'

    PRINT*, filename

    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )

    ! get raw data
    CALL check( nf90_inq_varid(ncid, 'LC'  , lcid    ) )
    CALL check( nf90_get_var  (ncid, lcid  , lcdata  ) )

    CALL check( nf90_inq_varid(ncid, 'PCTT', pcttid  ) )
    CALL check( nf90_get_var  (ncid, pcttid, pcttdata) )

    CALL check( nf90_inq_varid(ncid, 'PCTH', pcthid  ) )
    CALL check( nf90_get_var  (ncid, pcthid, pcthdata) )

    CALL check( nf90_inq_varid(ncid, 'PCTB', pctbid  ) )
    CALL check( nf90_get_var  (ncid, pctbid, pctbdata) )

    CALL check( nf90_inq_varid(ncid, 'BT'  , btid    ) )
    CALL check( nf90_get_var  (ncid, btid  , btdata  ) )

    CALL check( nf90_inq_varid(ncid, 'NT'  , ntid    ) )
    CALL check( nf90_get_var  (ncid, ntid  , ntdata  ) )

    CALL check( nf90_inq_varid(ncid, 'ET'  , etid    ) )
    CALL check( nf90_get_var  (ncid, etid  , etdata  ) )

    CALL check( nf90_inq_varid(ncid, 'DT'  , dtid    ) )
    CALL check( nf90_get_var  (ncid, dtid  , dtdata  ) )

    CALL check( nf90_inq_varid(ncid, 'KG'  , kgid    ) )
    CALL check( nf90_get_var  (ncid, kgid  , kgdata  ) )

    CALL check( nf90_inq_varid(ncid, 'LAI' , laiid   ) )
    CALL check( nf90_get_var  (ncid, laiid , laidata ) )

    CALL check( nf90_inq_varid(ncid, 'PREC', precid  ) )
    CALL check( nf90_get_var  (ncid, precid, precdata) )

    CALL check( nf90_inq_varid(ncid, 'TAVG', tavgid  ) )
    CALL check( nf90_get_var  (ncid, tavgid, tavgdata) )

    CALL check( nf90_inq_varid(ncid, 'TMAX', tmaxid  ) )
    CALL check( nf90_get_var  (ncid, tmaxid, tmaxdata) )

    CALL check( nf90_inq_varid(ncid, 'TMIN', tminid  ) )
    CALL check( nf90_get_var  (ncid, tminid, tmindata) )

    CALL check( nf90_inq_varid(ncid, 'HTOP', htopid  ) )
    CALL check( nf90_get_var  (ncid, htopid, htopdata) )

    CALL check( nf90_close(ncid) )

    PRINT*, 'raw data read completed'

    filename = TRIM(ROOT_DIR)//TRIM(CCI_DIR)//'RG_'//TRIM(creg(1))//'_'//&
               TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.CCI'//TRIM(year)//'.nc'

    PRINT*, filename

    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )

    CALL check( nf90_inq_varid(ncid, 'Tree_Broadleaf_Evergreen'  , tbeid  ) )
    CALL check( nf90_get_var  (ncid, tbeid                       , tbedata) )

    CALL check( nf90_inq_varid(ncid, 'Tree_Broadleaf_Deciduous'  , tbdid  ) )
    CALL check( nf90_get_var  (ncid, tbdid                       , tbddata) )

    CALL check( nf90_inq_varid(ncid, 'Tree_Needleleaf_Evergreen' , tneid  ) )
    CALL check( nf90_get_var  (ncid, tneid                       , tnedata) )

    CALL check( nf90_inq_varid(ncid, 'Tree_Needleleaf_Deciduous' , tndid  ) )
    CALL check( nf90_get_var  (ncid, tndid                       , tnddata) )

    CALL check( nf90_inq_varid(ncid, 'Shrub_Broadleaf_Evergreen' , sbeid  ) )
    CALL check( nf90_get_var  (ncid, sbeid                       , sbedata) )

    CALL check( nf90_inq_varid(ncid, 'Shrub_Broadleaf_Deciduous' , sbdid  ) )
    CALL check( nf90_get_var  (ncid, sbdid                       , sbddata) )

    CALL check( nf90_inq_varid(ncid, 'Shrub_Needleleaf_Evergreen', sneid  ) )
    CALL check( nf90_get_var  (ncid, sneid                       , snedata) )

    CALL check( nf90_inq_varid(ncid, 'Natural_Grass'             , ngid   ) )
    CALL check( nf90_get_var  (ncid, ngid                        , ngdata ) )

    CALL check( nf90_inq_varid(ncid, 'Managed_Grass'             , mgid   ) )
    CALL check( nf90_get_var  (ncid, mgid                        , mgdata ) )

    CALL check( nf90_inq_varid(ncid, 'Urban_areas'               , uaid   ) )
    CALL check( nf90_get_var  (ncid, uaid                        , uadata ) )

    CALL check( nf90_close(ncid) )

    PRINT*, 'CCI data read completed'
    !!!!!!!!!!
    !149-158
    !!!!!!!!!!
    ! convert nan value to 0
    WHERE(tbedata /= tbedata) tbedata=0.
    WHERE(tbddata /= tbddata) tbddata=0.
    WHERE(tnedata /= tnedata) tnedata=0.
    WHERE(tnddata /= tnddata) tnddata=0.
    WHERE(sbedata /= sbedata) sbedata=0.
    WHERE(sbddata /= sbddata) sbddata=0.
    WHERE(snedata /= snedata) snedata=0.
    WHERE(ngdata  /= ngdata)  ngdata =0.
    WHERE(mgdata  /= mgdata)  mgdata =0.
    WHERE(uadata  /= uadata)  uadata =0.

    ! initialization
    pcrop   (:,:) = 0
    purban  (:,:) = 0
    pwetland(:,:) = 0
    pice    (:,:) = 0
    pwater  (:,:) = 0
    pocean  (:,:) = 0

    ppft (:,:,:)  = 0
    lclai(:,:,:)  = 0
    lcsai(:,:,:)  = 0

    pftlai(:,:,:,:) = 0
    pftsai(:,:,:,:) = 0

    !PRINT*, '\n'
    PRINT*, 'Start to process region: ', reg

    ! loop for each small 500m grid
    DO i=1,xydim,1
       !PRINT*, i
       DO j=1,xydim,1
          !PRINT*, j
          ! get data
          !------------

          ! NOTE: 1km, 600x600
          j1 = CEILING(j*1./2)
          i1 = CEILING(i*1./2)

          ! get land cover and kg zone
          lc = lcdata(j,i)
          kg = kgdata(j1,i1)

      ! yuan, 1/1/2020: do not calculate ocean points
      ! not land, supposed to be ocean
          IF (kg == 0) THEN
             lcdata(j,i) = 0
             pocean(j,i) = 100
             CYCLE
          ENDIF

      ! yuan, 1/2/2020: set NA data to ocean
      ! inconsistant with soil data
      ! right now only process RG_25_150_25_155
      ! need to deal with other regions
      ! !!!NOTE: 需要重新run这个程序，使与land cover type保持一致
          ! set NA data to barren
          IF (lc /= lc) THEN
             !lc = 16
             !lcdata(j,i) = 16
             lcdata(j,i) = 0
             pocean(j,i) = 100.
             CYCLE
          ENDIF

          IF (lc == 11) THEN
             pwetland(j,i) = 100.
          ELSE IF (lc == 12) THEN
             pcrop(j,i) = 100.
          ELSE IF (lc == 13) THEN
             purban(j,i) = 100.
          ELSE IF (lc == 14) THEN
             pcrop(j,i) = 50.
          ELSE IF (lc == 15) THEN
             pice (j,i) = 100.
             CYCLE
          ELSE IF (lc == 17) THEN
             pwater(j,i) = 100.
             CYCLE
          ENDIF

          ! get additional data
          pctt   = pcttdata(j,i)
          pcth   = pcthdata(j,i)
          pctb   = pctbdata(j,i)
          lai(:) = laidata (j,i,:)*0.1


          ! NOTE: 1km, 600x600
          bt   = btdata(j1,i1)
          nt   = ntdata(j1,i1)
          et   = etdata(j1,i1)
          dt   = dtdata(j1,i1)

          prec(:) = precdata(j1,i1,:)
          tavg(:) = tavgdata(j1,i1,:)
          tmax(:) = tmaxdata(j1,i1,:)
          tmin(:) = tmindata(j1,i1,:)
          htop500(j,i) = htopdata(j1,i1)

          ! CCI data
          tbe  = tbedata(j,i)*100
          tbd  = tbddata(j,i)*100
          tne  = tnedata(j,i)*100
          tnd  = tnddata(j,i)*100
          sbe  = sbedata(j,i)*100
          sbd  = sbddata(j,i)*100
          sne  = snedata(j,i)*100
          ng   = ngdata (j,i)*100
          mg   = mgdata (j,i)*100
          ua   = uadata (j,i)*100

          ! ------------------------------------------------------------
          ! set crop/urban/water/glacier/wetland(CoLM) pecent
          ! basic steps/rules:
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! step 1: adjust due to water%, total soil reduced  -> (1-water%)
          !
          ! step 2: adjust due to ice% (additonal), and crop%
          !   ice% from bare%
          !   crop% from grass%
          !
          ! step 3: adjust due to urban% and wetland%
          !   urban and wetland have the same vegeta comp as natrural ones
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !
          ! Further update:
          ! the values above need to be updated by higher revolution
          ! data. 用更专业化的数据或更高分辨率的数据对水体、城市、湿地、冰川
          ! 和耕地进行更新。
          ! 注: 目前并没有以上更新，因此水体、城市、湿地和冰川根据 MODIS相应
          ! 的地表分类设置为100%。耕地则进行简单判断，与该格点草本覆盖率相关。
          ! ------------------------------------------------------------

          ! set vegetation fraction
          ! ------------------------------

          ! calcualte tree and herb fractional cover
          tree_cover = 0
          herb_cover = 0

          summ = pctt+pcth+pctb

          IF ( summ/=summ .or. summ==0) THEN  ! check NAN value of fortran
             tree_cover = tbe + tbd + tne + tnd + ua*0.05
             herb_cover = sbe + sbd + sne + ng  + mg + ua*0.15
          ELSE
             tree_cover = pctt/0.8
             herb_cover = pcth/0.8

             summ = tree_cover+herb_cover
             IF (summ > 100) THEN
                tree_cover = tree_cover * 100./summ
                herb_cover = herb_cover * 100./summ
             ENDIF
          ENDIF

          IF (tree_cover+herb_cover == 0) THEN
             ppft (j,i,1) = 100.
             pcrop(j,i)   = 0.
             CYCLE
          ENDIF

          ! bare soil fractional cover
          ppft(j,i,1) = max(0., (100-tree_cover-herb_cover))

          ! split tree
          IF ( tree_cover > 0.) THEN   ! check for end
             ! tropical
             IF (kg<=4 .or. kg==6) THEN
                ! assumed no ndt here
                bt = bt + nt
                nt = 0.

                IF (et+dt > 0) THEN
                   bet = tree_cover*et/(et+dt)
                   bdt = tree_cover*dt/(et+dt)
                ELSE
                   bet = tree_cover*0.5
                   bdt = tree_cover*0.5
                ENDIF

                ! dominated by evergreen broadleaf
                ! Min fradctional 80% of tree cover
                ! type, the same below
                IF (lc == 2) THEN
                   bet = max(tree_cover*0.8, bet)
                   bdt = tree_cover-bet
                ENDIF

                ! dominated by bdt
                IF (lc == 4) THEN
                   bdt = max(tree_cover*0.8, bdt)
                   bet = tree_cover-bdt
                ENDIF

                ppft(j,i,4+1) = bet
                ppft(j,i,6+1) = bdt
             !ENDIF

             ELSE IF (kg <= 16) THEN
             ! temperature
             ! assume no ndt
                IF (dt<=bt .and. nt<=et) THEN
                   bdt = dt
                   net = nt
                   bet = bt-dt
                ELSE
                   !PRINT*, 'Deceduous needleleaf in temperate, remove it!')
                   bdt = bt + (dt-bt)/2.
                   net = et + (nt-et)/2.
                   bet = 0.
                ENDIF

                !PRINT*, bdt, net, bet
                summ = bdt + net + bet
                !summ = sum(bdt + net + bet)
                IF (summ > 0.) THEN
                   bdt = tree_cover * bdt/summ
                   net = tree_cover * net/summ
                   bet = tree_cover * bet/summ
                ELSE
                   bdt = tree_cover * 1/2.
                   net = tree_cover * 1/2.
                   bet = 0.
                ENDIF

                ! dominated by net
                IF (lc == 1) THEN
                   net = max(tree_cover*0.8, net)
                   summ= bdt + bet

                   IF (summ > 0.) THEN
                      bdt = (tree_cover-net)*bdt/summ
                      bet = (tree_cover-net)*bet/summ
                   ENDIF
                ENDIF

                ! dominated by bet
                IF (lc == 2) THEN
                   bet = max(tree_cover*0.8, bet)
                   summ= bdt + net

                   IF (summ > 0.) THEN
                      bdt = (tree_cover-bet)*bdt/summ
                      net = (tree_cover-bet)*net/summ
                   ENDIF
                ENDIF

                ! dpminated by bdt
                IF (lc == 4) THEN
                   bdt = max(tree_cover*0.8, bdt)
                   summ= bet+net

                   IF (summ > 0.) THEN
                      bet = (tree_cover-bdt)*bet/summ
                      net = (tree_cover-bdt)*net/summ
                   ENDIF
                ENDIF

                ppft(j,i,1+1) = net
                ppft(j,i,5+1) = bet
                ppft(j,i,7+1) = bdt

             ELSE IF (kg>=17 .and. kg<=30) THEN
                ! boreal and polar
                ! assume no bet
                IF (bt<=dt .and. et<=nt) THEN
                   bdt = bt
                   net = et
                   ndt = dt - bt
                ELSE
                   !PRINT*, 'Evergreen broadleaf in bereal, remove it!'
                   bdt = dt + (bt-dt)/2.
                   net = nt + (et-nt)/2.
                   ndt = 0.
                ENDIF

                summ = bdt + net + ndt
                IF (summ > 0.) THEN
                   bdt = tree_cover * bdt/summ
                   net = tree_cover * net/summ
                   ndt = tree_cover * ndt/summ
                ELSE
                   bdt = tree_cover * 1/2.
                   net = tree_cover * 1/2.
                   ndt = 0.
                ENDIF

                ! dominated by net
                ! min fractional 80% of tree cover
                IF (lc == 1) THEN
                   net = max(tree_cover*0.8, net)
                   summ= bdt + ndt

                   IF (summ > 0.) THEN
                      bdt = (tree_cover-net)*bdt/summ
                      ndt = (tree_cover-net)*ndt/summ
                   ENDIF
                ENDIF

                ! dominated by ndt
                IF (lc == 3) THEN
                   ndt = max(tree_cover*0.8, ndt)
                   summ= bdt + net

                   IF (summ > 0.) THEN
                      bdt = (tree_cover-ndt)*bdt/summ
                      net = (tree_cover-ndt)*net/summ
                   ENDIF
                ENDIF

                ! dominated by bdt
                IF (lc == 4) THEN
                   bdt = max(tree_cover*0.8, bdt)
                   summ= ndt + net

                   IF (summ > 0.) THEN
                      ndt = (tree_cover-bdt)*ndt/summ
                      net = (tree_cover-bdt)*net/summ
                   ENDIF
                ENDIF

                ppft(j,i,2+1) = net
                ppft(j,i,3+1) = ndt
                ppft(j,i,8+1) = bdt
             ENDIF
          ELSE
             ppft(j,i,2:9) = 0.
          ENDIF

          ! split herbaceous using cci data
          IF (herb_cover > 0.) THEN
             IF (kg < 17) THEN
                ! tropical or temperature
                sbe = sbe + sne
                ng  = ng  + mg
                summ= sbe + sbd + ng

                IF (summ > 0.) THEN
                   sbe = herb_cover * sbe/summ
                   sbd = herb_cover * sbd/summ
                   ng  = herb_cover * ng /summ
                ELSE
                   sbe = herb_cover * 1/4.
                   sbd = herb_cover * 1/4.
                   ng  = herb_cover * 1/2.
                ENDIF

                ppft(j,i, 9+1) = sbe
                ppft(j,i,10+1) = sbd
                ppft(j,i,13+1) = ng

             ELSE IF (kg < 29) THEN
                ! boreal
                sbd = sbd + sbe + sne
                ng  = ng  + mg
                summ= sbd + ng

                IF (summ > 0.) THEN
                   sbd = herb_cover * sbd/summ
                   ng  = herb_cover * ng /summ
                ELSE
                   sbd = herb_cover * 1/2.
                   ng  = herb_cover * 1/2.
                ENDIF

                ppft(j,i,11+1) = sbd
                ppft(j,i,13+1) = ng
             ELSE
                ! polar
                ppft(j,i,12+1) = herb_cover
             ENDIF
          ELSE
             ppft(j,i,10:16) = 0.
          ENDIF

          pctpft(:) = ppft(j,i,:)
          ! calculate LAI
          ! ------------------------

          ! covert 8-dau lai to monthly lai
          ! loop for each month
          laitot = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12./)  !not define

          DO imonth=1,12,1
             IF (imonth==1 .or. imonth==4 .or. imonth==5 .or. imonth==10 .or. imonth==11) THEN
                laitot(imonth) = sum(lai(idx(imonth,1:4))*wgt(imonth,1:4))/ &
                                 sum(wgt(imonth,1:4))  !need check for nan
             ELSE
                laitot(imonth) = sum(lai(idx(imonth,:))*wgt(imonth,:))/ &
                                 sum(wgt(imonth,:))  !need check for nan
             ENDIF
          ENDDO

          ! calculate GDD and phi
          phi(:,:) = 1.

          IF (all(tavg==-3.4e+38)) THEN
             ! when T and P are missing
             PRINT*, 'NA temperature! Error!'
          ELSE
             loc1    = minloc(tavg) ! check difference between min and which.min(R)
             itmin   = loc1(1)
             mcnt    = 0
             dd2  (:)= 0.
             dd5  (:)= 0.
             gdd2 (:)= 0.
             gdd5 (:)= 0.
             rgdd2(:)= 0.
             rgdd5(:)= 0.
          ENDIF

          DO WHILE (mcnt < 12)

             imonth = mod((itmin+mcnt-1),12) + 1

             ! calculate gdd according to tmin, tmax, tavg, 2, 5degree
             ! ------------------------------

             ! calcualte day index of tavg
             IF ((tmax(imonth)-tmin(imonth)) > 0.) THEN
                davg = (tmax(imonth)-tavg(imonth))/ &
                      (tmax(imonth)-tmin(imonth))*dom(imonth)
             ELSE
                davg = 0
             ENDIF

             ! gdd2 in different cases
             IF (2. < tmin(imonth)) THEN
                dd2(imonth) = (tavg(imonth)-2.)* dom(imonth)
             ELSE IF (2. >= tmax(imonth)) THEN
                dd2(imonth) = 0.
             ELSE IF (2.>=tmin(imonth) .and. 2.<tavg(imonth) .and. tavg(imonth)/=tmin(imonth)) THEN
                dd2(imonth) = (tavg(imonth)-2.)/2.*(tavg(imonth)-2.)/ &
                            (tavg(imonth)-tmin(imonth))*davg + &
                            (tmax(imonth)+tavg(imonth)-4.)/2.* &
                            (dom(imonth)-davg)
             ELSE IF (2.>=tavg(imonth) .and. 2.<tmax(imonth) &
                   .and. tmax(imonth)/=tavg(imonth)) THEN
                dd2(imonth) = (tmax(imonth)-2.)/2.*(tmax(imonth)-2.)/ &
                            (tmax(imonth)-tavg(imonth))*(dom(imonth)-davg)
             ELSE
                dd2(imonth) = max(0., ((tavg(imonth)-2.)*dom(imonth)))
             ENDIF

             ! gdd5 in different cases
             IF (5. < tmin(imonth)) THEN
                dd5(imonth) = (tavg(imonth)-5.)* dom(imonth)
             ELSE IF (5. >= tmax(imonth)) THEN
                dd5(imonth) = 0.
             ELSE IF (5.>=tmin(imonth) .and. 5.<tavg(imonth) &
                   .and. tavg(imonth)/=tmin(imonth)) THEN
                dd5(imonth) = (tavg(imonth)-5.)/2.*(tavg(imonth)-5.)/ &
                            (tavg(imonth)-tmin(imonth))*davg + &
                            (tmax(imonth)+tavg(imonth)-10.)/2.* &
                            (dom(imonth)-davg)
             ELSE IF (5.>=tavg(imonth) .and. 5.<tmax(imonth) &
                   .and. tmax(imonth)/=tavg(imonth)) THEN
                dd5(imonth) = (tmax(imonth)-5.)/2.*(tmax(imonth)-5.)/ &
                            (tmax(imonth)-tavg(imonth))*(dom(imonth)-davg)
             ELSE
                dd5(imonth) = max(0., ((tavg(imonth)-5.)*dom(imonth)))
             ENDIF

             mcnt = mcnt + 1
          ENDDO

          ! split the tmin month gdd into 2 parts
          ! gdd, rgdd

          IF (dd2(itmin) == 0) THEN
             gdd2 (itmin) = 0.
             rgdd2(itmin) = 0.
          ELSE
             nextmonth = mod(itmin,12) + 1
             prevmonth = mod((10+itmin),12) + 1
             summ      = max(0., (tmin(nextmonth)-2.)) + &
                         max(0., (tmin(prevmonth)-2.))

             IF (summ > 0.) THEN
                gdd2 (itmin) = dd2(itmin)*max(0., (tmin(nextmonth)-2.))/summ
                rgdd2(itmin) = dd2(itmin)*max(0., (tmin(prevmonth)-2.))/summ
             ELSE
                IF ((tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin)) > 0) THEN
                   gdd2 (itmin) = dd2(itmin)*(tavg(nextmonth)-tavg(itmin))/ &
                                  (tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin))
                   rgdd2(itmin) = dd2(itmin)*(tavg(prevmonth)-tavg(itmin))/ &
                                  (tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin))
                ELSE
                   gdd2 (itmin) = dd2(itmin) * 0.5
                   rgdd2(itmin) = dd2(itmin) * 0.5
                ENDIF
             ENDIF
          ENDIF

          IF (dd5(itmin) == 0) THEN
             gdd5 (itmin) = 0.
             rgdd5(itmin) = 0.
          ELSE
             nextmonth = mod(itmin,12) + 1
             prevmonth = mod((itmin+10),12) + 1
             summ      = max(0., (tmin(nextmonth)-5.)) + &
                         max(0., (tmin(prevmonth)-5.))

             IF (summ > 0.) THEN
                gdd5 (itmin) = dd5(itmin)*max(0., tmin(nextmonth)-5.)/summ
                rgdd5(itmin) = dd5(itmin)*max(0., tmin(prevmonth)-5.)/summ
             ELSE
                IF ((tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin)) > 0) THEN
                   gdd5 (itmin) = dd5(itmin)*(tavg(nextmonth)-tavg(itmin))/ &
                                  (tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin))
                   rgdd5(itmin) = dd5(itmin)*(tavg(prevmonth)-tavg(itmin))/ &
                                  (tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin))
                ELSE
                   gdd5 (itmin) = dd5(itmin) * 0.5
                   rgdd5(itmin) = dd5(itmin) * 0.5
                ENDIF
             ENDIF
          ENDIF

          imonth_prev = itmin
          mcnt        = 1
          DO WHILE (mcnt < 12)

             imonth = mod((itmin+mcnt-1),12) + 1

             gdd2(imonth) = gdd2(imonth_prev) + dd2(imonth)
             gdd5(imonth) = gdd5(imonth_prev) + dd5(imonth)

             imonth_prev  = imonth
             mcnt         = mcnt + 1
          ENDDO

          imonth_prev = itmin
          mcnt        = 1
          DO WHILE (mcnt < 12)

             imonth = mod((12+itmin-mcnt-1),12) + 1

             rgdd2(imonth) = rgdd2(imonth_prev) + dd2(imonth)
             rgdd5(imonth) = rgdd5(imonth_prev) + dd5(imonth)

             imonth_prev = imonth
             mcnt        = mcnt + 1
          ENDDO

          gdd2 (itmin) = dd2(itmin)
          rgdd2(itmin) = dd2(itmin)
          gdd5 (itmin) = dd5(itmin)
          rgdd5(itmin) = dd5(itmin)

          ! calculate phi using gdd/sgdd
          !!!!!!
          DO imonth=1,12,1
             gdd2(imonth)   = min(gdd2(imonth), rgdd2(imonth))
             gdd5(imonth)   = min(gdd5(imonth), rgdd5(imonth))

             phi (4,imonth) = max(phimin(imonth), &
                                  min(1., gdd2(imonth)/sgdd(4)))
             phi (8,imonth) = max(phimin(imonth), &
                                  min(1., gdd5(imonth)/sgdd(8)))
             phi (9,imonth) = max(phimin(imonth), &
                                  min(1., gdd5(imonth)/sgdd(9)))
             phi (11,imonth)= max(phimin(imonth), &
                                  min(1., gdd5(imonth)/sgdd(11)))
             phi (12,imonth)= max(phimin(imonth), &
                                  min(1., gdd5(imonth)/sgdd(12)))
             phi (13,imonth)= max(phimin(imonth), &
                                  min(1., gdd5(imonth)/sgdd(13)))
             phi (14,imonth)= max(phimin(imonth), &
                                  min(1., gdd5(imonth)/sgdd(14)))
             phi (15,imonth)= max(phimin(imonth), &
                                  min(1., gdd5(imonth)/sgdd(15)))
             phi (16,imonth)= max(phimin(imonth), &
                                  min(1., gdd5(imonth)/sgdd(16)))
          ENDDO

          ! distribution LAI
          ! sum(lai*PFT%) = MODIS_LAI

          ! calculate initial LAI
          laiini(:,:) = 0.
          pctpft(:)   = ppft(j,i,:)/100.

          DO imonth=1,12,1
             sumwgt = sum(phi(:,imonth)*laimax(:)*pctpft(:))

             IF (sumwgt > 0.) THEN
                laiini(:,imonth) = phi(:,imonth)*laimax(:)/sumwgt*laitot(imonth)
             ELSE
                laiini(:,imonth) = 0.
             ENDIF
          ENDDO

          ! adjust LAI
          ! -----------------------------
          !laiinimax = apply(laiini, 1, max) !check
          loc2      = maxloc(laitot)
          laimaxloc = loc2(1)
          laiinimax(:) = laiini(:,laimaxloc)

          DO imonth=1,12,1

             ! evergreen tree phenology (Zeng et al., 2002)
             ! minimum fraction of maximum initial PFT LAI
             ! laiini[,imonth] = pmax(laiini[,imonth], laiinimax*minfr)

             ! max value limited
             DO iloop=1,16
                IF (laiini(iloop,imonth) > laiinimax(iloop)) THEN
                   laiini(iloop,imonth) = laiinimax(iloop)
                ENDIF
             ENDDO

             ! calculate for nonevergreen PFT
             !indx1  = (/2,3,5,6,10/)
             sumevg = sum(laiini(indx1(:),imonth)*pctpft(indx1(:)))
             laitot_nonevg = max(laitot(imonth)-sumevg, 0.)

             !indx2  = (/4,7,8,9,11,12,13,14,15,16/)

             sumnon = sum(phi(indx2(:),imonth)*laimax(indx2(:))*pctpft(indx2(:)))
             IF (sumnon > 0.) THEN
                DO inx=1,10
                   laiini(indx2(inx),imonth) = phi(indx2(inx),imonth) * &
                                       laimax(indx2(inx)) / sumnon * laitot_nonevg
                ENDDO
             ELSE
                sumnon = sum(laimax(indx2(:))*pctpft(indx2(:)))
                IF (sumnon > 0.) THEN
                   ! do not consider phenology
                   DO inx=1,10
                      laiini(indx2(inx),imonth) = laimax(indx2(inx)) / sumnon * laitot_nonevg
                   ENDDO
                ELSE
                   ! no percentage cover
                   DO inx=1,10
                      laiini(indx2(inx),imonth) = 0._r8
                   ENDDO
                ENDIF
             ENDIF
          ENDDO

          ! max LAI value limited (need more test and check)
          IF (count(laiini>10.)>0) THEN
             WHERE (laiini>10.) laiini=10.
          ENDIF

          ! split C3/C4 nature grass (Still et al., 2003)
          ! ------------------------------

          ! C4 case
          IF (count(tavg>22.)==12 .and. count(prec>25.)==12) THEN
             ppft(j,i,14+1)  = sum(ppft(j,i,13:14)) ! index check
             ppft(j,i,13:14) = 0.                   ! index check
          ELSE IF (herb_cover>75 .and. count(laiini(14,:)>0)>0) THEN
             ! minxed C3/C4 case
             laiup = 0.
             DO inx=1,12
                IF (tavg(inx)>22. .and. prec(inx)>25.) THEN
                   laiup = laiup + laiini(13+1,inx)
                ENDIF
             ENDDO
             frac_c4 = laiup/sum(laiini(14,:))

             ppft(j,i,14+1) = ppft(j,i,13+1) * frac_c4
             ppft(j,i,13+1) = ppft(j,i,13+1) * (1.-frac_c4)
          ENDIF

          ! calculate SAI
          ! --------------------
          DO ipft=1,pfts,1
             saimin1(ipft) = saimin(ipft)*maxval(laiini(ipft,:))
          ENDDO
          saimin1(2:16)= saimin1(2:16)/laimax(2:16)
          saiini (:,:) = 0.
          saiini1(:,:) = saiini(:,:)


          DO iloop=1,12
             saiini1(:,iloop) = saimin(:)
          ENDDO

          ! PFT LAI diff
          !laidiff(:,:) = laiini(:,1:12) - laiini(:,:)
          laidiff(:,1) = laiini(:,12) - laiini(:,1)
          DO iloop=1,11,1
             inx = iloop+1
             laidiff(:,inx) = laiini(:,iloop) - laiini(:,inx)
          ENDDO

          WHERE(laidiff<0.) laidiff=0.

          iloop = 1
          DO WHILE (iloop<13)
             saires(:) = sairtn*saiini1(:,mod((iloop+10),12)+1)
             DO inx=1,16
                x1 = (saires(inx)+laidiff(inx,iloop)*0.5)
                x2 = saimin1(inx)
                saiini(inx,iloop) = max(x1,x2) !maxval((saires(:)+laidiff(:,imonth)*0.5), saimin1())
             ENDDO
             iloop = iloop + 1

             IF (iloop == 13) THEN
                sum_judg = abs(sum(saiini-saiini1))
                IF (sum_judg > 1e-6) THEN
                   iloop  = 1
                   saiini1(:,:) = saiini(:,:)
                ENDIF
             ENDIF
          ENDDO

          ! adjust herb to cropland
          ! ---------------------------------

          ! crop exist
          IF (pcrop(j,i) > 0.) THEN
             ! if herb frac < crop, set it all to crop
             IF (herb_cover <= pcrop(j,i)) THEN
                ppft(j,i,16) = herb_cover
                ppft(j,i,10:15) = 0.
                laiini(10:15,:) = 0.
                saiini(10:15,:) = 0.
                pcrop(j,i) = herb_cover
             ELSE
                ! if herb > crop, remove the crop
                ! fraction fro herb
                ppft(j,i,16) = pcrop(j,i)
                ppft(j,i,10:15) = ppft(j,i,10:15)*(herb_cover-pcrop(j,i))/ &
                                  herb_cover
             ENDIF
          ENDIF

          ! use adjusted PCT for crop and C3/C4
          pctpft(:) = ppft(j,i,:)/100.
          DO ipft=1,16,1
             IF (pctpft(ipft) < 1e-6) THEN
                laiini(ipft,:) = 0.
                saiini(ipft,:) = 0.
             ENDIF
          ENDDO

          IF (count(laiini>20) > 0) THEN
             PRINT*, 'check for laiini!'
             STOP
          ENDIF

          ! max SAI value limited
          IF (count(saiini>3) > 0) THEN
             WHERE(saiini>3.) saiini = 3.
          ENDIF

          IF (abs(sum(pctpft)-1.) > 1e-5) THEN
             PRINT*, pctpft
             PRINT*, abs(sum(pctpft)-1.)
             PRINT*, 'Sum of area is not equle to 100%! STOP!'
             STOP
          ENDIF

          IF (abs(pcrop(j,i)-ppft(j,i,16)) > 1e-3) THEN
             PRINT*, 'Crop area is not conserved! STOP!'
             STOP
          ENDIF

          ! output data
          pftlai(j,i,:,:) = laiini(:,:)
          pftsai(j,i,:,:) = saiini(:,:)

          DO iloop=1,12
             laiini1(:,iloop) = pctpft(:)*laiini(:,iloop)
             saiini2(:,iloop) = pctpft(:)*saiini(:,iloop)
          ENDDO

          DO iloop=1,12
             lclai(j,i,iloop) = sum(laiini1(:,iloop))
             lcsai(j,i,iloop) = sum(saiini2(:,iloop))
          ENDDO
          IF (maxval(abs(lclai(j,i,:)-laitot(:))) > 1) THEN
             !PRINT*, 'LAI not conserved, set it to laitot'
             !PRINT*, laitot
             !PRINT*, lclai(j,i,:)

             lclai(j,i,:) = laitot(:)
             IF (abs(sum(pctpft)-1.) > 1e-5) THEN
                PRINT*, 'check'
                STOP
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    dll = (reg(4)-reg(2))*1./(xydim*1.)

    DO dims=1,xydim
       lons(dims) = reg(2) + dims*dll - dll/2
       lats(dims) = reg(1) - dims*dll + dll/2
    ENDDO

    filename = TRIM(OUT_DIR)//TRIM(SRF_DIR)//'RG_'//TRIM(creg(1))//'_'//&
               TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))//'.MOD'//TRIM(year)//'.nc'
    PRINT*, filename
    CALL check( nf90_create(filename, NF90_NETCDF4, ncid) )

    ! define dimensions
    CALL check( nf90_def_dim(ncid, "lat", xydim, lat_dimid) )
    CALL check( nf90_def_dim(ncid, "lon", xydim, lon_dimid) )
    CALL check( nf90_def_dim(ncid, "pft", pfts , pft_dimid) )
    CALL check( nf90_def_dim(ncid, "mon", mon  , mon_dimid) )

    ! define variables
    ! ---------------------------------------
    CALL check( nf90_def_var(ncid, "lat", NF90_FLOAT, lat_dimid, lat_vid, deflate_level=6) )
    CALL check( nf90_def_var(ncid, "lon", NF90_FLOAT, lon_dimid, lon_vid, deflate_level=6) )
    CALL check( nf90_def_var(ncid, "pft", NF90_INT  , pft_dimid, pft_vid, deflate_level=6) )
    CALL check( nf90_def_var(ncid, "mon", NF90_INT  , mon_dimid, mon_vid, deflate_level=6) )

    CALL check( nf90_put_att(ncid, lat_vid, "long_name", "Latitude"     ) )
    CALL check( nf90_put_att(ncid, lat_vid, "units"    , "degrees_north") )
    CALL check( nf90_put_att(ncid, lon_vid, "long_name", "Longitude"    ) )
    CALL check( nf90_put_att(ncid, lon_vid, "units"    , "degrees_east" ) )
    CALL check( nf90_put_att(ncid, pft_vid, "long_name", "Index of PFT" ) )
    CALL check( nf90_put_att(ncid, pft_vid, "units"    , "-"            ) )
    CALL check( nf90_put_att(ncid, mon_vid, "long_name", "Month of year") )
    CALL checK( nf90_put_att(ncid, mon_vid, "units"    , "month"        ) )

    CALL check( nf90_put_var(ncid, lat_vid, lats  ) )
    CALL check( nf90_put_var(ncid, lon_vid, lons  ) )
    CALL check( nf90_put_var(ncid, pft_vid, pftnum) )
    CALL check( nf90_put_var(ncid, mon_vid, mons  ) )
    ! land cover data

    shortvalue = 255
    XY2D = (/lon_dimid, lat_dimid/)
    XY3D = (/lon_dimid, lat_dimid, mon_dimid/)
    XY3F = (/lon_dimid, lat_dimid, pft_dimid/)
    XY4D = (/lon_dimid, lat_dimid, pft_dimid, mon_dimid/)

    CALL check( nf90_def_var(ncid, "LC" , NF90_SHORT  , XY2D, lc_id, deflate_level=6) )
    CALL checK( nf90_put_att(ncid, lc_id, "units"     , "-"        ) )
    CALL check( nf90_put_att(ncid, lc_id, "long_name" , "MODIS Land Cover Type (LC_Type1) data product, MCD12Q1 V006") )
    CALL check( nf90_put_att(ncid, lc_id, "_FillValue", shortvalue ) )

    fillvalue = -999.
    CALL check( nf90_def_var(ncid, "MONTHLY_LC_LAI", NF90_FLOAT  , XY3D, clai_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, clai_id         , "units"     , "m^2/m^2"    ) )
    CALL check( nf90_put_att(ncid, clai_id         , "long_name" , "Monthly landcover LAI values") )
    CALL check( nf90_put_att(ncid, clai_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "MONTHLY_LC_SAI", NF90_FLOAT  , XY3D, csai_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, csai_id         , "units"     , "m^2/m^2"    ) )
    CALL check( nf90_put_att(ncid, csai_id         , "long_name" , "Monthly landcover SAI values") )
    CALL check( nf90_put_att(ncid, csai_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "MONTHLY_LAI"   , NF90_FLOAT  , XY4D, plai_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, plai_id         , "units"     , "m^2/m^2"    ) )
    CALL check( nf90_put_att(ncid, plai_id         , "long_name" , "Monthly PFT LAI values") )
    CALL check( nf90_put_att(ncid, plai_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "MONTHLY_SAI"   , NF90_FLOAT  , XY4D, psai_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, psai_id         , "units"     , "m^2/m^2"    ) )
    CALL check( nf90_put_att(ncid, psai_id         , "long_name" , "Monthly PFT SAI values") )
    CALL check( nf90_put_att(ncid, psai_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "PCT_CROP"      , NF90_FLOAT  , XY2D, crop_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, crop_id         , "units"     , "%"          ) )
    CALL check( nf90_put_att(ncid, crop_id         , "long_name" , "Percent crop cover") )
    CALL check( nf90_put_att(ncid, crop_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "PCT_URBAN"     , NF90_FLOAT  , XY2D, urbn_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, urbn_id         , "units"     , "%"          ) )
    CALL check( nf90_put_att(ncid, urbn_id         , "long_name" , "Percent urban cover") )
    CALL check( nf90_put_att(ncid, urbn_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "PCT_WETLAND"   , NF90_FLOAT  , XY2D, wetl_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, wetl_id         , "units"     , "%"          ) )
    CALL check( nf90_put_att(ncid, wetl_id         , "long_name" , "Percent wetland cover") )
    CALL check( nf90_put_att(ncid, wetl_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "PCT_GLACIER"   , NF90_FLOAT  , XY2D, gice_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, gice_id         , "units"     , "%"          ) )
    CALL check( nf90_put_att(ncid, gice_id         , "long_name" , "Percent glacier/ice cover") )
    CALL check( nf90_put_att(ncid, gice_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "PCT_WATER"     , NF90_FLOAT  , XY2D, watr_id, deflate_level=6) )
    CALL checK( nf90_put_att(ncid, watr_id         , "units"     , "%"          ) )
    CALL check( nf90_put_att(ncid, watr_id         , "long_name" , "Percent water body cover") )
    CALL check( nf90_put_att(ncid, watr_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "PCT_OCEAN"     , NF90_FLOAT  , XY2D, ocea_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, ocea_id         , "units"     , "%"          ) )
    CALL check( nf90_put_att(ncid, ocea_id         , "long_name" , "Percent ocean cover") )
    CALL check( nf90_put_att(ncid, ocea_id         , "_FillValue", fillvalue    ) )

    CALL check( nf90_def_var(ncid, "PCT_PFT"       , NF90_FLOAT  , XY3F, ppft_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, ppft_id         , "units"     , "%"          ) )
    CALL check( nf90_put_att(ncid, ppft_id         , "long_name" , "Percent PFT cover") )
    CALL checK( nf90_put_att(ncid, ppft_id         , "_FillValue", fillvalue    ) )

    ! tree height
    shortvalue = 255
    CALL check( nf90_def_var(ncid, "HTOP"          , NF90_SHORT  , XY2D, htop_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, htop_id         , "units"     , "m"          ) )
    CALL check( nf90_put_att(ncid, htop_id         , "long_name" , "Global forest canopy height") )
    CALL check( nf90_put_att(ncid, htop_id         , "_FillValue", shortvalue   ) )

    ! put vars
    CALL check( nf90_inq_varid(ncid, "LC"            , lc_id   ) )
    CALL check( nf90_put_var  (ncid, lc_id           , lcdata  ) )
    CALL check( nf90_inq_varid(ncid, "MONTHLY_LC_LAI", clai_id ) )
    CALL check( nf90_put_var  (ncid, clai_id         , lclai   ) )
    CALL check( nf90_inq_varid(ncid, "MONTHLY_LC_SAI", csai_id ) )
    CALL check( nf90_put_var  (ncid, csai_id         , lcsai   ) )
    CALL check( nf90_inq_varid(ncid, "MONTHLY_LAI"   , plai_id ) )
    CALL check( nf90_put_var  (ncid, plai_id         , pftlai  ) )
    CALL check( nf90_inq_varid(ncid, "MONTHLY_SAI"   , psai_id ) )
    CALL check( nf90_put_var  (ncid, psai_id         , pftsai  ) )
    CALL check( nf90_inq_varid(ncid, "PCT_CROP"      , crop_id ) )
    CALL check( nf90_put_var  (ncid, crop_id         , pcrop   ) )
    CALL check( nf90_inq_varid(ncid, "PCT_URBAN"     , urbn_id ) )
    CALL check( nf90_put_var  (ncid, urbn_id         , purban  ) )
    CALL check( nf90_inq_varid(ncid, "PCT_WETLAND"   , wetl_id ) )
    CALL check( nf90_put_var  (ncid, wetl_id         , pwetland) )
    CALL check( nf90_inq_varid(ncid, "PCT_GLACIER"   , gice_id ) )
    CALL check( nf90_put_var  (ncid, gice_id         , pice    ) )
    CALL check( nf90_inq_varid(ncid, "PCT_WATER"     , watr_id ) )
    CALL check( nf90_put_var  (ncid, watr_id         , pwater  ) )
    CALL check( nf90_inq_varid(ncid, "PCT_OCEAN"     , ocea_id ) )
    CALL check( nf90_put_var  (ncid, ocea_id         , pocean  ) )
    CALL check( nf90_inq_varid(ncid, "PCT_PFT"       , ppft_id ) )
    CALL check( nf90_put_var  (ncid, ppft_id         , ppft    ) )
    CALL check( nf90_inq_varid(ncid, "HTOP"          , htop_id ) )
    CALL check( nf90_put_var  (ncid, htop_id         , htop500 ) )

    CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Title'  , Title  ) )
    CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Authors', Authors) )
    CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Adderss', Address) )
    CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Email'  , Email  ) )

    CALL check( nf90_close(ncid) )

 ENDDO

 100 CLOSE(11)
 101 CLOSE(12)

CONTAINS

   SUBROUTINE check(status)
      INTEGER, intent(in) :: status

      IF (status /= nf90_noerr) THEN
         print *, trim( nf90_strerror(status))
         stop 2
      ENDIF
   END SUBROUTINE check
END PROGRAM mkmod
