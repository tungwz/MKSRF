
! ======================================================
! aggreate high-resolution land surface dataset to
! lower resolutioin
!
! History:
!   2019/06: Hua Yuan, initial version
!   2022/03: Hua Yuan, simplify PC aggregation (Mark 03/05/2022)
! ======================================================

PROGRAM MakeSurface

   USE netcdf
   USE omp_lib

   IMPLICIT NONE

   INTEGER, parameter :: r8 = selected_real_kind(12)

   CHARACTER (len=*), parameter :: REGFILE = "reg_5x5"
   CHARACTER (len=*), parameter :: DATASRC = "MOD"
   CHARACTER (len=*), parameter :: SRF_DIR = "/home/yuanhua/tera02/mksrf/srf_5x5/"
   CHARACTER (len=*), parameter :: OUT_DIR = "/home/yuanhua/tera02/mksrf/srf_0.5x0.5/"
   CHARACTER (len=4)            :: year    = "2005"

   CHARACTER (len=*), parameter :: Title   = "Land surface model input vagetation data"
   CHARACTER (len=*), parameter :: Authors = "Yuan et al."
   CHARACTER (len=*), parameter :: Address = "School of Atmospheric Sciences, Sun Yat-sen University, Guangzhou, China"
   CHARACTER (len=*), parameter :: Email   = "yuanh25@mail.sysu.edu.cn"

   INTEGER, parameter :: nxy   = 1200
   INTEGER, parameter :: nxo   = 720
   INTEGER, parameter :: nyo   = 360
   INTEGER, parameter :: nxy1  = 10
   INTEGER, parameter :: nmon  = 12
   INTEGER, parameter :: npft  = 16

   ! define input variables
   REAL(r8), dimension(:,:,:,:), allocatable :: pftlai
   REAL(r8), dimension(:,:,:,:), allocatable :: pftsai
   REAL(r8), dimension(:,:,:)  , allocatable :: lclai
   REAL(r8), dimension(:,:,:)  , allocatable :: lcsai
   REAL(r8), dimension(:,:,:)  , allocatable :: ppft
   REAL(r8), dimension(:,:)    , allocatable :: lcdata
   REAL(r8), dimension(:,:)    , allocatable :: purban
   REAL(r8), dimension(:,:)    , allocatable :: pcrop
   REAL(r8), dimension(:,:)    , allocatable :: pwetland
   REAL(r8), dimension(:,:)    , allocatable :: pglacier
   REAL(r8), dimension(:,:)    , allocatable :: pwater
   REAL(r8), dimension(:,:)    , allocatable :: htop

   ! define output variables
   REAL(r8), dimension(:,:)      , allocatable :: area
   REAL(r8), dimension(:,:)      , allocatable :: pct_land
   REAL(r8), dimension(:,:)      , allocatable :: pct_urban
   REAL(r8), dimension(:,:)      , allocatable :: pct_crop
   REAL(r8), dimension(:,:)      , allocatable :: pct_wetland
   REAL(r8), dimension(:,:)      , allocatable :: pct_glacier
   REAL(r8), dimension(:,:)      , allocatable :: pct_water
   REAL(r8), dimension(:,:,:)    , allocatable :: pct_lc
   REAL(r8), dimension(:,:,:)    , allocatable :: pct_pft
   REAL(r8), dimension(:,:,:)    , allocatable :: htop_lc
   REAL(r8), dimension(:,:,:)    , allocatable :: htop_pft
   REAL(r8), dimension(:,:,:,:)  , allocatable :: pct_pc
   REAL(r8), dimension(:,:,:)    , allocatable :: lai_grid
   REAL(r8), dimension(:,:,:)    , allocatable :: sai_grid
   REAL(r8), dimension(:,:,:,:)  , allocatable :: lai_lc
   REAL(r8), dimension(:,:,:,:)  , allocatable :: sai_lc
   REAL(r8), dimension(:,:,:,:)  , allocatable :: lai_pft
   REAL(r8), dimension(:,:,:,:)  , allocatable :: sai_pft
   REAL(r8), dimension(:,:,:,:,:), allocatable :: lai_pc
   REAL(r8), dimension(:,:,:,:,:), allocatable :: sai_pc

   ! define other variables
   ! -----------------------------------------------

   CHARACTER (len=256) :: FILE_NAME, reg1, reg2, reg3, reg4

   INTEGER :: iostatus, reg(4)
   INTEGER :: ncid, times

   ! variable ids
   INTEGER :: lon_vid, lat_vid, lc_vid, pft_vid, mon_vid
   INTEGER :: lclai_vid, lcsai_vid, lai_vid, sai_vid
   INTEGER :: purban_vid, pcrop_vid, pwater_vid, pglacier_vid, pwetland_vid
   INTEGER :: ppft_vid, htop_vid
   INTEGER :: lat_dimid, lon_dimid, lc_dimid, pft_dimid, mon_dimid

   ! output variables/ids
   INTEGER :: varea, vpct_land, vlai_grid, vsai_grid
   INTEGER :: vlai_lc, vsai_lc, vpct_lc, vhtop_lc
   INTEGER :: vpct_pft, vlai_pft, vsai_pft, vhtop_pft
   INTEGER :: vpct_urban, vpct_crop, vpct_wetland, vpct_glacier, vpct_water
   INTEGER :: vpct_pc, vlai_pc, vsai_pc
   INTEGER :: months(nmon), pfts(npft), lcs(22)

   REAL(r8) :: pi, deg2rad, re, dx, dy, sumarea, wgt, sumpct
   REAL(r8) :: sarea(nxy,nxy)
   REAL(r8) :: dll, lone(nxy), lonw(nxy), latn(nxy), lats(nxy)
   REAL(r8) :: lone1(nxy1), lonw1(nxy1), latn1(nxy1), lats1(nxy1), lons1(nxy1)
   REAL(r8) :: latso(nyo), lonso(nxo)
   REAL(r8) :: glai_lc, glai_pft, glai_pc, gsai_lc, gsai_pft, gsai_pc

   INTEGER :: nlc
   INTEGER :: lc, ip, il, im
   INTEGER :: i, j, k, i1, j1, io, jo
   INTEGER :: URBAN, WETLAND, CROP, WATER, GLACIER

   INTEGER :: XY2D(2), GRID3d(3), LC3D(3), PFT3D(3), LC4D(4), PFT4D(4), PC4D(4), PC5D(5)
   INTEGER :: argn

   pi = 4.*atan(1.)
   deg2rad = pi/180.
   re = 6.37122e6 * 0.001
   times = nxy/nxy1

   IF (DATASRC == "MOD") THEN
      nlc = 17
   ELSE
      nlc = 22
   ENDIF

   ! MODIS IGBP land cover TYPE
   IF (DATASRC == "MOD") THEN
      WETLAND = 11
      URBAN   = 13
      GLACIER = 15
      WATER   = 17
   ENDIF

   ! ESA land cover TYPE
   IF (DATASRC == "ESA") THEN
      URBAN   = 19
      WATER   = 21
      GLACIER = 22
   ENDIF

   ! allocate memory
   allocate( lcdata       (nxy, nxy) )
   allocate( purban       (nxy, nxy) )
   allocate( pcrop        (nxy, nxy) )
   allocate( pwetland     (nxy, nxy) )
   allocate( pglacier     (nxy, nxy) )
   allocate( pwater       (nxy, nxy) )
   allocate( htop         (nxy, nxy) )
   allocate( ppft         (nxy, nxy, npft) )
   allocate( lclai        (nxy, nxy, nmon) )
   allocate( lcsai        (nxy, nxy, nmon) )
   allocate( pftlai       (nxy, nxy, npft, nmon) )
   allocate( pftsai       (nxy, nxy, npft, nmon) )
   allocate( area         (nxo, nyo) )
   allocate( pct_land     (nxo, nyo) )
   allocate( pct_urban    (nxo, nyo) )
   allocate( pct_crop     (nxo, nyo) )
   allocate( pct_wetland  (nxo, nyo) )
   allocate( pct_glacier  (nxo, nyo) )
   allocate( pct_water    (nxo, nyo) )
   allocate( pct_lc       (nxo, nyo, nlc) )
   allocate( htop_lc      (nxo, nyo, nlc) )
   allocate( pct_pft      (nxo, nyo, npft) )
   allocate( htop_pft     (nxo, nyo, npft) )
   allocate( pct_pc       (nxo, nyo, npft, nlc) )
   allocate( lai_grid     (nxo, nyo, nmon) )
   allocate( sai_grid     (nxo, nyo, nmon) )
   allocate( lai_lc       (nxo, nyo, nlc, nmon) )
   allocate( sai_lc       (nxo, nyo, nlc, nmon) )
   allocate( lai_pft      (nxo, nyo, npft, nmon) )
   allocate( sai_pft      (nxo, nyo, npft, nmon) )
   allocate( lai_pc       (nxo, nyo, npft, nlc, nmon) )
   allocate( sai_pc       (nxo, nyo, npft, nlc, nmon) )

   DO i = 1, 22
      lcs(i) = i
   ENDDO

   DO i = 1, npft
      pfts(i) = i - 1
   ENDDO

   DO i = 1, nmon
      months(i) = i
   ENDDO

   ! initialization
   pct_land    (:,:) = 0.
   pct_lc    (:,:,:) = 0.
   lai_grid  (:,:,:) = 0.
   sai_grid  (:,:,:) = 0.
   lai_lc  (:,:,:,:) = 0.
   sai_lc  (:,:,:,:) = 0.
   htop_lc   (:,:,:) = 0.
   lai_pft (:,:,:,:) = 0.
   sai_pft (:,:,:,:) = 0.
   pct_pft   (:,:,:) = 0.
   pct_crop    (:,:) = 0.
   pct_glacier (:,:) = 0.
   pct_water   (:,:) = 0.
   pct_urban   (:,:) = 0.
   pct_wetland (:,:) = 0.
   htop_pft  (:,:,:) = 0.
   lai_pc(:,:,:,:,:) = 0.
   sai_pc(:,:,:,:,:) = 0.
   pct_pc  (:,:,:,:) = 0.
   area        (:,:) = 0.

   ! get args from command line
   argn = IARGC()
   IF (argn > 0) THEN
      CALL getarg(1, year)
   ENDIF

   ! open the region file
   open(unit=11, file=REGFILE, form='formatted', status='old', action='read')

   ! loop for each region
   DO
      read(11, *, iostat=iostatus) reg

      IF (iostatus /= 0) THEN
         EXIT
      ENDIF

      ! get region file name and open nc file
      write(reg1, "(i4)") reg(1)
      write(reg2, "(i4)") reg(2)
      write(reg3, "(i4)") reg(3)
      write(reg4, "(i4)") reg(4)

      FILE_NAME = SRF_DIR//'RG_' &
         //trim(adjustL(reg1))//'_' &
         //trim(adjustL(reg2))//'_' &
         //trim(adjustL(reg3))//'_' &
         //trim(adjustL(reg4))//'.' &
         //DATASRC//trim(year)//'.nc'

      print *,">>> Processing file ", trim(FILE_NAME), "..."
      CALL check( nf90_open(trim(FILE_NAME), nf90_nowrite, ncid) )

      ! get the raw data
      CALL check( nf90_inq_varid(ncid, "LC"            , lc_vid      ) )
      CALL check( nf90_inq_varid(ncid, "MONTHLY_LC_LAI", lclai_vid   ) )
      CALL check( nf90_inq_varid(ncid, "MONTHLY_LC_SAI", lcsai_vid   ) )
      CALL check( nf90_inq_varid(ncid, "MONTHLY_LAI"   , lai_vid     ) )
      CALL check( nf90_inq_varid(ncid, "MONTHLY_SAI"   , sai_vid     ) )
      CALL check( nf90_inq_varid(ncid, "PCT_URBAN"     , purban_vid  ) )
      CALL check( nf90_inq_varid(ncid, "PCT_CROP"      , pcrop_vid   ) )
      CALL check( nf90_inq_varid(ncid, "PCT_WETLAND"   , pwetland_vid) )
      CALL check( nf90_inq_varid(ncid, "PCT_WATER"     , pwater_vid  ) )
      CALL check( nf90_inq_varid(ncid, "PCT_GLACIER"   , pglacier_vid) )
      CALL check( nf90_inq_varid(ncid, "PCT_PFT"       , ppft_vid    ) )
      CALL check( nf90_inq_varid(ncid, "HTOP"          , htop_vid    ) )

      CALL check( nf90_get_var(ncid, lc_vid   , lcdata) )
      CALL check( nf90_get_var(ncid, lclai_vid, lclai ) )
      CALL check( nf90_get_var(ncid, lcsai_vid, lcsai ) )
      CALL check( nf90_get_var(ncid, ppft_vid , ppft  ) )
      CALL check( nf90_get_var(ncid, lai_vid  , pftlai) )
      CALL check( nf90_get_var(ncid, sai_vid  , pftsai) )

      CALL check( nf90_get_var(ncid, purban_vid  , purban  ) )
      CALL check( nf90_get_var(ncid, pcrop_vid   , pcrop   ) )
      CALL check( nf90_get_var(ncid, pwetland_vid, pwetland) )
      CALL check( nf90_get_var(ncid, pglacier_vid, pglacier) )
      CALL check( nf90_get_var(ncid, pwater_vid  , pwater  ) )
      CALL check( nf90_get_var(ncid, htop_vid    , htop    ) )

      ! close file
      CALL check( nf90_close(ncid) )

      ! pre-process of raw data
      pcrop   (:,:) = pcrop   (:,:)/100.
      purban  (:,:) = purban  (:,:)/100.
      pwetland(:,:) = pwetland(:,:)/100.
      pglacier(:,:) = pglacier(:,:)/100.
      pwater  (:,:) = pwater  (:,:)/100.
      ppft  (:,:,:) = ppft  (:,:,:)/100.

      ! calculate the edge of small grids
      dll = (reg(4)-reg(2))*1./nxy  ! INTEGER to REAL
      DO i = 1, nxy
         lonw(i) = reg(2) + i*dll - dll
         lone(i) = reg(2) + i*dll
         latn(i) = reg(1) - i*dll + dll
         lats(i) = reg(1) - i*dll
      ENDDO

      ! calculate the area size of small grids
      DO i = 1, nxy
         dx = (lone(1)-lonw(1))*deg2rad
         dy = sin(latn(i)*deg2rad) - sin(lats(i)*deg2rad)
         sarea(:,i) = dx*dy*re*re
      ENDDO

      ! output grid edge
      dll = (reg(4)-reg(2))*1./nxy1  ! INTEGER to REAL
      DO i = 1, nxy1
         lonw1(i) = reg(2) + i*dll - dll
         lone1(i) = reg(2) + i*dll
         latn1(i) = reg(1) - i*dll + dll
         lats1(i) = reg(1) - i*dll
      ENDDO

      ! loop for output grids
!$OMP PARALLEL DO NUM_THREADS(92) &
!$OMP PRIVATE(k,i1,j1,io,jo,dx,dy,sumarea,i,j,lc,wgt,sumpct,ip,il,im) &
!$OMP PRIVATE(glai_lc,glai_pft,glai_pc,gsai_lc,gsai_pft,gsai_pc)
      DO k = 1, nxy1*nxy1

         ! calculate i1, j1
         i1 = (k-1)/nxy1 + 1
         j1 = mod(k-1, nxy1) + 1
         jo = int((reg(2) +180)/0.5) + j1
         io = int((90 - reg(1))/0.5) + i1

         ! calcualte area size
         dx = (lone1(j1)-lonw1(j1))*deg2rad
         dy = sin(latn1(i1)*deg2rad) - sin(lats1(i1)*deg2rad)
         sumarea = dx*dy*re*re

         ! loop for each small grid for aggregation
         DO i = (i1-1)*times+1, i1*times
            DO j = (j1-1)*times+1, j1*times

               ! PFT level aggregation
               pct_water  (jo,io) = pct_water  (jo,io) + pwater  (j,i)*sarea(j,i)
               pct_glacier(jo,io) = pct_glacier(jo,io) + pglacier(j,i)*sarea(j,i)
               pct_wetland(jo,io) = pct_wetland(jo,io) + pwetland(j,i)*sarea(j,i)
               pct_urban  (jo,io) = pct_urban  (jo,io) + purban  (j,i)*sarea(j,i)
               pct_crop   (jo,io) = pct_crop   (jo,io) + pcrop   (j,i)*sarea(j,i)

               ! total area include ocean
               area(jo,io) = area(jo,io) + sarea(j,i)

               lc = lcdata(j,i)
               IF (lc /= 0) THEN

                  ! aggregate land
                  pct_land(jo,io) = pct_land(jo,io) + sarea(j,i)

                  ! aggregate on lc level
                  ! need more check for htop_lc: 是否需要排除非植被区域
                  pct_lc (jo,io,lc)   = pct_lc (jo,io,lc)   + sarea(j,i)
                  lai_lc (jo,io,lc,:) = lai_lc (jo,io,lc,:) + sarea(j,i)*lclai(j,i,:)
                  sai_lc (jo,io,lc,:) = sai_lc (jo,io,lc,:) + sarea(j,i)*lcsai(j,i,:)
                  htop_lc(jo,io,lc)   = htop_lc(jo,io,lc)   + sarea(j,i)* htop(j,i)

                  ! aggregate on pc level
! yuan, 03/05/2022: a simple way to aggregate PC
                  DO ip = 1, npft
                     pct_pc(jo,io,ip,lc) = pct_pc(jo,io,ip,lc) + &
                        sarea(j,i)*ppft(j,i,ip)
                     lai_pc(jo,io,ip,lc,:) = lai_pc(jo,io,ip,lc,:) + &
                        sarea(j,i)*ppft(j,i,ip)*pftlai(j,i,ip,:)
                     sai_pc(jo,io,ip,lc,:) = sai_pc(jo,io,ip,lc,:) + &
                        sarea(j,i)*ppft(j,i,ip)*pftsai(j,i,ip,:)
                  ENDDO

                  ! aggregate on PFT level
                  ! 隐含的假设: 水体、冰川、城市、湿地和作物可由专业化的数据或高分辨率数据获得
                  ! 下面的代码具有兼容性，但目前并没有进行如上更新
                  ! exclude water, ice, urban, wetland, and crop
                  wgt = max(0., 1.-pwater(j,i)-pglacier(j,i)-purban(j,i)-pwetland(j,i)-pcrop(j,i))

                  ! adjust nature PFT's total percentage to 100%
                  sumpct = sum(ppft(j,i,1:(npft-1)))
                  IF (sumpct > 0.) THEN
                     ppft(j,i,1:(npft-1)) = ppft(j,i,1:(npft-1))/sumpct
                  ENDIF

                  DO ip = 1, npft-1
                     pct_pft(jo,io,ip) = pct_pft(jo,io,ip) + &
                        sarea(j,i)*wgt*ppft(j,i,ip)
                     htop_pft(jo,io,ip) = htop_pft(jo,io,ip) + &
                        sarea(j,i)*wgt*ppft(j,i,ip)*htop(j,i)
                     lai_pft(jo,io,ip,:) = lai_pft(jo,io,ip,:) + &
                        sarea(j,i)*wgt*ppft(j,i,ip)*pftlai(j,i,ip,:)
                     sai_pft(jo,io,ip,:) = sai_pft(jo,io,ip,:) + &
                        sarea(j,i)*wgt*ppft(j,i,ip)*pftsai(j,i,ip,:)
                  ENDDO

                  ! for crop (npft=16, is crop)
                  pct_pft(jo,io,npft) = pct_pft(jo,io,npft) + &
                     sarea(j,i)*ppft(j,i,npft)
                  htop_pft(jo,io,npft) = htop_pft(jo,io,npft) + &
                     sarea(j,i)*ppft(j,i,npft)*htop(j,i)
                  lai_pft(jo,io,npft,:) = lai_pft(jo,io,npft,:) + &
                     sarea(j,i)*ppft(j,i,npft)*pftlai(j,i,npft,:)
                  sai_pft(jo,io,npft,:) = sai_pft(jo,io,npft,:) + &
                     sarea(j,i)*ppft(j,i,npft)*pftsai(j,i,npft,:)
               ENDIF
            ENDDO
         ENDDO

         ! calculate LAI
         ! NOTE: Currently, LAI is not conserved to the LC LAIs
         ! ----------------------------------

         ! LC level
         DO il = 1, nlc
            IF (pct_lc(jo,io,il) > 0) THEN
                lai_lc(jo,io,il,:) = lai_lc(jo,io,il,:) / pct_lc(jo,io,il)
                sai_lc(jo,io,il,:) = sai_lc(jo,io,il,:) / pct_lc(jo,io,il)
               htop_lc(jo,io,il)   = htop_lc(jo,io,il)  / pct_lc(jo,io,il)
            ENDIF
         ENDDO

         ! PFT level
         DO ip = 1, npft
            IF (pct_pft(jo,io,ip) > 0) THEN
                lai_pft(jo,io,ip,:) =  lai_pft(jo,io,ip,:)/pct_pft(jo,io,ip)
                sai_pft(jo,io,ip,:) =  sai_pft(jo,io,ip,:)/pct_pft(jo,io,ip)
               htop_pft(jo,io,ip)   = htop_pft(jo,io,ip)  /pct_pft(jo,io,ip)
            ENDIF
         ENDDO

         ! PC level
         DO il = 1, nlc
            DO ip = 1, npft
               IF (pct_pc(jo,io,ip,il) > 0)  THEN
                  lai_pc(jo,io,ip,il,:) = lai_pc(jo,io,ip,il,:)/pct_pc(jo,io,ip,il)
                  sai_pc(jo,io,ip,il,:) = sai_pc(jo,io,ip,il,:)/pct_pc(jo,io,ip,il)
               ENDIF
            ENDDO
         ENDDO

         ! calculate fractional cover
         ! ----------------------------------

         ! PC level
         DO il = 1, nlc
            IF (pct_lc(jo,io,il) > 0) THEN
               pct_pc(jo,io,:,il) = pct_pc(jo,io,:,il)/pct_lc(jo,io,il) * 100.
            ENDIF
         ENDDO

         IF (pct_land(jo,io) > 0) THEN
            ! PFT level
            pct_pft    (jo,io,:) = pct_pft    (jo,io,:) / pct_land(jo,io) * 100.
            pct_urban  (jo,io)   = pct_urban  (jo,io)   / pct_land(jo,io) * 100.
            pct_wetland(jo,io)   = pct_wetland(jo,io)   / pct_land(jo,io) * 100.
            pct_crop   (jo,io)   = pct_crop   (jo,io)   / pct_land(jo,io) * 100.
            pct_glacier(jo,io)   = pct_glacier(jo,io)   / pct_land(jo,io) * 100.
            pct_water  (jo,io)   = pct_water  (jo,io)   / pct_land(jo,io) * 100.

            ! LC level
            pct_lc (jo,io,:) = pct_lc (jo,io,:) / pct_land(jo,io) * 100.
         ENDIF

         IF (abs(sumarea-area(jo,io))/sumarea > 1e-5) THEN
            print *, "Calculate area error! STOP!"
         ENDIF

         ! land fractional cover
         pct_land(jo,io) = pct_land(jo,io) / area(jo,io) * 100.

         ! check
         sumpct = sum(pct_pft(jo,io,:))+ &
            pct_urban(jo,io)+pct_wetland(jo,io)+ &
            pct_glacier(jo,io)+pct_water(jo,io)
         IF (sumpct > 1e-6 .and. abs(sumpct-100) > 1e-3) THEN
            print *, sumpct
            print *, pct_pft(jo,io,:)
            print *, sum(pct_pft(jo,io,:))
            print *, "Sum of pct_pft+urban+wetland+glacier+crop+water not equal 1! STOP!"
         ENDIF

         sumpct = sum(pct_lc(jo,io,:))
         IF (sumpct > 1e-6 .and. abs(sumpct-100) > 1e-3) THEN
            print *, sumpct
            print *, "Sum of pct_lc not equal 1! STOP!"
         ENDIF

         DO il = 1, nlc
            IF (DATASRC=="MOD" .and. (il==WATER .or. il==GLACIER .or. il==URBAN .or. il==WETLAND)) THEN
               CYCLE
            ENDIF
            IF (DATASRC=="ESA" .and. (il==WATER .or. il==GLACIER .or. il==URBAN)) THEN
               CYCLE
            ENDIF
            sumpct = sum(pct_pc(jo,io,:,il))
            IF (sumpct > 1e-6 .and. abs(sumpct-100) > 1e-3) THEN
               print *, sumpct
               print *, "Sum of pct_pc not equal 1! STOP!"
            ENDIF
         ENDDO

         ! make LAI/SAI conserved to the LC LAI/SAI
         DO im = 1, nmon

            ! initialize for grid lai from pc
            glai_pc = 0.
            gsai_pc = 0.

            glai_lc   = sum(pct_lc(jo,io,:) * lai_lc(jo,io,:,im))
            gsai_lc   = sum(pct_lc(jo,io,:) * sai_lc(jo,io,:,im))
            lai_grid(jo,io,im) = glai_lc/100.
            sai_grid(jo,io,im) = gsai_lc/100.
            glai_pft  = sum(pct_pft(jo,io,:)*lai_pft(jo,io,:,im))
            gsai_pft  = sum(pct_pft(jo,io,:)*sai_pft(jo,io,:,im))

            ! add LAI of urban and wetland
            IF (DATASRC=="MOD") THEN
               glai_pft = glai_pft + pct_lc(jo,io,URBAN  )*sum(pct_pc(jo,io,:,URBAN  )*lai_pc(jo,io,:,URBAN  ,im)/100.)
               gsai_pft = gsai_pft + pct_lc(jo,io,URBAN  )*sum(pct_pc(jo,io,:,URBAN  )*sai_pc(jo,io,:,URBAN  ,im)/100.)
               glai_pft = glai_pft + pct_lc(jo,io,WETLAND)*sum(pct_pc(jo,io,:,WETLAND)*lai_pc(jo,io,:,WETLAND,im)/100.)
               gsai_pft = gsai_pft + pct_lc(jo,io,WETLAND)*sum(pct_pc(jo,io,:,WETLAND)*sai_pc(jo,io,:,WETLAND,im)/100.)
            ENDIF

            ! add LAI of urban
            IF (DATASRC=="ESA") THEN
               glai_pft = glai_pft + pct_lc(jo,io,URBAN)*sum(pct_pc(jo,io,:,URBAN)*lai_pc(jo,io,:,URBAN,im)/100.)
               gsai_pft = gsai_pft + pct_lc(jo,io,URBAN)*sum(pct_pc(jo,io,:,URBAN)*sai_pc(jo,io,:,URBAN,im)/100.)
            ENDIF

! yuan, 01/05/2020: BUG!!!!
! 把il写成nlc,没有考虑pct_lc%,需要除以100
            DO il = 1, nlc
               glai_pc = glai_pc + pct_lc(jo,io,il)*sum(pct_pc(jo,io,:,il)*lai_pc(jo,io,:,il,im)/100.)
               gsai_pc = gsai_pc + pct_lc(jo,io,il)*sum(pct_pc(jo,io,:,il)*sai_pc(jo,io,:,il,im)/100.)
            ENDDO

            ! adjust PFT LAI/SAI
            IF (glai_pft > 0. .and. pct_land(jo,io)>1.) THEN
               lai_pft(jo,io,:,im) = lai_pft(jo,io,:,im) * min(2., glai_lc/glai_pft)
               WHERE (lai_pft(jo,io,:,im) > 10.) lai_pft(jo,io,:,im) = 10.
            ENDIF

            IF (gsai_pft > 0. .and. pct_land(jo,io)>1.) THEN
               sai_pft(jo,io,:,im) = sai_pft(jo,io,:,im) * min(2., gsai_lc/gsai_pft)
               WHERE (sai_pft(jo,io,:,im) > 3.) sai_pft(jo,io,:,im) = 3.
            ENDIF

            ! adjust pc LAI
            IF (glai_pc > 0. .and. pct_land(jo,io)>1.) THEN
               lai_pc(jo,io,:,:,im) = lai_pc(jo,io,:,:,im) * min(2., glai_lc/glai_pc)
               WHERE (lai_pc(jo,io,:,:,im) > 10.) lai_pc(jo,io,:,:,im) = 10.
            ENDIF

            IF (gsai_pc > 0. .and. pct_land(jo,io)>1.) THEN
               sai_pc(jo,io,:,:,im) = sai_pc(jo,io,:,:,im) * min(2., gsai_lc/gsai_pc)
               WHERE (sai_pc(jo,io,:,:,im) > 3.) sai_pc(jo,io,:,:,im) = 3.
            ENDIF

         ENDDO

      ENDDO  ! loop k
!$OMP END PARALLEL DO

   ENDDO

   ! create NC file
   FILE_NAME = OUT_DIR//'global_0.5x0.5.'//DATASRC//trim(year)//"_V5.nc"

   CALL check( nf90_create(FILE_NAME, NF90_NETCDF4, ncid) )

   ! Define the dimensions.
   CALL check( nf90_def_dim(ncid, "lat",  nyo , lat_dimid) )
   CALL check( nf90_def_dim(ncid, "lon",  nxo , lon_dimid) )
   CALL check( nf90_def_dim(ncid, "lc" ,  nlc , lc_dimid ) )
   CALL check( nf90_def_dim(ncid, "pft",  npft, pft_dimid) )
   CALL check( nf90_def_dim(ncid, "mon",  nmon, mon_dimid) )

   ! Define the coordinate variables.
   CALL check( nf90_def_var(ncid, "lat" , NF90_FLOAT, lat_dimid , lat_vid) )
   CALL check( nf90_def_var(ncid, "lon" , NF90_FLOAT, lon_dimid , lon_vid) )
   CALL check( nf90_def_var(ncid, "lc"  , NF90_INT  , lc_dimid  , lc_vid ) )
   CALL check( nf90_def_var(ncid, "pft" , NF90_INT  , pft_dimid , pft_vid) )
   CALL check( nf90_def_var(ncid, "mon" , NF90_INT  , mon_dimid , mon_vid) )

   ! Assign units attributes to coordinate variables.
   CALL check( nf90_put_att(ncid, lat_vid , "long_name", "Latitude"      ) )
   CALL check( nf90_put_att(ncid, lat_vid , "units"    , "degrees_north" ) )
   CALL check( nf90_put_att(ncid, lon_vid , "long_name", "Longitude"     ) )
   CALL check( nf90_put_att(ncid, lon_vid , "units"    , "degrees_east"  ) )
   CALL check( nf90_put_att(ncid, lc_vid  , "long_name", "LC index"      ) )
   CALL check( nf90_put_att(ncid, pft_vid , "long_name", "PFT index"     ) )
   CALL check( nf90_put_att(ncid, mon_vid , "long_name", "Month"         ) )
   CALL check( nf90_put_att(ncid, mon_vid , "units"    , "month"         ) )

   ! define output variables
   XY2D = (/ lon_dimid, lat_dimid /)
   CALL check( nf90_def_var(ncid, "AREA"       , NF90_FLOAT, XY2D, varea       , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "LANDFRAC"   , NF90_FLOAT, XY2D, vpct_land   , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "PCT_URBAN"  , NF90_FLOAT, XY2D, vpct_urban  , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "PCT_CROP"   , NF90_FLOAT, XY2D, vpct_crop   , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "PCT_WETLAND", NF90_FLOAT, XY2D, vpct_wetland, deflate_level=6) )
   CALL check( nf90_def_var(ncid, "PCT_GLACIER", NF90_FLOAT, XY2D, vpct_glacier, deflate_level=6) )
   CALL check( nf90_def_var(ncid, "PCT_WATER"  , NF90_FLOAT, XY2D, vpct_water  , deflate_level=6) )

   LC3D   = (/ lon_dimid, lat_dimid, lc_dimid  /)
   GRID3D = (/ lon_dimid, lat_dimid, mon_dimid /)
   PFT3D  = (/ lon_dimid, lat_dimid, pft_dimid /)
   CALL check( nf90_def_var(ncid, "PCT_LC"      , NF90_FLOAT, LC3D  , vpct_lc  , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "HTOP_LC"     , NF90_FLOAT, LC3D  , vhtop_lc , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "PCT_PFT"     , NF90_FLOAT, PFT3D , vpct_pft , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "HTOP_PFT"    , NF90_FLOAT, PFT3D , vhtop_pft, deflate_level=6) )
   CALL check( nf90_def_var(ncid, "MONTHLY_LAI" , NF90_FLOAT, GRID3D, vlai_grid, deflate_level=6) )
   CALL check( nf90_def_var(ncid, "MONTHLY_SAI" , NF90_FLOAT, GRID3D, vsai_grid, deflate_level=6) )

   PC4D  = (/ lon_dimid, lat_dimid, pft_dimid, lc_dimid  /)
   LC4D  = (/ lon_dimid, lat_dimid, lc_dimid , mon_dimid /)
   PFT4D = (/ lon_dimid, lat_dimid, pft_dimid, mon_dimid /)
   CALL check( nf90_def_var(ncid, "PCT_PC"         , NF90_FLOAT, PC4D , vpct_pc , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "MONTHLY_LC_LAI" , NF90_FLOAT, LC4D , vlai_lc , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "MONTHLY_LC_SAI" , NF90_FLOAT, LC4D , vsai_lc , deflate_level=6) )
   CALL check( nf90_def_var(ncid, "MONTHLY_PFT_LAI", NF90_FLOAT, PFT4D, vlai_pft, deflate_level=6) )
   CALL check( nf90_def_var(ncid, "MONTHLY_PFT_SAI", NF90_FLOAT, PFT4D, vsai_pft, deflate_level=6) )

   PC5D = (/ lon_dimid, lat_dimid, pft_dimid, lc_dimid, mon_dimid /)
   CALL check( nf90_def_var(ncid, "MONTHLY_PC_LAI", NF90_FLOAT, PC5D, vlai_pc, deflate_level=6) )
   CALL check( nf90_def_var(ncid, "MONTHLY_PC_SAI", NF90_FLOAT, PC5D, vsai_pc, deflate_level=6) )

   ! Assign units attributes to the netCDF variables.
   CALL check( nf90_put_att(ncid, varea       , "units"    , "km^2"                    ) )
   CALL check( nf90_put_att(ncid, varea       , "long_name", "Area of grid"            ) )
   CALL check( nf90_put_att(ncid, vpct_land   , "units"    , "%"                       ) )
   CALL check( nf90_put_att(ncid, vpct_land   , "long_name", "Land fractioanl coverage") )
   CALL check( nf90_put_att(ncid, vpct_urban  , "units"    , "%"                       ) )
   CALL check( nf90_put_att(ncid, vpct_urban  , "long_name", "Percent urban cover"     ) )
   CALL check( nf90_put_att(ncid, vpct_crop   , "units"    , "%"                       ) )
   CALL check( nf90_put_att(ncid, vpct_crop   , "long_name", "Percent crop cover"      ) )
   CALL check( nf90_put_att(ncid, vpct_wetland, "units"    , "%"                       ) )
   CALL check( nf90_put_att(ncid, vpct_wetland, "long_name", "Percent wetland cover"   ) )
   CALL check( nf90_put_att(ncid, vpct_glacier, "units"    , "%"                       ) )
   CALL check( nf90_put_att(ncid, vpct_glacier, "long_name", "Percent glacier cover"   ) )
   CALL check( nf90_put_att(ncid, vpct_water  , "units"    , "%"                       ) )
   CALL check( nf90_put_att(ncid, vpct_water  , "long_name", "Percent water cover"     ) )

   CALL check( nf90_put_att(ncid, vpct_lc  , "units"    , "%"                            ) )
   CALL check( nf90_put_att(ncid, vpct_lc  , "long_name", "Percent land cover type cover") )
   CALL check( nf90_put_att(ncid, vpct_pft , "units"    , "%"                            ) )
   CALL check( nf90_put_att(ncid, vpct_pft , "long_name", "Percent PFT cover"            ) )
   CALL check( nf90_put_att(ncid, vhtop_lc , "units"    , "m"                            ) )
   CALL check( nf90_put_att(ncid, vhtop_lc , "long_name", "LC tree top height"           ) )
   CALL check( nf90_put_att(ncid, vhtop_pft, "units"    , "m"                            ) )
   CALL check( nf90_put_att(ncid, vhtop_pft, "long_name", "PFT tree top height"          ) )

   CALL check( nf90_put_att(ncid, vpct_pc, "units"    , "%"                         ) )
   CALL check( nf90_put_att(ncid, vpct_pc, "long_name", "Percent extended PFT cover") )

   CALL check( nf90_put_att(ncid, vlai_grid, "units"    , "m^2/m^2"                 ) )
   CALL check( nf90_put_att(ncid, vlai_grid, "long_name", "Monthly GRID LAI values" ) )
   CALL check( nf90_put_att(ncid, vsai_grid, "units"    , "m^2/m^2"                 ) )
   CALL check( nf90_put_att(ncid, vsai_grid, "long_name", "Monthly GRID SAI values" ) )

   CALL check( nf90_put_att(ncid, vlai_lc  , "units"    , "m^2/m^2"                 ) )
   CALL check( nf90_put_att(ncid, vlai_lc  , "long_name", "Monthly LC LAI values"   ) )
   CALL check( nf90_put_att(ncid, vsai_lc  , "units"    , "m^2/m^2"                 ) )
   CALL check( nf90_put_att(ncid, vsai_lc  , "long_name", "Monthly LC SAI values"   ) )
   CALL check( nf90_put_att(ncid, vlai_pft , "units"    , "m^2/m^2"                 ) )
   CALL check( nf90_put_att(ncid, vlai_pft , "long_name", "Monthly PFT LAI values"  ) )
   CALL check( nf90_put_att(ncid, vsai_pft , "units"    , "m^2/m^2"                 ) )
   CALL check( nf90_put_att(ncid, vsai_pft , "long_name", "Monthly PFT SAI values"  ) )

   CALL check( nf90_put_att(ncid, vlai_pc, "units"    , "m^2/m^2"              ) )
   CALL check( nf90_put_att(ncid, vlai_pc, "long_name", "Monthly pc LAI values") )
   CALL check( nf90_put_att(ncid, vsai_pc, "units"    , "m^2/m^2"              ) )
   CALL check( nf90_put_att(ncid, vsai_pc, "long_name", "Monthly pc SAI values") )

   CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Title'  , Title  ) )
   CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Authors', Authors) )
   CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Adderss', Address) )
   CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Email'  , Email  ) )

   ! End define mode.
   CALL check( nf90_enddef(ncid) )

   ! Write the coordinate variable data. This will put the latitudes
   ! and longitudes of our data grid into the netCDF file.
   DO i = 1, nxo
      lonso(i) = -180. + i*dll - 0.5*dll
   ENDDO
   DO i = 1, nyo
      latso(i) =   90. - i*dll + 0.5*dll
   ENDDO
   CALL check( nf90_put_var(ncid, lat_vid,  latso      ) )
   CALL check( nf90_put_var(ncid, lon_vid,  lonso      ) )

   CALL check( nf90_put_var(ncid, lc_vid ,  lcs(1:nlc) ) )
   CALL check( nf90_put_var(ncid, pft_vid,  pfts       ) )
   CALL check( nf90_put_var(ncid, mon_vid,  months     ) )

   ! put variables
   CALL check( nf90_put_var(ncid, varea       , area       ) )
   CALL check( nf90_put_var(ncid, vpct_land   , pct_land   ) )
   CALL check( nf90_put_var(ncid, vpct_lc     , pct_lc     ) )
   CALL check( nf90_put_var(ncid, vlai_grid   , lai_grid   ) )
   CALL check( nf90_put_var(ncid, vsai_grid   , sai_grid   ) )
   CALL check( nf90_put_var(ncid, vlai_lc     , lai_lc     ) )
   CALL check( nf90_put_var(ncid, vsai_lc     , sai_lc     ) )
   CALL check( nf90_put_var(ncid, vhtop_lc    , htop_lc    ) )
   CALL check( nf90_put_var(ncid, vlai_pft    , lai_pft    ) )
   CALL check( nf90_put_var(ncid, vsai_pft    , sai_pft    ) )
   CALL check( nf90_put_var(ncid, vpct_pft    , pct_pft    ) )
   CALL check( nf90_put_var(ncid, vpct_crop   , pct_crop   ) )
   CALL check( nf90_put_var(ncid, vpct_urban  , pct_urban  ) )
   CALL check( nf90_put_var(ncid, vpct_wetland, pct_wetland) )
   CALL check( nf90_put_var(ncid, vpct_glacier, pct_glacier) )
   CALL check( nf90_put_var(ncid, vpct_water  , pct_water  ) )
   CALL check( nf90_put_var(ncid, vhtop_pft   , htop_pft   ) )
   CALL check( nf90_put_var(ncid, vpct_pc     , pct_pc     ) )
   CALL check( nf90_put_var(ncid, vlai_pc     , lai_pc     ) )
   CALL check( nf90_put_var(ncid, vsai_pc     , sai_pc     ) )

   ! Close the file. This causes netCDF to flush all buffers and make
   ! sure your data are really written to disk.
   CALL check( nf90_close(ncid) )

   print *,"*** SUCCESS writing file ", trim(FILE_NAME), "!"

   ! deallocate memory
   deallocate( lcdata       )
   deallocate( lclai        )
   deallocate( lcsai        )
   deallocate( ppft         )
   deallocate( pftlai       )
   deallocate( pftsai       )
   deallocate( purban       )
   deallocate( pcrop        )
   deallocate( pwetland     )
   deallocate( pglacier     )
   deallocate( pwater       )
   deallocate( htop         )
   deallocate( area         )
   deallocate( pct_land     )
   deallocate( pct_urban    )
   deallocate( pct_crop     )
   deallocate( pct_wetland  )
   deallocate( pct_glacier  )
   deallocate( pct_water    )
   deallocate( pct_lc       )
   deallocate( htop_lc      )
   deallocate( pct_pft      )
   deallocate( htop_pft     )
   deallocate( pct_pc       )
   deallocate( lai_grid     )
   deallocate( sai_grid     )
   deallocate( lai_lc       )
   deallocate( sai_lc       )
   deallocate( lai_pft      )
   deallocate( sai_pft      )
   deallocate( lai_pc       )
   deallocate( sai_pc       )

CONTAINS

   SUBROUTINE check(status)
      INTEGER, intent(in) :: status

      IF (status /= nf90_noerr) THEN
         print *, trim(nf90_strerror(status))
         STOP 2
      ENDIF
   END SUBROUTINE check

END PROGRAM MakeSurface
