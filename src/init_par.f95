SUBROUTINE init_params
! Loads parameters from the namelists (projname).in
! Allocates matrices and set time variables.

   USE mod_param
   USE mod_seed
   USE mod_grid
   USE mod_name
   USE mod_time 
   USE mod_domain
   USE mod_vel
   USE mod_traj
 !  USE mod_dens
   USE mod_tempsalt
   USE mod_streamfunctions
   USE mod_tracer
   USE mod_getfile
   USE mod_write
   
#if defined diffusion || turb 
   USE mod_diffusion
#endif
#ifdef sediment
   USE mod_orbital
   USE mod_sed
#endif
   IMPLICIT NONE

!!----------------------------------------------------------------------------
   
   INTEGER                                    ::  gridVerNum

! Setup namelists
   namelist /INIT_NAMELIST_VERSION/ gridVerNum
   namelist /INIT_GRID_DESCRIPTION/ GCMname, GCMsource, gridName, gridSource,&
                                    griddesc, inDataDir!, topoDataDir !!joakim
   namelist /INIT_CASE_DESCRIPTION/ caseNum, caseDesc
   namelist /INIT_GRID_SIZE/        imt, jmt, km, nst, subGrid, subGridImin, &
                                    subGridImax, subGridJmin, subGridJmax,   &
                                    subGridKmin, subGridKmax, SubGridFile,   &
                                    subGridID, nperio
! Define the time parameters:
! ngcm    : number of hours between GCM output
! dtgcm   : time step in seconds of the GCM simulation
! intmin  : initial time step number of GCM, the simulation will start from this
! intspin : seeding time in steps (of output of GCM simulation)
! endJD   : number of days of the particle tracking simulation
! iter    : number of time steps between GCM output
! skipseed: seed every skipseed GCM outputs (if 1, seed all the time)
   namelist /INIT_TIME/             ngcm, dtgcm, intmin, intspin, endJD, iter, &
                                    skipseed
   namelist /INIT_WRITE_TRAJS/      twritetype, kriva, outDataDir, outDataFile, &
                                    outdircase, intminInOutFile, intpsi, outdirdate
          
   namelist /INIT_SEEDING/          nff, isec, idir, nqua, partQuant,        &
                                    loneparticle, SeedType, ist1,  &
                                    ist2, jst1, jst2, kst1, kst2, tst1, tst2,&
                                    seedDir, seedFile, varSeedFile, seedTime,&
                                    seedAll, seedPos, seedparts, seedpart_id,&
                                    seedsubints
   namelist /INIT_KILLZONES/        nend, ienw, iene, jens, jenn, timax
   namelist /INIT_TEMP_SALT/        tmin0, tmax0, smin0, smax0, rmin0, rmax0,&
                                    tmine, tmaxe, smine, smaxe, rmine, rmaxe
#if defined diffusion || turb 
   namelist /INIT_DIFFUSION/        ah, av
#endif
#ifdef sediment
   namelist /INIT_SEDIMENT/         partdiam, rhos, cwamp, twave, critvel
#endif

!!--------------------------------------------------------------------------  
   
   Project  = PROJECT_NAME
   Case     = CASE_NAME

   IF ((IARGC() > 0) )  THEN
      CALL getarg(1,project)
   END IF
   IF ((IARGC() > 1) )  THEN
      CALL getarg(2, Case)
   END IF
   
   if (len(trim(projdir)) == 0) then
      if (len(trim(ormdir)) .ne. 0) then
         projdir = trim(ormdir)//'/'//'projects/'//trim(Project)//'/'
      else
         projdir = 'projects/'//trim(Project)
      end if
   end if

   OPEN (8,file=trim(projdir)//'/'//trim(Case)//'.in',    &
        & status='OLD', delim='APOSTROPHE')
   ! -- Check if the namefiles has correct version number. 
   READ (8,nml=INIT_NAMELIST_VERSION)
   IF (gridVerNum < 6) THEN
      PRINT *,'                     ERROR                     '
      PRINT *,'Your namefile out of date. The latest version is described at:'
      PRINT *,'http://docs.tracmass.org/namelist.html'
      PRINT *,'Change gridVerNum to 6 when done.'
      STOP
   END IF
   READ (8,nml=INIT_GRID_DESCRIPTION)
   READ (8,nml=INIT_CASE_DESCRIPTION)
   READ (8,nml=INIT_GRID_SIZE)
   READ (8,nml=INIT_TIME)
   READ (8,nml=INIT_WRITE_TRAJS)
   READ (8,nml=INIT_SEEDING)
   READ (8,nml=INIT_KILLZONES)
   READ (8,nml=INIT_TEMP_SALT)
#if defined diffusion || turb 
   READ (8,nml=INIT_DIFFUSION)
#endif
#ifdef sediment
   READ (8,nml=INIT_SEDIMENT)   
#endif
   CLOSE (8)

   print *,'Run file    : ',trim(projdir)//'/'//trim(Case)//'.in'
   OPEN (8,file=trim(projdir)//'/'//trim(Case)//'.in',     &
        & status='OLD', delim='APOSTROPHE')
   READ (8,nml=INIT_NAMELIST_VERSION)
   IF (gridVerNum < 6) THEN
      PRINT *,'                     ERROR                     '
      PRINT *,'Your namefile out of date. The latest version is described at:'
      PRINT *,'http://docs.tracmass.org/namelist.html'
      PRINT *,'Change gridVerNum to 6 when done.'
      STOP
   END IF
   READ (8,nml=INIT_GRID_DESCRIPTION)
   READ (8,nml=INIT_CASE_DESCRIPTION)
   READ (8,nml=INIT_GRID_SIZE)
   READ (8,nml=INIT_TIME)
   READ (8,nml=INIT_WRITE_TRAJS)
   READ (8,nml=INIT_SEEDING)
   READ (8,nml=INIT_KILLZONES)
   READ (8,nml=INIT_TEMP_SALT)
#if defined diffusion || turb 
   READ (8,nml=INIT_DIFFUSION)
#endif
#ifdef sediment
   READ (8,nml=INIT_SEDIMENT)   
#endif
   CLOSE (8)
      
   SELECT CASE (subGrid)
   CASE (0)          
      PRINT *,'Sub-grid    : Use the Full grid.'     
      subGridImin =   1 
      subGridJmin =   1
      subGridKmin =   1
      subGridImax = imt
      subGridJmax = jmt 
      subGridKmax = km 
   CASE (1)
      PRINT *,'Sub-grid    : ', subGridImin ,subGridImax, &
           &   subGridJmin ,subGridJmax
      imt = subGridImax-subGridImin
      jmt = subGridJmax-subGridJmin

      
	if (subGridKmax == 0) subGridKmax = km      
#if !defined(explicit_w) && !defined(twodim)
      if ((subGridKmax-subGridKmin+1) < km) then
      	print *, subGridKmax-subGridKmin+1, km
         print *, 'ERROR!'
         print *, 'subGridKmin and subGridKmax requires -Dtwodim  or -Dexplicit_w'
         print *, 'to be selected in the project Makefile.'
         stop
      end if
#endif
      km  = subGridKmax-subGridKmin+1
   CASE default
      PRINT *,'==================== ERROR ===================='
      PRINT *,'This subGrid selection is not implemented yet.'
      PRINT *,'subGrid = ' ,subGrid
      STOP
   END SELECT
   start1d  = [subGridKmin]
   count1d  = [subGridKmax]
   start2d  = [1, 1,           subGridImin, subGridJmin]
   count2d  = [1, 1,           imt,         jmt        ]
   start3d  = [1, subGridImin, subGridJmin, subGridKmin]
   count3d  = [1, imt,         jmt,         km         ]
   
   IF (endJD < startJD) then
      PRINT *,'==================== ERROR ===================='
      PRINT *,'End JD must be positive, number of days to run.'
      PRINT *,'endJD = ' , endJD
      STOP
   END IF

#ifdef timeanalyt
      iter=1
#endif
   timax    =  24.*3600.*timax ! convert time lengths from days to seconds
   dstep    =  1.d0/dble(iter)
   ! number of GCM steps to run
   intrun   =  jd2ints(endJD - startJD)
   
   if (nff == 1) then
      intmax = intmin + intrun
   else
      intmax = intmin - intrun
   end if
   
   tseas = ngcm*3600.d0
   dtmin =  dstep * tseas

   ! --- ist -1 to imt ---
   IF ( ist1 == -1) THEN 
      ist1 = IMT
   END IF
   IF ( ist2 == -1) THEN 
      ist2 = IMT
   END IF
   ! --- jst -1 to jmt ---
   IF ( jst1 == -1) THEN
      jst1=jmt
   END IF
   IF ( jst2 == -1) THEN
      jst2=jmt
   END IF 
   ! --- kst -1 to km ---
   IF ( kst1 == -1) THEN
      kst1 = KM
   END IF
   IF ( kst2 == -1) THEN 
      kst2 = KM
   END IF

   call setup_outdatadir
      
   if (outDataFile == '')  outdataFile = Case

!!---------------------------------------------------------------------------
!!------------------------ A L L O C A T I O N S ----------------------------
!!---------------------------------------------------------------------------

      ! --- Allocate information about the coordinates and grid ---

      ALLOCATE ( csu (0:jmt), cst(jmt)  ) 
      ALLOCATE ( phi(0:jmt),   zlev(0:km) ) 
      ALLOCATE ( dyt(jmt), dxv(imt,0:jmt), dyu(0:imt,jmt) ) 
      ALLOCATE ( kmtb(imt, jmt) ) ! To allow for fractional bottom cells
      ALLOCATE ( mask(imt,jmt) )
      mask = 1
      dyt = 0
      dxv = 0
      dyu = 0

#if  zgrid3D
      ALLOCATE ( dzt(imt,jmt,km,nst) )
      dzt = 0
#endif /*zgrid3D*/
#ifdef varbottombox
      ALLOCATE ( dztb(imt,jmt,nst) )
#endif /*varbottombox*/
      ALLOCATE ( dxdy(imt,jmt) )   
      ALLOCATE ( kmt(imt,jmt), dz(km) )
    
      ! --- Allocate velocity fields, temperature, salinity, density, --- 
      ! --- sea-surface height, and trajectory data                   ---
      ALLOCATE ( vflux(imt,0:jmt,km,nst) )
      ALLOCATE ( uflux(0:imt,jmt,km,nst) )
      ALLOCATE ( hs(imt,jmt,nst) )
#if defined explicit_w || full_wflux
      ALLOCATE ( wflux(imt ,jmt ,0:km,NST) )
#else
      ALLOCATE ( wflux(0:km,NST) )
#endif

      hs    = 0.
      uflux = 0.
      vflux = 0.
      wflux = 0.d0
      ALLOCATE ( uvel(imt,jmt,km) ,vvel(imt,jmt,km) ,wvel(imt,jmt,km) )
      
      ! === Init mod_traj ===
      ntracmax = INT(FLOOR(MIN(DBLE(intspin), 24.0*endJD/ngcm) / skipseed))
      print*, "ntracmax: ", ntracmax
      ALLOCATE ( trj(NTRJ,ntracmax), nrj(NNRJ,ntracmax) )
      ALLOCATE ( nexit(NEND) ) 
      nrj = 0
      trj = 0.d0
      nexit = 0
      ntractot = 0
      numseedsubints = max(count(seedsubints /= -1), 1)
      if (nqua == 5 .AND. seedsubints(1)==-1) then
         print *,  "Error! "
         print *,  "At least one element in numseedsubints must be " // &
                   "given when nqua=5 is used." 
         stop
      elseif  (nqua /= 5) then
         seedsubints(1) = 0
      end if

#ifdef tempsalt
      ALLOCATE ( tem(imt,jmt,km,nst) ) 
      ALLOCATE ( sal(imt,jmt,km,nst) )
      ALLOCATE ( rho(imt,jmt,km,nst) )
      tem = 0.
      sal = 0.
      rho = 0.
#endif

      ! --- Allocate Lagrangian stream functions ---
#ifdef streamxy
      ALLOCATE ( stxyy(imt,jmt,nend), stxyx(imt,jmt,nend) )
      stxyy=0.
      stxyx=0.
#endif
#ifdef streamv
      ALLOCATE ( stxz(imt,km,nend), styz(jmt,km,nend) )
      stxz=0.
      styz=0.
#endif
#ifdef streamr
      ALLOCATE ( stxr(imt,mr,nend,lov), styr(jmt,mr,nend,lov), stzr(km,mr,nend,lov) )
      stxr=0.
      styr=0.
      stzr=0.
#endif
#ifdef stream_thermohaline
      ALLOCATE ( psi_ts(MR,MR,2,nend) )
      psi_ts=0.
#endif

      ! --- Allocate tracer data ---
#ifdef tracer
      ALLOCATE ( tra(imt,jmt,km) )
      tra=0.
#endif

      ! --- Allocate sedimentation data ---
#ifdef sediment
      ALLOCATE (orb(km) )
      nsed = 0
      nsusp = 0
#endif

END SUBROUTINE init_params

!!----------------------------------------------------------------------------
!!----------------------------------------------------------------------------
