SUBROUTINE readfields

  USE netcdf
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  
#ifdef tempsalt
  USE mod_dens
#endif
  
  IMPLICIT none
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  CHARACTER                                 :: dates(62)*17
  INTEGER, SAVE                             :: nread,ndates
  INTEGER                                   :: timeStepNumber
  CHARACTER(LEN=10), SAVE                   :: fstamp,rfilv,rfilh,rfilr
  logical around
  
  ! = Loop variables
  INTEGER                                   :: t ,i ,j ,k
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                       :: gridfile
  INTEGER                                   :: start1d, count1d
  INTEGER, DIMENSION(2)                     :: start2d, count2d
  INTEGER, DIMENSION(3)                     :: start3d, count3d
  INTEGER                                   :: ierr

  ! = ECCO Grid fields
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)     :: gridDRC ,gridDRF
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: gridDXC ,gridDXG
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: gridDYC ,gridDYG
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: hFacW   ,hFacS
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: gridRAC
  
  ! = Input fields from GCM
  !REAL,       ALLOCATABLE, DIMENSION(:,:,:) :: fieldx ,fieldy ,fieldw
 
  ! ===   ===   ===
  
  alloCondGrid: if(.not. allocated (gridDRC)) then
     allocate ( gridDRC(km)       ,gridDRF(km)       )
     allocate ( gridDXC(imt,jmt)  ,gridDXG(imt,jmt)  )
     allocate ( gridDYC(imt,jmt)  ,gridDYG(imt,jmt)  )
     allocate ( gridRAC(imt,jmt)                     )
     allocate ( hFacW(imt,jmt,km) ,hFacS(imt,jmt,km) )
  end if alloCondGrid
  
  !alloCondUVW: if(.not. allocated (fieldx)) then
  !   allocate ( fieldx(imt,jmt,km) ,fieldy(imt,jmt,km) )
  !   allocate ( fieldw(imt,jmt,km) )
  !end if alloCondUVW
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  call datasetswap !Copy field(t+1) to field(t).

  timeStepNumber = ints*72
  fstamp='0000000000'
  write (fstamp,'(i10.10)') timeStepNumber
  write (*,*) "filestamp ", fstamp
  
  ! === initialise ===
  !print *,'ints=',ints,intstart
  initCond: if(ints.eq.intstart) then
     ! call coordinat
     hs    = 0.
     uflux = 0.
     vflux = 0.
#ifdef tempsalt
     tem   = 0.
     sal   = 0.
     rho   = 0.
#endif
     ndates=0

     ! ======================================================
     !    ===  Set up the grid ===
     ! ======================================================
     ! Vertical beginning/end indices
     start1d  = 1
     count1d  = 4
     gridfile = trim(inDataDir) // 'DRC.data'
     gridDRC  = get1dfield()
     gridfile = trim(inDataDir) // 'DRF.data'
     gridDRF  = get1dfield()
     
     ! Horizontal beginning/end indices
     start2d  = [   1,  1]
     count2d  = [  60, 60]
     ! 3D beginning/end indices
     start3d  = [   1,  1, 1]
     count3d  = [  60, 60, 4]
     gridfile = trim(inDataDir) // 'DXC.data'
     gridDXC  = get2dfield()
     gridfile = trim(inDataDir) // 'RAC.data'
     gridRAC  = get2dfield()
     gridfile = trim(inDataDir) // 'DXG.data'
     gridDXG  = get2dfield()
     gridfile = trim(inDataDir) // 'DYC.data'
     gridDYC  = get2dfield()
     gridfile = trim(inDataDir) // 'DYG.data'
     gridDYG  = get2dfield() 
     gridfile = trim(inDataDir) // 'hFacW.data'
     hFacW    = get3dfield()
     gridfile = trim(inDataDir) // 'hFacS.data'
     hFacS    = get3dfield()

     dxdy     = gridRAC
     dz       = gridDRC(km:1:-1)
     kmt      = sum(ceiling(hFacW),3)
  endif initCond   ! === End init section ===
  
  start3d  = [   1,  1, 1]
  count3d  = [  60, 60, 4]
  gridfile = trim(inDataDir)//'UVEL.'//fstamp//'.data'
  uvel     = get3dfield()
  gridfile = trim(inDataDir)//'VVEL.'//fstamp//'.data'
  vvel     = get3dfield()
  gridfile = trim(inDataDir)//'WVEL.'//fstamp//'.data'
  wvel   = get3dfield()

  !CAwrite (*,*) "u-vel min/max", MINVAL(uvel), MAXVAL(uvel)
  !CAwrite (*,*) "v-vel min/max", MINVAL(vvel), MAXVAL(vvel)
  !CAwrite (*,*) "w-vel min/max", MINVAL(wvel), MAXVAL(wvel)
  !Density not included
  !SSH not included
  
  kloop: do k=1,km
     uflux(:,:,km-k+1,2)=uvel(:,:,k)*gridDYG*gridDRF(k)*hFacW(:,:,k)
     vflux(:,1:imt,km-k+1,2)=vvel(:,:,k)*gridDXG*gridDRF(k)*hFacS(:,:,k)
#ifdef explicit_w
     wflux(1:imt,1:jmt,km-k+1,2)=wvel(1:imt,1:jmt,k)*gridRAC(:,:)
#endif
  end do kloop

  return











  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  !    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  !    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###











contains
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  function get1dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:)       :: get1dfield
    INTEGER                                 :: d
    INTEGER                                 :: rl
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count1d+start1d-1
    allocate ( get1dfield(d) )
    rl=d*4
    open(unit=3001,file=gridfile, access='direct', recl=rl, iostat=ierr)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,gridfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    read(3001, rec=1) get1dfield
    close (3001)
  end function get1dfield

  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  
  function get2dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:)     :: get2dfield
    INTEGER,             DIMENSION(2)       :: d
    INTEGER                                 :: rl
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count2d+start2d-1
    allocate ( get2dfield(d(1),d(2)) )
    rl=product((d))*4
    open(unit=3001,file=gridfile, access='direct', recl=rl, iostat=ierr)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,gridfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    read(3001, rec=1) get2dfield
    close (3001)
  end function get2dfield

  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  
  function get3dfield ()
    REAL, ALLOCATABLE,   DIMENSION(:,:,:)   :: get3dfield
    INTEGER,             DIMENSION(3)       :: d
    INTEGER                                 :: rl
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
    d=count3d+start3d-1
    allocate ( get3dfield(d(1),d(2),d(3)) )
    rl=product((d))*4
    open(unit=3001,file=gridfile, access='direct', recl=rl, iostat=ierr)
    fileError: if(ierr.ne.0) then
       print *,'Error when trying to open the file'
       print *,'   ' ,gridfile
       print *,'    Error code: ' , ierr
       stop 3001
    end if fileError
    read(3001, rec=1) get3dfield
    close (3001)
  end function get3dfield

  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  
  !function get4dfield ()
  !  REAL, ALLOCATABLE,   DIMENSION(:,:,:,:) :: get4dfield
  !  INTEGER,             DIMENSION(4)       :: d
  !end function get4dfield
end subroutine readfields
