SUBROUTINE readfields

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
  INTEGER, SAVE                             :: nread,ndates
  INTEGER                                   :: timeStepNumber
  CHARACTER(LEN=10), SAVE                   :: fstamp
  
  ! = Loop variables
  INTEGER                                   :: t ,i ,j ,k
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                       :: gridfile
  INTEGER                                   :: start1d, count1d
  INTEGER, DIMENSION(2)                     :: start2d, count2d
  INTEGER, DIMENSION(3)                     :: start3d, count3d
  INTEGER                                   :: ierr

  ! ===   ===   ===

  call datasetswap !Copy field(t+1) to field(t).

  write(*,*) "ints: ", ints
  timeStepNumber = ints*180 ! 180 is the number of MITgcm timesteps in one output unit
  fstamp='0000000000'
  write (fstamp,'(i10.10)') timeStepNumber
  
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
   endif initCond
  
  start3d  = [   1,  1, 1]
  count3d  = [  60, 60, 4]
  gridfile = trim(inDataDir)//'UVEL.'//fstamp//'.data'
  uvel     = get3dfield()
  gridfile = trim(inDataDir)//'VVEL.'//fstamp//'.data'
  vvel     = get3dfield()
  gridfile = trim(inDataDir)//'WVEL.'//fstamp//'.data'
  wvel   = get3dfield()

  write (*,*) '<file>.'//fstamp//'.data'
  !CAwrite (*,*) "u-vel min/max", MINVAL(uvel), MAXVAL(uvel)
  !CAwrite (*,*) "v-vel min/max", MINVAL(vvel), MAXVAL(vvel)
  !CAwrite (*,*) "w-vel min/max", MINVAL(wvel), MAXVAL(wvel)
  !CAwrite (*,*) "dyu min/max", MINVAL(dyu), MAXVAL(dyu)
  !CAwrite (*,*) "dxv min/max", MINVAL(dxv), MAXVAL(dxv)
  !CAwrite (*,*) "dzu min/max", MINVAL(dzu), MAXVAL(dzu)
  !CAwrite (*,*) "dzv min/max", MINVAL(dzv), MAXVAL(dzv)
  !Density not included
  !SSH not included
  
  kloop: do k=1,km
     uflux(:,:,km-k+1,2)     = uvel(:,:,k)*dyu*dzu(:,:,k,1)
     vflux(:,1:imt,km-k+1,2) = vvel(:,:,k)*dxv*dzv(:,:,k,1)
#ifdef explicit_w
     wflux(1:imt,1:jmt,km-k+1,2) = wvel(:,:,k)*dxdy
#endif
  end do kloop

  return



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

end subroutine readfields
