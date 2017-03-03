SUBROUTINE setupgrid
  
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel

  IMPLICIT none
  ! =============================================================
  !    ===  Set up the grid ===
  ! =============================================================
  ! Subroutine for defining the grid of MITgcm.
  ! Run once before the loop starts.
  ! -------------------------------------------------------------
  ! The following arrays has to be populated:
  !
  !  dxdy - Horizontal area of cells (T points)
  !  dz   - Thickness of standard level (T point) 
  !  dzt  - Time-invariant thickness of level (T point)
  !  dzu  - Time-invariant thickness of level (U point)
  !  dzv  - Time-invariant thickness of level (V point)
  !  kmt  - Number of levels from surface to seafloor (T point)
  !  kmu  - Number of levels from surface to seafloor (U point)
  !  kmv  - Number of levels from surface to seafloor (V point)
  ! -------------------------------------------------------------


  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
  ! = Variables for filename generation 
  INTEGER, SAVE                             :: nread,ndates
  CHARACTER(LEN=10), SAVE                   :: fstamp
  
  ! = Loop variables
  INTEGER                                   :: i ,j ,k
  
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                       :: gridfile
  INTEGER                                   :: start1d, count1d
  INTEGER, DIMENSION(2)                     :: start2d, count2d
  INTEGER, DIMENSION(3)                     :: start3d, count3d
  INTEGER                                   :: ierr

  ! = MITgcm Grid fields
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: hFacW   ,hFacS, hFacC
  
  allocate ( hFacW(imt,jmt,km) ,hFacS(imt,jmt,km) )
  allocate ( hFacC(imt,jmt,km)                    )
  
  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  ! ======================================================
  !    ===  Set up the grid ===
  ! ======================================================
  ! Vertical beginning/end indices
  start1d  = 1
  count1d  = 4
  ! Horizontal beginning/end indices
  start2d  = [   1,  1]
  count2d  = [  60, 60]
  ! 3D beginning/end indices
  start3d  = [   1,  1, 1]
  count3d  = [  60, 60, 4]

  ! load grid
  gridfile = trim(inDataDir) // 'DRF.data'
  dz       = get1dfield()
  dz       = dz(km:1:-1)
  gridfile = trim(inDataDir) // 'RAC.data'
  dxdy     = get2dfield()
  gridfile = trim(inDataDir) // 'DXG.data'
  dxv      = get2dfield()
  gridfile = trim(inDataDir) // 'DYG.data'
  dyu      = get2dfield() 
  gridfile = trim(inDataDir) // 'hFacW.data'
  hFacW    = get3dfield()
  gridfile = trim(inDataDir) // 'hFacS.data'
  hFacS    = get3dfield()
  gridfile = trim(inDataDir) // 'hFacC.data'
  hFacC    = get3dfield()

  kmt      = sum(ceiling(hFacC),3)
  allocate ( kmu(imt,jmt), kmv(imt,jmt) )
  kmu      = sum(ceiling(hFacW),3)
  kmv      = sum(ceiling(hFacS),3)

  allocate ( dzu(imt,jmt,km,1),dzv(imt,jmt,km,1) )

  kloop: do k=1,km
     dzt(:,:,km-k+1,1) = dz(k)*hFacC(:,:,k)
     dzu(:,:,km-k+1,1) = dz(k)*hFacW(:,:,k)
     dzv(:,:,km-k+1,1) = dz(k)*hFacS(:,:,k)
  end do kloop
  dzt(:,:,:,2) = dzt(:,:,:,1)

  !
  ! Ensure thickness is zero in invalid points
  !
  do k=1,km
     where (k .le. (km-kmt(1:imt,1:jmt)))
        dzt(:,:,k,1) = 0
        dzu(:,:,k,1) = 0
        dzv(:,:,k,1) = 0
     end where
  enddo 



!
! Copied from readfields
!
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

end SUBROUTINE setupgrid
