MODULE mod_getMITgcm
  
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel

  IMPLICIT none
!
! Copied from readfields
!
contains
  ! ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  function get1dfield (gridfile, start1d, count1d)
    REAL, ALLOCATABLE,   DIMENSION(:)       :: get1dfield
    CHARACTER (len=200), INTENT(IN)         :: gridfile
    INTEGER, INTENT(IN)                     :: start1d, count1d
    INTEGER                                 :: d, rl, ierr
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
  function get2dfield (gridfile, start2d, count2d)
    REAL, ALLOCATABLE,   DIMENSION(:,:)     :: get2dfield
    CHARACTER (len=200), INTENT(IN)         :: gridfile
    INTEGER, DIMENSION(2), INTENT(IN)       :: start2d, count2d
    INTEGER,             DIMENSION(2)       :: d
    INTEGER                                 :: rl, ierr
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
  function get3dfield (gridfile, start3d, count3d)
    REAL, ALLOCATABLE,   DIMENSION(:,:,:)   :: get3dfield
    CHARACTER (len=200), INTENT(IN)         :: gridfile
    INTEGER, DIMENSION(3), INTENT(IN)       :: start3d, count3d
    INTEGER,             DIMENSION(3)       :: d
    INTEGER                                 :: rl, ierr
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

END MODULE mod_getMITgcm
