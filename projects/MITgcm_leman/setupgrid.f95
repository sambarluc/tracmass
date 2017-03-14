SUBROUTINE setupgrid
  
  USE mod_getMITgcm
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
  !  kmtb - Number of fractional levels from surface to seafloor
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
  count1d  = km
  ! Horizontal beginning/end indices
  start2d  = [   1,  1]
  count2d  = [ imt,jmt]
  ! 3D beginning/end indices
  start3d  = [   1,  1, 1]
  count3d  = [ imt,jmt,km]

  ! load grid
  gridfile = trim(inDataDir) // 'DRF.data'
  dz       = get1dfield(gridfile, start1d, count1d)
  gridfile = trim(inDataDir) // 'RAC.data'
  dxdy     = get2dfield(gridfile, start2d, count2d)
  gridfile = trim(inDataDir) // 'DXG.data'
  dxv(:,0:jmt-1) = get2dfield(gridfile, start2d, count2d)
  gridfile = trim(inDataDir) // 'DYG.data'
  dyu(0:imt-1,:) = get2dfield(gridfile, start2d, count2d) 
  gridfile = trim(inDataDir) // 'hFacW.data'
  hFacW    = get3dfield(gridfile, start3d, count3d)
  gridfile = trim(inDataDir) // 'hFacS.data'
  hFacS    = get3dfield(gridfile, start3d, count3d)
  gridfile = trim(inDataDir) // 'hFacC.data'
  hFacC    = get3dfield(gridfile, start3d, count3d)

  kmt      = sum(ceiling(hFacC),3)
  kmtb     = sum(hFacC,3)

  allocate ( dzu(imt,jmt,km,1),dzv(imt,jmt,km,1) )

  ! Change the ordering of the vertical axis
  kloop: do k=1,km
     dzt(:,:,km-k+1,1) = dz(k)*hFacC(:,:,k)
     dzu(:,:,km-k+1,1) = dz(k)*hFacW(:,:,k)
     dzv(:,:,km-k+1,1) = dz(k)*hFacS(:,:,k)
  end do kloop
  dzt(:,:,:,2) = dzt(:,:,:,1)

  ! Invert order of dz to be consistent with tracmass ordering
  dz       = dz(km:1:-1)

  return
end SUBROUTINE setupgrid
