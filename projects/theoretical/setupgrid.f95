SUBROUTINE setupgrid
  
  USE mod_param
  USE mod_vel
  
  USE mod_time
  USE mod_grid
  USE mod_name
  USE mod_vel
  
  USE mod_getfile
  
  IMPLICIT none
  ! =============================================================
  !    ===  Set up the grid ===
  ! =============================================================
  ! Subroutine for defining the grid of the GCM. Run once
  ! before the loop starts.
  ! -------------------------------------------------------------
  ! The following arrays has to be populated:
  !
  !  dxdy - Area of horizontal cell-walls.
  !  dz   - Height of k-cells in 1 dim. |\
  !  dzt  - Height of k-cells i 3 dim.  |- Only one is needed
  !  kmt  - Number of k-cells from surface to seafloor.
  !
  ! The following might be needed to calculate
  ! dxdy, uflux, and vflux
  !
  !  dzu - Height of each u-gridcell.
  !  dzv - Height of each v-gridcell.
  !  dxu -
  !  dyu -
  ! -------------------------------------------------------------



  ! === Init local variables for the subroutine ===
  INTEGER                                    :: i ,j ,k ,kk

  integer :: idx, idy, idz, idlon, idlat, iddep, iddx, iddy, iddz
  integer, dimension(1) :: dim1d
  integer, dimension(2) :: dim2d
  integer, dimension(3) :: dim3d
  real, allocatable, dimension(:,:) :: lon, lat
  real, allocatable, dimension(:) :: depth
  logical :: lwrite_nc, lread_nc

lwrite_nc = .true.
lread_nc  = .true.

allocate( lon(imt,jmt), lat(imt,jmt) )
allocate( depth(km) )

kmt=KM ! flat bottom

!dxdeg=dx*deg
!dydeg=dy*deg

! Nicoletta Fabboni velocities, which have analytical solutions
dxv (:,:) = 250. 
dyu (:,:) = dxv(:,:)
dxdy(:,:) = dxv(:imt,:) * dyu(:imt,:)
dz  (:)   = 10.
dzt(:,:,:,:) = 10.

mask(:,:) = 1

do j = 1, jmt
   do i = 1, imt
      lon(i,j) = i * dxv(i,j)
      lat(i,j) = j * dyu(i,j)
   end do 
end do

depth(1) = dz(1)/2.
do k = 2, km
   depth(k) = SUM(dz(1:k-1)) + dz(k)/2.
end do


end SUBROUTINE setupgrid
