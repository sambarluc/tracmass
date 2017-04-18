PROGRAM TRACMASS
!
! Some details on the grid definitions of tracmass
! The grid indexing for the C-cell is:
!
!      -------v(i,j)-------      y=j
!      |                  |
!      |                  |
!      |                  |
!    u(i-1,j)  T(i,j)   u(i,j)   y=i-0.5
!      |                  |
!      |                  |
!      |                  |
!      -------v(i,j-1)-----      y=i-1
!
!    x=i-1    x=i-0.5    x=i
!
! Note this is not the same as in MITgcm, so when reading GCM
! fields care must be taken in order to load the GCM field at
! the right indices.
!
  USE mod_seed, only:  nqua, nff, num
  USE mod_grid, only:  dyu, dxv
  USE mod_time, only:  intmin, intstart, intend, intspin, intrun, intmax, &
       partQuant, voltr
  USE mod_write, only: open_outfiles, close_outfiles
  USE mod_print, only: print_header_main, writesetup_main
#ifdef diffusion
  USE mod_diffusion
#endif
  
  IMPLICIT none
  INTEGER                                    :: i,j,n

  call print_header_main
  call init_params
  call writesetup_main
  
  if(nff == 1) then ! forward 
     intstart =  intmin          
     intend   =  intmax
  elseif(nff == 2) then ! backward
     intstart =  intmin
     intend   =  intmax
     intspin  = -intspin
     intrun   = -intrun
     nff      =  -1
  end if

  call setupgrid

  if (minval(dxv) < 0) then
     print *, " "
     print *, " === Error! === "
     print *, "The array dxv contains negative values."
     print *, "Please check your setupgrid.f95 file."
     stop
  end if
  if (minval(dyu) < 0) then
     print *, " "
     print *, " === Error! === "
     print *, "The array dyu contains negative values."
     print *, "Please check your setupgrid.f95 file."
     stop
  end if

  call init_seed

  if(nqua.eq.1) then ! number of trajectories (per time resolution)
     ! num=NTRACMAX
     num = partQuant
  elseif(nqua.eq.2) then 
     voltr = partQuant 
  elseif(nqua.eq.3) then 
     voltr = partQuant
  end if
  
  call open_outfiles
  call loop
  call close_outfiles

end PROGRAM TRACMASS
