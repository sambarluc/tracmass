SUBROUTINE readfields

  USE mod_getMITgcm
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
  INTEGER, DIMENSION(3)                     :: start3d, count3d

  ! ===   ===   ===

  call datasetswap !Copy field(t+1) to field(t).

  timeStepNumber = intmin + (ints-intmin)*int(ngcm*3600.0/dtgcm)
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
  uvel     = get3dfield(gridfile, start3d, count3d)
  gridfile = trim(inDataDir)//'VVEL.'//fstamp//'.data'
  vvel     = get3dfield(gridfile, start3d, count3d)
  gridfile = trim(inDataDir)//'WVEL.'//fstamp//'.data'
  wvel   = get3dfield(gridfile, start3d, count3d)

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
    ! In MITgcm, uvel(1,:) holds the boundary condition
    ! while tracmass expects it at 0. Furthermore, tracmass has a different
    ! grid indexing than MITgcm, with velocities entering cell (i,j)
    ! having indices (i-1,j) and (i, j-1). dyu and dxv have been defined
    ! taking this index shift into account. Uvel and dzu are MITgcm-style,
    ! thus here we combine the two into the overall flux which must follow
    ! tracmass style.
     uflux(0:imt-1,:,km-k+1,2) = uvel(:,:,k)*dyu(0:imt-1,:)*dzu(:,:,k,1)
    !In MITgcm, vvel(:,1) holds the boundary condition
     vflux(:,0:jmt-1,km-k+1,2) = vvel(:,:,k)*dxv(:, 0:jmt-1)*dzv(:,:,k,1)
#ifdef explicit_w
     wflux(1:imt,1:jmt,km-k+1,2) = wvel(:,:,k)*dxdy
#endif
  end do kloop

  return
end subroutine readfields
