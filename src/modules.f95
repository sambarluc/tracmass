

MODULE mod_precdef		! Precision definitions
   !integer, parameter                       :: P4 = selected_real_kind(6, 37)
   integer, parameter                       :: DP = selected_real_kind(15, 307)
   integer, parameter                       :: QP = selected_real_kind(33, 4931)
ENDMODULE mod_precdef


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_param
  USE mod_precdef
  INTEGER                                   :: jmax, ntracmax
  INTEGER, PARAMETER                        :: MR=501 ! or 1001
  INTEGER                                   :: ncoor,kriva,iter,ngcm
  REAL(DP)                                  :: dtgcm
  REAL(DP), PARAMETER                       :: UNDEF=1.d20 
  REAL(DP), PARAMETER                       :: EPS=1.d-7 ! the small epsilon

  REAL(DP), PARAMETER                       :: grav = 9.81
  REAL(DP), PARAMETER                       :: PI = 3.14159265358979323846d0
  REAL(DP), PARAMETER                       :: radius = 6371229.d0 
  REAL(DP), PARAMETER                       :: radian = pi/180.d0  
  REAL(DP), PARAMETER                       :: deg=radius*radian   
  REAL(DP), PARAMETER                       :: tday=24.d0 * 3600.d0
ENDMODULE mod_param
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_loopvars
  USE mod_precdef
  REAL(DP)                                  :: ds, dsmin
  REAL(DP)                                  :: dse, dsw, dsn, dss
  REAL(DP)                                  :: dsu, dsd, dsc
  LOGICAL                                   :: scrivi
  INTEGER                                   :: niter
  REAL(DP)                                  :: ss0
  INTEGER                                   :: lbas
  REAL(DP)                                  :: subvol
ENDMODULE mod_loopvars
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

MODULE mod_traj
  USE mod_precdef
  ! Variables connected to particle positions.
  INTEGER, PARAMETER                        :: NNRJ=8, NTRJ=7
  INTEGER                                   :: nend
  INTEGER                                   :: ntrac, ntractot=0
  ! === Particle arrays ===
  REAL(DP), ALLOCATABLE,  DIMENSION(:,:)    :: trj
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: nrj 
  ! === Particle counters ===
  INTEGER                                   :: nout=0, nloop=0, nerror=0, nrh0=0
  INTEGER, ALLOCATABLE,DIMENSION(:)         :: nexit
  ! === Particle positions ===
  INTEGER                                   :: ia, ja, ka, iam
  INTEGER                                   :: ib, jb, kb, ibm
  REAL(DP)                                  :: x0, y0, z0
  REAL(DP)                                  :: x1, y1, z1
ENDMODULE mod_traj
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_tempsalt
  USE mod_precdef
  REAL(DP)                                  :: rmin, tmin, smin
  REAL(DP)                                  :: rmax, tmax, smax
  REAL(DP)                                  :: dr ,dtemp ,dsalt
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:) :: tem,sal,rho
  REAL*4                                  :: tmin0 ,tmax0
  REAL*4                                  :: smin0 ,smax0
  REAL*4                                  :: rmin0 ,rmax0
  REAL*4                                  :: tmine ,tmaxe
  REAL*4                                  :: smine ,smaxe
  REAL*4                                  :: rmine ,rmaxe
end MODULE mod_tempsalt

  ! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_grid
  USE mod_param, only: pi, undef, iter
  USE mod_precdef
  IMPLICIT NONE

  INTEGER                                   :: imt, jmt, km
  INTEGER                                   :: nst=2
  INTEGER                                   :: nsm=1,  nsp=2
  INTEGER                                   :: wnsm=1, wnsp=2
  INTEGER                                   :: nperio=1
  REAL(DP)                                  :: dx,dy
  REAL(DP)                                  :: dxdeg,dydeg,stlon1,stlat1
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)   :: hs
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)   :: botbox
  REAL*4, ALLOCATABLE, DIMENSION(:,:)     :: dxv, dyu, ang
  REAL(DP), ALLOCATABLE, DIMENSION(:)       :: dz
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)     :: dxdy
  REAL(DP)                                  :: dxyz
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: mask
  REAL(DP), ALLOCATABLE, DIMENSION(:)       :: csu,cst,dyt,phi

  ! === Vertical grids ===
  REAL(DP), ALLOCATABLE, DIMENSION(:)       :: zlev
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: z_r, z_w
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)     :: dzt, dzu, dzv
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dzt0, dzu0, dzv0
  REAL, ALLOCATABLE, DIMENSION(:,:)         :: dzt0surf,dzu0surf,dzv0surf
#ifdef varbottombox 
  REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: dztb
#endif /*varbottombox*/
#ifdef ifs
  REAL(DP), ALLOCATABLE, DIMENSION(:)       :: aa, bb
#endif
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: kmt, kmu, kmv
  INTEGER                                   :: subGrid     ,subGridID
  INTEGER                                   :: subGridImin ,subGridImax
  INTEGER                                   :: subGridJmin ,subGridJmax
  INTEGER                                   :: subGridKmin=1 ,subGridKmax=0
  CHARACTER(LEN=200)                        :: SubGridFile 
  INTEGER                                   :: degrade_space=0

#ifdef ifs
  REAL(DP), PARAMETER                       :: R_d = 287.05d0
  REAL(DP), PARAMETER                       :: L_v = 2.5d0 * 1e+6   
  REAL(DP), PARAMETER                       :: c_d = 1004.d0
#endif

CONTAINS
  function l2d(lon1,lon2,lat1,lat2)
    real                                    :: lon1,lon2,lat1,lat2,l2d
    real                                    :: rlon1,rlon2,rlat1,rlat2
    real                                    :: dlon,dlat,a,c
    dlon = (lon2 - lon1)/180*pi
    rlat1 = lat1 /180.*pi
    rlat2 = lat2 /180.*pi
    dlat = rlat2 - rlat1
    a = (sin(dlat/2))**2 + cos(rlat1) * cos(rlat2) * (sin(dlon/2))**2
    c = 2 * asin(min(1.0,sqrt(a)))
    l2d = 6367 * c * 1000
  end function l2d

  subroutine calc_dxyz(intrpr, intrpg)

    use mod_traj, only: ib,jb,kb
    IMPLICIT NONE
    REAL(DP)                                    :: intrpr, intrpg

    ! T-box volume in m3
#ifdef zgrid3D
    dxyz = intrpg * dzt(ib,jb,kb,nsp) + intrpr * dzt(ib,jb,kb,nsm)
#else
    dxyz =dz(kb)
#ifdef varbottombox
    if(kb == KM+1-kmt(ib,jb) ) dxyz=dztb(ib,jb,1)
#endif /*varbottombox*/
#ifdef freesurface
    if(kb == KM) dxyz=dxyz+intrpg*hs(ib,jb,nsp)+intrpr*hs(ib,jb,nsm)
#endif /*freesurface*/
#endif /*zgrid3D*/
    dxyz=dxyz*dxdy(ib,jb)
    if (dxyz<0) then
       print *,'=========================================================='
       print *,'ERROR: Negative box volume                                '
       print *,'----------------------------------------------------------'
       print *,'dxdy = ', dxdy(ib,jb)
       print *,'ib  = ', ib, ' jb  = ', jb, ' kb  = ', kb 
       print *,'----------------------------------------------------------'
       print *,'The run is terminated'
       print *,'=========================================================='
       stop
    end if
  end subroutine calc_dxyz

ENDMODULE mod_grid
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_time
  USE mod_precdef
  ! Variables and routines for timekeeping

  !Timestep increasing with one for each new velocity field
  INTEGER                                   :: ints      ,intstart ,intend
  INTEGER                                   :: intrun    ,intspin
  INTEGER                                   :: intmin    ,intmax
  INTEGER                                   :: nff=1

  ! === JD when the run starts
  REAL(DP)                                  :: startJD=0, ttpart
  REAL(DP)                                  :: endJD=-999
  ! === Current JD
  REAL(DP)                                  :: currJDtot
  ! === Looping time
  INTEGER                                   :: loopints, loopintstart
  REAL(DP)                                  :: loopJD
  INTEGER                                   :: loopDay
  INTEGER                                   :: loopHour, loopMin, loopSec

  INTEGER*8                                 :: ntime
  ! === Time-interpolation variables in loop ===
  REAL(DP)                                  :: dt, t0
  REAL(DP)                                  :: dtreg
  REAL(DP)                                  :: tseas, tyear, dtmin,voltr
  REAL(DP)                                  :: tstep, dstep, tss, partQuant
  REAL(DP)                                  :: ts, tt
  REAL(DP)                                  :: intrpr, intrpg
  REAL(DP)                                  :: intrpb, intrpbg
CONTAINS

  subroutine updateClock  
    USE mod_param, only: ngcm
    IMPLICIT NONE
    ttpart = anint((anint(tt,8)/tseas-floor(anint(tt,8)/tseas))*tseas)/tseas

    currJDtot = (ints+ttpart)*(dble(ngcm)/24.)
    
    loopints = ints - intstart

    loopJD = (loopints + ttpart)*(dble(ngcm)/24)
    loopDay  = int(loopJD)
    loopJD = (loopJD - dble(loopDay)) * 24.0
    loopHour = int(loopJD,8)
    loopJD = (loopJD - dble(loopHour)) * 60
    LoopMin  = int(loopJD,8)
    loopSec  = int((loopJD - dble(loopMin)) * 60,8)
  end subroutine updateClock

  INTEGER function jd2ints(jd)
    USE mod_param, only: ngcm
    REAL(DP) :: jd
    jd2ints = nint((jd)/(real(ngcm)/24.))
    return
  end function jd2ints
  
  subroutine calc_time
    USE mod_loopvars, only: ds, dsc, dsmin
    use mod_grid, only: dxyz
    USE mod_param, only: iter
    IMPLICIT NONE

#ifdef regulardt
    if(ds == dsmin) then ! transform ds to dt in seconds
       !            dt=dt  ! this makes dt more accurate
    else
       dt = ds * dxyz 
    endif
#else
    if(ds == dsmin) then ! transform ds to dt in seconds
       dt=dtmin  ! this makes dt more accurate
    else
       dt = ds * dxyz 
    endif
#endif /*regulardt*/
    if(dt.lt.0.d0) then
       Print *,"Error! dt is less than zero."
       print *,'dt=',dt,"ds=",ds,"dxyz=",dxyz,"dsmin=",dsmin
       stop 4968
    endif
    ! === if time step makes the integration ===
    ! === exceed the time when fields change ===
    if(tss+dt/tseas*dble(iter).ge.dble(iter)) then
       dt=dble(int(ts,8)+1)*tseas-tt
       tt=dble(int(ts,8)+1)*tseas
       ts=dble(int(ts,8)+1)
       tss=dble(iter)
       ds=dt/dxyz
       dsc=ds
    else
       tt=tt+dt
#if defined regulardt
       if(dt == dtmin) then
          ts=ts+dstep
          tss=tss+1.d0
       elseif(dt == dtreg) then  
          ts=nint((ts+dtreg/tseas)*dble(iter), 8)/dble(iter)
          !                 ts=ts+dtreg/tseas
          tss=dble(nint(tss+dt/dtmin))
       else
          ts=ts+dt/tseas
          tss=tss+dt/dtmin
       endif
#else
       if(dt == dtmin) then
          ts=ts+dstep
          tss=tss+1.d0
       else
          ts =ts +dt/tseas
          tss=tss+dt/tseas*dble(iter)
          !                 tss=tss+dt/dtmin
       endif
#endif /*regulardt*/
    end if
    ! === time interpolation constant ===
#if defined fixedtimestep 
    intrpbg=0.d0 ! mimics Ariane's lack of linear interpolation of the velocity fields
#else    
    intrpbg=dmod(ts,1.d0) 
#endif
    intrpb =1.d0-intrpbg
  end subroutine calc_time

ENDMODULE mod_time
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_buoyancy
  USE mod_precdef
ENDMODULE mod_buoyancy
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_domain
  USE mod_precdef
  INTEGER, DIMENSION(10)                :: ienw ,iene, jens ,jenn
  REAL*4                                :: timax
ENDMODULE mod_domain
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===        
MODULE mod_dens
  USE mod_precdef
ENDMODULE mod_dens
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===        


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_vel
  USE mod_grid, only: nsm, nsp, dzt
  USE mod_precdef
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)    :: uflux, vflux
#if defined explicit_w || full_wflux
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:)    :: wflux
#else
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)        :: wflux
#endif
  REAL,   ALLOCATABLE, DIMENSION(:,:,:)      :: uvel ,vvel ,wvel 
  REAL(DP)                                     :: ff
  INTEGER                                    :: degrade_time=0
    integer, save                            :: degrade_counter = 0

CONTAINS
 
  subroutine datasetswap

    USE  mod_grid, only      : nsm,nsp,hs
    USE  mod_tempsalt, only  : tem,sal,rho

    IMPLICIT NONE

    hs(:,:,nsm)      = hs(:,:,nsp)
    uflux(:,:,:,nsm) = uflux(:,:,:,nsp)
    vflux(:,:,:,nsm) = vflux(:,:,:,nsp)
!    wflux(:,nsm) = wflux(:,nsp) 
#if  zgrid3D
    dzt(:,:,:,nsm)   = dzt(:,:,:,nsp)
#endif
#if defined explicit_w || full_wflux
    wflux(:,:,:,nsm) = wflux(:,:,:,nsp)
#endif
#ifdef tempsalt
    tem(:,:,:,nsm)   = tem(:,:,:,nsp)
    sal(:,:,:,nsm)   = sal(:,:,:,nsp)
    rho(:,:,:,nsm)   = rho(:,:,:,nsp)
#endif
  end subroutine datasetswap


ENDMODULE mod_vel
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_name
  CHARACTER(LEN=200)                         :: inDataDir
  CHARACTER(LEN=200)                         :: topoDataDir
  CHARACTER(LEN=200)                         :: projDesc
  CHARACTER(LEN=200)                         :: GCMname   ,GCMsource
  CHARACTER(LEN=200)                         :: gridName  ,gridSource
  CHARACTER(LEN=200)                         :: gridDesc
  CHARACTER (LEN=23)                         :: Project, Case
  CHARACTER(LEN=200)                         :: caseName  ,caseDesc
ENDMODULE mod_name
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_streamfunctions
  USE mod_precdef
#ifdef streamxy
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)        :: stxyy, stxyx
#endif
#ifdef streamv
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:)        :: stxz, styz
#endif
#ifdef streamr
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)      :: stxr,styr, stzr
#endif
#ifdef stream_thermohaline
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:)      :: psi_ts
#endif
#ifdef tracer_convergence
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)      :: converg
#endif
  INTEGER                                    :: intpsi=120 
  ! to be read by the xxx.in files in future
#ifdef streamts
  INTEGER, PARAMETER                        :: LOV=3
#else
  INTEGER, PARAMETER                        :: LOV=1
#endif
ENDMODULE mod_streamfunctions
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_tracer
#ifdef tracer
  REAL, ALLOCATABLE, DIMENSION(:,:,:)        :: tra
#endif
ENDMODULE mod_tracer
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
#if defined diffusion || turb 
!#ifdef diffusion
MODULE mod_diffusion
  REAL                                       :: ah, av
ENDMODULE mod_diffusion
#endif
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===


! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
MODULE mod_sed
#ifdef sediment
  !  REAL :: wsed,rhos,D,critvel,T,cwamp,kincrit
  REAL                                       :: wsed, partdiam
  REAL                                       :: rhos, cwamp, twave
  REAL                                       :: critvel, kincrit

  INTEGER                                    :: nsed=0, nsusp=0
  LOGICAL                                    :: res
  INTEGER                                    :: numseedsubtimes=0
  REAL, ALLOCATABLE, DIMENSION(:)            :: seedsubtimes


#endif
ENDMODULE mod_sed
MODULE mod_orbital
#ifdef sediment
  REAL, ALLOCATABLE, DIMENSION(:)            :: orb
#endif
ENDMODULE mod_orbital
! ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===
