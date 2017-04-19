
module mod_write

  USE mod_name, only: case, caseNum, Project
  USE mod_time 
#ifdef ncwrite
  USE mod_ncwrite
#endif

  IMPLICIT NONE
  INTEGER                                    :: intminInOutFile
  CHARACTER(LEN=200)                         :: outDataDir, outDataFile
  CHARACTER (LEN=200)                        ::  projdir="", ormdir=""
  CHARACTER(LEN=200)                         :: inargstr1='', inargstr2=''
  INTEGER                                    :: twritetype = 0
  INTEGER                                    :: fileseq = 0
  CHARACTER(LEN=20)                          :: rankstamp=''
  LOGICAL                                    :: outdirdate = .true.
  LOGICAL                                    :: outdircase = .true.
  CHARACTER (LEN=10)                         :: yearstr
  
CONTAINS


  subroutine setup_outdatadir

    if (len(trim(outDataDir)) == 0) then
       outDataDir = 'projects/'
    end if
    if (outdircase .eqv. .true.) THEN
       outDataDir = trim(outDataDir) // trim(Project) // '/'
       outDataDir = trim(outDataDir) // 'results/' // trim(Case)
       yearstr = "xxxx"
       write (yearstr(:),'(I4.4)') caseNum
       outDataDir = trim(outDataDir) // trim(yearstr) // '/'
    END IF
    if (outdirdate .eqv. .true.) then
       yearstr = 'XXXXXXXXXX'
       write (yearstr(:),'(I10.10)') int(intmin)
       outDataDir = trim(outDataDir)//trim(yearstr) // '/'
    end if    
    call system('mkdir -p ' // trim(outDataDir))
    
  end subroutine setup_outdatadir

  
  subroutine open_outfiles

    IMPLICIT NONE
    CHARACTER(LEN=200)                         :: fullWritePref
    CHARACTER(LEN=20)                          :: intminstamp='', partstamp=''
    
    if ((intminInOutFile.eq.1) .or. (intminInOutFile.eq.3)) then
       write (intminstamp, '(A,i8.8)') '_t', intstart
    end if
    if ((intminInOutFile.eq.2) .or. (intminInOutFile.eq.3)) then
         write (partstamp, '(A,i6.6)') '_p', max(ints-intstart,0)+1
      end if
      
    fullWritePref =  trim(outDataDir)  // trim(outDataFile) //    &
                     trim(inargstr1)   // trim(inargstr2)   //    & 
                     trim(intminstamp) // trim(partstamp)   //    &
                     trim(rankstamp)

#if defined textwrite
    open(56, file=trim(fullWritePref)//'_run.asc')    
    open(57, file=trim(fullWritePref)//'_out.asc')  
    open(58, file=trim(fullWritePref)//'_ini.asc')   
    open(59, file=trim(fullWritePref)//'_err.asc')
#endif

#if defined binwrite
    open(unit=75 ,file=trim(fullWritePref)//'_out.bin', &  
         access='direct' ,form='unformatted' ,recl=20 ,status='replace')
    open(unit=76 ,file=trim(fullWritePref)//'_run.bin', &  
         access='direct' ,form='unformatted' ,recl=20 ,status='replace')
    open(unit=77 ,file=trim(fullWritePref)//'_kll.bin', &
         access='direct' ,form='unformatted' ,recl=20 ,status='replace')
    open(unit=78 ,file=trim(fullWritePref)//'_ini.bin', &  
         access='direct' ,form='unformatted' ,recl=20 ,status='replace')
    open(unit=79 ,file=trim(fullWritePref)//'_err.bin', &  
         access='direct' ,form='unformatted' ,recl=20 ,status='replace')
#endif

#if defined ncwrite
    CALL create_ncfile( trim(fullWritePref)//'_out.nc', ncid_out, NCtime_out)
    CALL create_ncfile( trim(fullWritePref)//'_run.nc', ncid_run, NCtime_run)
    CALL create_ncfile( trim(fullWritePref)//'_kll.nc', ncid_kll, NCtime_kll)
    CALL create_ncfile( trim(fullWritePref)//'_ini.nc', ncid_ini, NCtime_ini)
    CALL create_ncfile( trim(fullWritePref)//'_err.nc', ncid_err, NCtime_err)
#endif

#if defined csvwrite
    open(unit=85, file=trim(fullWritePref)//'_out.csv', status='replace')
    open(unit=86, file=trim(fullWritePref)//'_run.csv', status='replace')
    open(unit=87, file=trim(fullWritePref)//'_kll.csv', status='replace')
    open(unit=88, file=trim(fullWritePref)//'_ini.csv', status='replace')
    open(unit=89, file=trim(fullWritePref)//'_err.csv', status='replace')
#endif


#ifdef streamxy
    open(51,file=trim(fullWritePref)//'_psi_xy_yx.bin',form='unformatted')
#endif
#if defined streamv
    open(52,file=trim(fullWritePref)//'_psi_yz_xz.bin',form='unformatted')
#endif
#if defined streamr 
    open(53,file=trim(fullWritePref)//'_psi_xr_yr.bin',form='unformatted')
#endif
#ifdef stream_thermohaline
    open(54,file=trim(fullWritePref)//'_psi_ts.bin',form='unformatted')
#endif

#ifdef rerun
    open(67, file=trim(fullWritePref)//'_rerun.asc')
#endif




  end subroutine open_outfiles

  subroutine close_outfiles

#if defined textwrite
    close(56)
    close(57)
    close(58)
    close(59)
#endif
#if defined binwrite
    close(75)
    close(76)
    close(77)
    close(78)
    close(79)
#endif
#if defined ncwrite
    CALL check( nf90_close(ncid_out) )
    CALL check( nf90_close(ncid_run) )
    CALL check( nf90_close(ncid_kll) )
    CALL check( nf90_close(ncid_ini) )
    CALL check( nf90_close(ncid_err) )
#endif
#if defined csvwrite
    close(75)
    close(76)
    close(77)
    close(78)
    close(79)
#endif
  end subroutine close_outfiles

  subroutine writedata(sel)
    USE mod_time
    USE mod_pos
    USE mod_traj
    USE mod_loopvars
    USE mod_name

    IMPLICIT NONE

    REAL                                 :: vort
    INTEGER                              :: sel ,xf ,yf ,zf ,n
    INTEGER*8, SAVE                      :: recPosIn=0  ,recPosOut=0
    INTEGER*8, SAVE                      :: recPosRun=0 ,recPosErr=0
    INTEGER*8, SAVE                      :: recPosKll=0
    REAL                                 :: x14 ,y14 ,z14
    REAL*8                               :: twrite
    ! === Variables to interpolate fields ===
    REAL                                       :: temp, salt, dens
    REAL                                       :: temp2, salt2, dens2

566 format(2i8,3f10.4,2f12.4 &
         ,f12.0,f6.1,f6.2,f6.2,f6.0,8e8.1 )
    
    xf   = floor(x1)
    yf   = floor(y1)
    zf   = floor(z1)
    
    !if ((sel .ne. 19) .and. (sel.ne.40)) then
       ! this requires too much memory
       !       vort = (vvel(xf+1,yf,zf)-vvel(xf-1,yf,zf))/4000 - &
       !            (uvel(xf,yf+1,zf)-uvel(xf,yf-1,zf))/4000   
    !end if
    
subvol =  trj(5,ntrac)
t0     =  trj(7,ntrac)
#if defined tempsalt
    call interp2(ib,jb,kb,temp,salt,dens)
#endif

#if defined textwrite 
    select case (sel)
    case (10)
       write(58,566) ntrac,niter,x1,y1,z1,tt/tday,t0/tday,subvol,temp,salt,dens
    case (11)
       if(  (kriva == 1 .AND. nrj(4,ntrac) == niter-1   ) .or. &
            (kriva == 2 .AND. scrivi                    ) .or. &
            (kriva == 3                                 ) .or. &
            (kriva == 4 .AND. niter == 1                ) .or. &
            (kriva == 5 .AND.                                  &
          &  MOD((REAL(tt)-REAL(t0))*REAL(NGCM)/REAL(ITER), 3600.) == 0.d0 ) .or. &
            (kriva == 6 .AND. .not.scrivi                  ) ) then
#if defined biol
          write(56,566) ntrac,ints,x1,y1,z1,tt/3600.,t0/3600.
#else
#if defined tempsalt
          write(56,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol,temp,salt,dens
#else
          write(56,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol
#endif        
#endif        
       endif
    case (13)
       ! === write sed pos ===
       write(57,566) ntrac,niter,x1,y1,z1, &
            tt/tday,t0/tday,subvol,temp,salt,dens 
    case (14)
       write(56,566) ntrac,ints,x1,y1,z1, &
            tt/60.,t0/3600.,subvol,temp,salt,dens
    case (15)
       write(57,566) ntrac,ints,x1,y1,z1, &
            tt/tday,t0/tday,subvol,temp,salt,dens
    case (16)
       if(kriva.ne.0 ) then
#if defined tempsalt
           !call interp(ib,jb,kb,x1,y1,z1,temp,salt,dens,1) 
           call interp2(ib,jb,kb,temp,salt,dens)
#endif
          write(56,566) ntrac,ints,x1,y1,z1, &
               tt/tday,t0/tday,subvol,temp,salt,dens
       end if
    case (17)
       write(57,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens  
    case (19)
       ! === write last sedimentation positions ===
       open(34,file=trim(outDataDir)//trim(outDataFile)//'_sed.asc') 
       do n=1,ntracmax
        if(nrj(1,n).ne.0) then
         write(34,566) n,nrj(4,n),trj(1,n),trj(2,n),trj(3,n),trj(4,n)/tday,trj(7,n)/tday
      endif
       enddo
       close(34)
    case (40)
       write(59,566) ntrac,ints,x1,y1,z1,tt/tday,t0/tday,subvol &
            ,temp,salt,dens  
    case (99) !switch
       
    end select
#endif 
   
#if defined binwrite 
    x14=real(x1,kind=4)
    y14=real(y1,kind=4)
    z14=real(z1,kind=4)
    twrite = tt/tday
    select case (sel)       
    case (10) !in
       recPosIn = recPosIn + 1
       write(unit=78 ,rec=recPosIn)  real(ntrac,kind=4),real(twrite,kind=4),x14,y14,z14
       return
    case (11)
       if(  (kriva == 1 .and. nrj(4,ntrac)  ==  niter-1 ) .or. &
            (kriva == 2 .and. scrivi                    ) .or. &
            (kriva == 3                                 ) .or. &
            (kriva == 4 .and. niter == 1                ) .or. &
            (kriva == 5 .and. abs(dmod(tt-t0,9.d0)) < 1e-5 ) .or. &
            (kriva == 6 .and. .not.scrivi                  ) ) then
#if defined tempsalt
          call interp(ib,jb,kb,x1,y1,z1,temp, salt,  dens,1)
 !         call interp(ib,jb,kb,x1,y1,z1,temp2,salt2, dens2,2)
          !z14=real(salt*rb+salt2*(1-rb),kind=4)
#endif
          recPosRun = recPosRun+1
          write(unit=76 ,rec=recPosRun) real(ntrac,kind=4),real(twrite,kind=4),x14,y14,z14
       end if
    case (13)
       recPosKll = recPosKll + 1
       write(unit=77 ,rec=recPosKll) real(ntrac,kind=4),real(twrite,kind=4),x14,y14,z14
    case (17) !out
       recPosOut = recPosOut + 1
       write(unit=77 ,rec=recPosOut) real(ntrac,kind=4),real(twrite,kind=4),x14,y14,z14
    case (19) !end
       recPosOut = recPosOut + 1
       write(unit=75 ,rec=recPosOut) real(ntrac,kind=4),real(twrite,kind=4),x14,y14,z14
    case (40) !error
       recPosErr=recPosErr + 1    
       write(unit=79 ,rec=recPosErr) real(ntrac,kind=4),real(twrite,kind=4),x14,y14,z14
    case (99) !switch
       if ((recPosRun > 50000000).and.(intminInOutFile.eq.2)) then
          call close_outfiles
          call open_outfiles
          recPosRun = 0
          recPosIn  = 0
          recPosOut = 0
          recPosErr = 0
          print *, "Switched run file" 
       end if
    end select
!CA we write either in binary or netcdf, not both
#else ifdef ncwrite 
    SELECT CASE (sel)       
    CASE (10)
       CALL write_ncpos(ncid_ini, ntrac, x1, y1, z1, tt, NCtime_ini)
       return
    CASE (11)
       if(  (kriva == 1 .and. nrj(4,ntrac)  ==  niter-1 ) .or. &
            (kriva == 2 .and. scrivi                    ) .or. &
            (kriva == 3                                 ) .or. &
            (kriva == 4 .and. niter == 1                ) .or. &
            (kriva == 5 .and. abs(dmod(tt-t0,9.d0)) < 1e-5 ) .or. &
            (kriva == 6 .and. .not.scrivi                  ) ) then
#if defined tempsalt
          STOP "UNIMPLEMENTED"
#endif
          CALL write_ncpos(ncid_run, ntrac, x1, y1, z1, tt, NCtime_run)
       end if
    CASE (13)
       CALL write_ncpos(ncid_kll, ntrac, x1, y1, z1, tt, NCtime_kll)
    CASE (17) !out
       CALL write_ncpos(ncid_kll, ntrac, x1, y1, z1, tt, NCtime_kll)
    CASE (19) !end
       CALL write_ncpos(ncid_out, ntrac, x1, y1, z1, tt, NCtime_out)
    CASE (40) !error
       CALL write_ncpos(ncid_err, ntrac, x1, y1, z1, tt, NCtime_err)
    END SELECT
#endif /* ncwrite */

#if defined csvwrite 
    x14=real(x1,kind=4)
    y14=real(y1,kind=4)
    z14=real(z1,kind=4)
    if (twritetype==1) then
       twrite = tt
    else if (twritetype==2) then
       call updateclock
       twrite = currJDtot
    else
       twrite = real(ints,kind=8)
    end if
    select case (sel)       
    case (10)
       write(88,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
       return
    case (11)
       if(  (kriva == 1 .and. nrj(4,ntrac)  ==  niter-1 ) .or. &
            (kriva == 2 .and. scrivi                    ) .or. &
            (kriva == 3                                 ) .or. &
            (kriva == 4 .and. niter == 1                ) .or. &
            (kriva == 5 .and. abs(dmod(tt-t0,9.d0)) < 1e-5 ) .or. &
            (kriva == 6 .and. .not.scrivi                  )  ) then
          !!!! CALL FIELD-INTERP !!!!
          write(86,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
       end if
    case (13)
       write(87,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (15)
       write(86,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (17)
       write(87,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (19)
       write(85,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (40)
       write(89,"(I0,4(',',F0.5))")  ntrac, twrite, x14, y14, z14
    case (99) !switch
       
    end select
#endif   
  end subroutine writedata


end module mod_write
