!=======================================================================
!     Adam Peplinski; 2015.10.20
!     Set of subroutines to read user and module parameters using
!     namelists.
!
!=======================================================================
!***********************************************************************
!     read parameters from the file
      subroutine uprm_read
      implicit none

      include 'SIZE'            ! NID
      include 'INPUT'           ! REAFLE

!     local variables
      integer len, ierr
      integer iunit

      character*132 fname

!     functions
      integer ltrunc
!-----------------------------------------------------------------------
!     Open parameter file and read contents
      ierr=0
      if (NID.eq.0) then
!     find free unit
         call IO_freeid(iunit, ierr)
!     get file name and open file
         if(ierr.eq.0) then
            call blank(fname,132)
            len = ltrunc(REAFLE,132) - 4
            call chcopy(fname,REAFLE,len)
            fname(len+1:len+6)='.upar'
            write(6,*) 'Openning parameter file: ',trim(fname)
            open (unit=iunit,file=fname,status='old',action='read',
     $           iostat=ierr)
         endif
      endif
      call err_chk(ierr,'Error opening .upar file.$')

!     place to call module _param_in routines
      call uprm_in(iunit)

!     close the file
      ierr=0
      if (NID.eq.0) close(unit=iunit,iostat=ierr)
      call err_chk(ierr,'Error closing .upar file.$')

!     stamp logs
      if (NIO.eq.0) write(*,*) 'User parameter list'
      call uprm_out(6)
      if (NIO.eq.0) write(*,*)

      return
      end
!***********************************************************************
!     read parameters
      subroutine uprm_in(iunit)
      implicit none

!     argument list
      integer iunit
!-----------------------------------------------------------------------
!     place to call module _param_in routines

      rewind(iunit)
!     user parameters
      call user_param_in(iunit)

      rewind(iunit)
!     restart
      call chkpt_param_in(iunit)

!      rewind(iunit)
!     sponge
!      call spng_param_in(iunit)

!      rewind(iunit)
!     SFD
!      call sfd_param_in(iunit)

      rewind(iunit)
!     RTFILTER
      call rtfil_param_in(iunit)

      rewind(iunit)
!     STATS

      return
      end
!***********************************************************************
!     output parameters
      subroutine uprm_out(iunit)
      implicit none

!     argument list
      integer iunit
!-----------------------------------------------------------------------
!     place to call module _param_in routines

!     user parameters
      call user_param_out(iunit)

!     restart
      call chkpt_param_out(iunit)

!     sponge
!      call spng_param_out(iunit)

!     SFD
!      call sfd_param_out(iunit)

!     RTFILTER
      call rtfil_param_out(iunit)

!     STATS

      return
      end
!***********************************************************************
!     read user parameters
      subroutine user_param_in(fid)
      implicit none

      include 'SIZE'            !
      include 'PARALLEL'        ! ISIZE, WDSIZE, LSIZE,CSIZE
      include 'USERPAR'         !

!     argument list
      integer fid

!     local variables
      integer ierr

!     namelists; cannot be empty
      namelist /USERPAR/ UPRM_PRB,L2FREQ,FIXGEOM,NEW_DT
!-----------------------------------------------------------------------
!     default values
! ! !       UPRM_PRB           = 10
! ! !       L2FREQ             = 1000
! ! !       FIXGEOM            = 0
! ! !       NEW_DT             = 1.0E-6
!     read the file
      ierr=0
! ! !       if (NID.eq.0) then
! ! !          read(unit=fid,nml=USERPAR,iostat=ierr)
!     add exceptions in the future and check
!     if it is parameter and compiler independent
!     iostat
!     84 - NAMELIST group header not found in external file
!     85 - NAMELIST group header not found in internal file

      return
      end
!***********************************************************************
!     write user parameters
      subroutine user_param_out(fid)
      implicit none

      include 'SIZE'            !
      include 'USERPAR'         !

!     argument list
      integer fid               ! file id

!     local variables
      integer ierr

!     namelists; cannot be empty
      namelist /USERPAR/ UPRM_PRB
!-----------------------------------------------------------------------
      ierr=0
      if (NID.eq.0) then
         write(unit=fid,nml=USERPAR,iostat=ierr)
      endif
      call err_chk(ierr,'Error writing USERPAR parameters.$')

      return
      end
!***********************************************************************

!
      subroutine uprm_read_MA

      implicit none

      include 'SIZE'
      include 'GEOM'                    ! xm1, ym1, zm1
      include 'SOLN'                    ! T
      include 'MASS'                    !BM1 for lambda2
      include 'TSTEP'                   ! ISTEP
      include 'INPUT'                   ! PARAM(12) (DT)
      include 'USERPAR'                 ! l2freq, FIXGEOM, NEW_DT
      include 'RTFILTER'                ! Diagnostic spectra only. Can be removed later.
      include 'STATS'
      include 'CHKPOINT'

!-----------------------------------------------------------------------

      UPRM_PRB = INT(PARAM(69))
      L2FREQ = INT(PARAM(70))

      NEW_DT = PARAM(72)

!
      CHKPTSTEP = INT(PARAM(75))
!
      IFCHKPTRST = .FALSE.
      if ((INT(PARAM(76))).gt.0) then
         IFCHKPTRST = .TRUE.
      endif
!
      rt_kut = INT(PARAM(78))
      rt_kai = PARAM(79)
      rt_wght = PARAM(80)
!
      rt_ifboyd     =.TRUE.
      if ((INT(PARAM(81)).le.0)) then
            rt_ifboyd = .false.
      endif
!
      stat_comp = INT(PARAM(87))
      stat_outp = INT(PARAM(88))


      return
      end
!***********************************************************************
