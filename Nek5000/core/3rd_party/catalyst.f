      subroutine catalyst_init()
#ifdef CATALYST
      include "mpif.h"
      CHARACTER*4 core
      integer par_rank, err
      common /PARALLELS/ par_rank, core
      save /PARALLELS/
c      implicit none
      call MPI_COMM_RANK(MPI_COMM_WORLD, par_rank, err)
      write(core, 10) par_rank
 10   format (I4)
      core = adjustl(trim(core))
      open(unit=44, file='perf/'//'core_'//adjustl(trim(core))//'.csv')
#endif
      call coprocessorinitialize()
      ! Add user defined pipelines
      call catalyst_usrpipe()

      end

      subroutine catalyst_end
#ifdef CATALYST
      close(unit=44)
#endif
      call coprocessorfinalize()

      end

      subroutine catalyst_process()
      implicit none
      include "mpif.h"
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
#ifdef CATALYST
      double precision before, after, lafter, cat, sim
      common /PARALLELS/ before, after, cat , sim
      save /PARALLELS/
c
      integer flag, dim

      before = MPI_Wtime()
#endif
c      lafter=0.0
      call requestdatadescription(istep, time, flag)
      if (flag .ne. 0) then
         call needtocreategrid(flag)
         dim = 2
         if (IF3D) dim = 3
         call creategrid(xm1, ym1, zm1, lx1, ly1, lz1, lelt, dim)
         call add_scalar_field(pr, "pressure"//char(0))
         call add_vector_field(vx, vy, vz, dim, "velocity"//char(0))
         call add_scalar_field(t, "temperature"//char(0))
         call coprocess()
      end if
#ifdef CATALYST
      lafter = after
      after = MPI_WTIME()

      cat = after-before
      sim = before-lafter

c      lafter = after

      WRITE (44, *)  before, ',', after, ',', cat, ',', sim
#endif

      end
