      subroutine catalyst_init
      
      call coprocessorinitialize()

      ! Add user defined pipelines
      call catalyst_usrpipe()

      end

      subroutine catalyst_end
      
      call coprocessorfinalize()

      end

      subroutine catalyst_process()
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'

      integer flag, dim

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
      
      end
         

