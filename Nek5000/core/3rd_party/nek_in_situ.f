c-----------------------------------------------------------------------
      subroutine in_situ_init()
#ifdef VISIT
      call visit_init()
#elif CATALYST
      call catalyst_init()
#endif
      end
c-----------------------------------------------------------------------
      subroutine in_situ_check()
#ifdef VISIT
      call visit_check()
#elif CATALYST
      call catalyst_process()
#endif
      end
c-----------------------------------------------------------------------
      subroutine in_situ_end()
#ifdef VISIT
      call visit_end()
#elif CATALYST
      call catalyst_end()
#endif
      end
c-----------------------------------------------------------------------

