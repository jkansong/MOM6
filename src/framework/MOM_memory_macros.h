!***********************************************************************
! This is a header file to define macros for static and dynamic memory *
! allocation.  Define STATIC_MEMORY_ in MOM_memory.h for static memory *
! allocation.  Otherwise dynamic memory allocation will be assumed.    *
!***********************************************************************

#ifdef STATIC_MEMORY_
#  define DEALLOC_(x)
#  define ALLOC_(x)
#  define ALLOCABLE_
#  define PTR_
#  define TO_NULL_

!   NIMEM and NJMEM are the maximum number of grid points in the
! x- and y-directions on each processsor.
#  define NIMEM_ (((NIGLOBAL_-1)/NIPROC_)+1+2*NIHALO_)
#  define NJMEM_ (((NJGLOBAL_-1)/NJPROC_)+1+2*NJHALO_)

#  ifdef SYMMETRIC_MEMORY_
#    define NIMEMB_   0:NIMEM_
#    define NJMEMB_   0:NJMEM_
#  else
#    define NIMEMB_   NIMEM_
#    define NJMEMB_   NJMEM_
#  endif
#  define NIMEMB_PTR_ NIMEMB_
#  define NJMEMB_PTR_ NJMEMB_
#  define NIMEMB_SYM_ 0:NIMEM_
#  define NJMEMB_SYM_ 0:NJMEM_
#  define NKMEM_      NK_
#  define NKMEM0_     0:NK_
#  define NK_INTERFACE_ NK_+1
#  define C1_         1
#  define C2_         2
#  define C3_         3
#  define SZI_(G)     NIMEM_
#  define SZJ_(G)     NJMEM_
#  define SZK_(G)     NK_
#  define SZK0_(G)    0:NK_
#  define SZIB_(G)    NIMEMB_
#  define SZJB_(G)    NJMEMB_
#  define SZIBS_(G)   0:NIMEM_
#  define SZJBS_(G)   0:NJMEM_
! These macros are for block structure.
#  ifndef NIBLOCK_
#    define NIBLOCK_ 1
#  endif
#  ifndef NJBLOCK_
#    define NJBLOCK_ 1
#  endif
#  define NIMEM_BK_ ((((NIGLOBAL_-1)/NIPROC_)/NIBLOCK_)+1+2*NIHALO_)
#  define NJMEM_BK_ ((((NJGLOBAL_-1)/NJPROC_)/NJBLOCK_)+1+2*NJHALO_)
#  ifdef SYMMETRIC_MEMORY_
#    define NIMEMB_BK_   0:NIMEM_BK_
#    define NJMEMB_BK_   0:NJMEM_BK_
#  else
#    define NIMEMB_BK_   NIMEM_BK_
#    define NJMEMB_BK_   NJMEM_BK_
#  endif
#  define SZI_BK_(G)     NIMEM_BK_
#  define SZJ_BK_(G)     NJMEM_BK_
#  define SZIB_BK_(G)    NIMEMB_BK_
#  define SZJB_BK_(G)    NJMEMB_BK_

#else
! Dynamic memory allocation

#  define DEALLOC_(x) deallocate(x)
#  define ALLOC_(x)   allocate(x)
#  define ALLOCABLE_ ,allocatable
#  define PTR_       ,pointer
#  define TO_NULL_   =>NULL()
#  define NIMEM_      :
#  define NJMEM_      :
#  define NIMEMB_PTR_ :
#  define NJMEMB_PTR_ :
#  ifdef SYMMETRIC_MEMORY_
#    define NIMEMB_   0:
#    define NJMEMB_   0:
#  else
#    define NIMEMB_   :
#    define NJMEMB_   :
#  endif
#  define NIMEMB_SYM_ 0:
#  define NJMEMB_SYM_ 0:
#  define NKMEM_      :
#  define NKMEM0_     0:
#  define NK_INTERFACE_  :
#  define C1_         :
#  define C2_         :
#  define C3_         :
#  define SZI_(G)     G%isd:G%ied
#  define SZJ_(G)     G%jsd:G%jed
#  define SZK_(G)     G%ks:G%ke
#  define SZK0_(G)    G%ks-1:G%ke
#  define SZIB_(G)    G%IsdB:G%IedB
#  define SZJB_(G)    G%JsdB:G%JedB
#  define SZIBS_(G)   G%isd-1:G%ied
#  define SZJBS_(G)   G%jsd-1:G%jed

! These macros are for block structure.
#  define NIMEM_BK_       :
#  define NJMEM_BK_       :
#  ifdef SYMMETRIC_MEMORY_
#    define NIMEMB_BK_   0:
#    define NJMEMB_BK_   0:
#  else
#    define NIMEMB_BK_   :
#    define NJMEMB_BK_  :
#  endif
#  define SZI_BK_(G)     G%isd_bk:G%ied_bk
#  define SZJ_BK_(G)     G%jsd_bk:G%jed_bk
#  define SZIB_BK_(G)    G%isdB_bk:G%iedB_bk
#  define SZJB_BK_(G)    G%jsdB_bk:G%jedB_bk
#endif
