
      integer     MAXNX,     MAXNY,     MAXN
      parameter ( MAXNX=128, MAXNY=128, MAXN=MAXNX*MAXNY )
 
      integer     MAX_FILL
      parameter ( MAX_FILL=64 )

      integer     MAX_SZP1
      parameter ( MAX_SZP1=MAXN+1 )

      integer     STENCIL_SZ
      parameter ( STENCIL_SZ=5 )

      integer     MAX_FACTR
      parameter ( MAX_FACTR=(2*MAX_FILL+STENCIL_SZ+1)*MAX_SZP1 )

      integer     LIPAR
      parameter ( LIPAR=3+STENCIL_SZ+4*MAX_SZP1+2*MAX_FACTR )

      integer     LRPAR
      parameter ( LRPAR=5+STENCIL_SZ+MAX_SZP1+2*MAX_FACTR )

