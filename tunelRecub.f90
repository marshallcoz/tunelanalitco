      module datos
      save
      real*8, dimension(2), parameter :: rho = (/1.45, 1.0/)
      real*8, dimension(2), parameter :: nu = (/0.25, 0.20/)
      real*8, dimension(2), parameter :: bet = (/2.405, 1.700/) * 10
      real*8, dimension(2) :: alf
      integer,parameter :: nRes = 2
      real*8, dimension(2+nRes) :: radios ! a, b 
      real, parameter :: Qq = 1000.0, periodicdamper = 1.0, TW = 1.0
      real*8, parameter :: DFREC = 0.05
      integer, parameter :: NFREC = 100, nplot = 100
      integer, parameter :: nmax = 500, nfracs = 1
      real*8, parameter :: Dns = 1.
      complex*16, parameter :: UI = cmplx(0.0d0,1.0d0,8), &
                               UR = cmplx(1.0d0,0.0d0,8), &
                               Z0 = cmplx(0.0d0,0.0d0,8)
      real*8, parameter :: PI = real(4.0d0*ATAN(1.0d0),8)
      contains
      subroutine set_radios
      radios(1) = 3.00 ! a
      radios(2) = radios(1)*1.1 ! b
      end subroutine set_radios
      end module datos
      module vars_func_of_w
      save
      complex*16 :: cOME
      real*8 :: FREC,OME,OMEI
      complex*16, dimension(2) :: lambda,amu
      
      complex*16, dimension(2,2),target :: VEL !( velocidad, region)
      !    cp(1)  cp(2) 
      !    cs(2)  cs(2)
      complex*16, dimension(:), pointer :: cp,cs 
      
      ! compressional and shear wave numbers:
      complex*16, dimension(2,2),target :: w_c  ! w/c ( velocidad, region)
      !    w/cp(1)  w/cp(2) 
      !    w/cs(2)  w/cs(2)
      complex*16, dimension(:), pointer :: alfa,beta 
      
      !ratios of the medium and the liner
      real*8 :: eta
      complex*16 :: muR,gammaP,gammaS
      
      !wavenumbers (squared)
      complex*16, dimension(2) :: s2,p2
      end module vars_func_of_w
      module RES
      use datos, only : NFREC,nRes
      type receptor
      real*8 :: r,th ! coordenadas
      integer :: reg ! 1:region de afuera
                     ! 2:inclusion
      complex*16,dimension(NFREC) :: & 
        s_rr_1,s_tt_1,s_rt_1,u_r_1,u_t_1, &
        s_rr_2,s_tt_2,s_rt_2,u_r_2,u_t_2
      end type receptor
      type(receptor), dimension(nRes),target :: Rw
      complex*16,dimension(NFREC) :: abscisa
      contains
      subroutine setup_resu
      use datos, only : pi,radios,z0
      integer :: i
        Rw(1)%r = radios(2) !r=b in the medium at the boundary of the liner
        Rw(1)%th = pi!/2.
        Rw(1)%reg = 1 ! in the medium
       
        Rw(2)%r = radios(1) !r=a at the inner side of the cylinder
        Rw(2)%th = pi!/2.
        Rw(2)%reg = 2 ! in the liner
        do i = 1,nRes
         radios(2+i) = Rw(i)%r
        ! reg = 1
          Rw(i)%s_rr_1 = z0
          Rw(i)%s_tt_1 = z0
          Rw(i)%s_rt_1 = z0
          Rw(i)%u_r_1 = z0
          Rw(i)%u_t_1 = z0
         ! reg = 2
          Rw(i)%s_rr_2 = z0
          Rw(i)%s_tt_2 = z0
          Rw(i)%s_rt_2 = z0
          Rw(i)%u_r_2 = z0
          Rw(i)%u_t_2 = z0
        end do
      end subroutine setup_resu
      
      subroutine setup_resu_sabana
      use datos, only : pi
        integer :: i
        real*8,parameter :: r_ini = 0.1
        real*8,parameter :: r_delta = 0.2
        real*8,parameter :: th_ini = pi/2.
        real*8,parameter :: th_delta = 0.0
        do i = 1,nRes
          Rw(i)%r = r_ini + (i-1)* r_delta
          Rw(i)%th = th_ini + (i-1)* th_delta
        end do
      end subroutine setup_resu_sabana
      end module RES
      module debug
      contains
      subroutine showMNmatrixZ(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n,outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(A,E9.2,A)',advance='no') "(",REAL(MAT(i,j)),","
          write(outpf,'(E9.2,A)',advance='no') AIMAG(MAT(i,j)),"i) "
        end do
        write(outpf,'(A)',advance='yes')''
      end do
      end subroutine
      end module debug
      module plotter 
      contains
      subroutine plot(i_Px,i_nom,reg, abscisa ,pt_RES,nplot)
      use dislin
      implicit none
      integer :: i_Px,i_nom,reg,nplot
      character(len=4),dimension(5) :: nom
      complex*16,dimension(nplot) :: abscisa,pt_RES
      character(len=100) :: tx
      nom(1) = "s_rr"
      nom(2) = "s_tt"
      nom(3) = "s_rt"
      nom(4) = "u_r_"
      nom(5) = "u_t_"
      ! graficar
      CALL METAFL('PDF')
      call filmod('DELETE')
      CALL SETPAG('DA4P')
      write(tx,'(a,i0,a,i0,a)') "RES_",i_Px,nom(i_nom),reg,".pdf"
      CALL SETFIL(trim(adjustl(tx)))
      
      call qplcrv(real(abscisa,4), & 
                 real(real(pt_RES),4), nplot,'FIRST')
      call color('RED')
      call qplcrv(real(abscisa,4), & 
                 real(real(pt_RES),4), nplot,'NEXT')
      call color('BLUE')
      call qplcrv(real(abscisa,4), & 
                 real(aimag(pt_RES),4), nplot,'LAST')
      
      write(tx,'(a,i0,a,i0,a)') "RES_",i_Px,nom(i_nom),reg,"abs.pdf"
      CALL SETFIL(trim(adjustl(tx)))
      CALL color('BLACK')
      call qplot(real(abscisa,4), & 
                 real(abs(pt_RES),4), nplot)
      end subroutine plot
      end module plotter
      module Hank
      use datos, only : nmax, nfracs
      use RES, only : nRes
      type argvelocidad
        complex*16, dimension (2) :: c ! 1 para alfa; 2 para beta
      end type argvelocidad
      
      type argradiovelocidad
        type (argvelocidad), dimension(2+nRes) :: r ! 1 para a;  2 para b; 3... para receptores
      end type argradiovelocidad
      
      type Bessel
        type(argradiovelocidad),dimension(1:4,-1:nmax*nfracs+nfracs) :: JYH1H2 ! -1,0,1,2,...,nmax
      ! i = 1 :: Jn
      ! i = 2 :: Yn
      ! i = 3 :: H(1)n
      ! i = 4 :: H(2)n
      end type Bessel
      
      type(Bessel), dimension(2),target :: Bess ! 1 para afuera; 2 para recubrimiento
      
      contains
      subroutine makeBessels(imprimo)
      USE SPECFUN
      use vars_func_of_w, only : w_c !( velocidad, region)
      use datos, only : nmax,nfracs,Dns,radios,ui,z0
      use RES, only : nRes,Rw
      use dislin
      use debug
      use plotter
      implicit none
      logical :: imprimo
      integer :: reg,c,r,imec,i,ns
      complex*16 :: z!,H20,H21,H22
      real*8 :: v,vm
      complex*16, dimension(0:nmax) :: CBJ,CDJ,CBY,CDY
      character(Len=3), dimension(4) :: nom
      character(len=100) :: tx
      real*4,dimension(nmax) :: x, pt_RES 
      nom(1) = "__J"
      nom(2) = "__Y"
      nom(3) = "_H1"
      nom(4) = "_H2"
      
      do reg=1,2  ! afuera 1 y adentro 2
      do r=1,2+nRes !a,b
      if (r .ge. 3) then
        if (Rw(r-2)%reg .ne. reg) cycle
      end if
      do c=1,2 !alfa y beta 
      Z = w_c(c,reg) * radios(r) 
      do i = 1,nfracs
      V = Dns * real(i-1,8) + real(nmax,8)
!     print*,V,Z
      Bess(reg)%JYH1H2(1,i-1:(nmax*nfracs)+(i-1):nfracs)%r(r)%c(c) = z0
      call CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
      ns = int(vm)
      ! J
      Bess(reg)%JYH1H2(1,i-1:(ns*nfracs)+(i-1):nfracs)%r(r)%c(c) = CBJ(0:ns)
      ! Y
      Bess(reg)%JYH1H2(2,i-1:(ns*nfracs)+(i-1):nfracs)%r(r)%c(c) = CBY(0:ns)
      ! H1
      Bess(reg)%JYH1H2(3,i-1:(ns*nfracs)+(i-1):nfracs)%r(r)%c(c) = CBJ(0:ns) + UI*CBY(0:ns)
      ! H2
      Bess(reg)%JYH1H2(4,i-1:(ns*nfracs)+(i-1):nfracs)%r(r)%c(c) = CBJ(0:ns) - UI*CBY(0:ns)
      end do ! i
      Bess(reg)%JYH1H2(1:4,-1)%r(r)%c(c) = - Bess(reg)%JYH1H2(1:4,1)%r(r)%c(c)
      
!     do i = 0,nmax* nfracs
!     print*,i,Bess(reg)%JYH1H2(3,i)%r(r)%c(c)
!     end do ! i
!     stop 199
      
!     ! Hankel of the second kind 
!     ! n = 0,1,2
!     z = w_c(c,reg) * radios(r) ! w_c(c,reg) * radios(r)
!     call hankels(z,H20,H21)
!     H22 = -H20 + 2./z * H21
!     Bess(reg)%JYH1H2(4,-1)%r(r)%c(c) = -H21 !-1
!     Bess(reg)%JYH1H2(4,0)%r(r)%c(c) = H20   !-0
!     Bess(reg)%JYH1H2(4,1)%r(r)%c(c) = H21   !+1
!     Bess(reg)%JYH1H2(4,2)%r(r)%c(c) = H22   !+2
!     ! n = 3,4,...nmax
!     do n=2,nmax
!     Bess(reg)%JYH1H2(4,n+1)%r(r)%c(c) = & 
!     2* n / z * Bess(reg)%JYH1H2(4,n)%r(r)%c(c) - Bess(reg)%JYH1H2(4,n-1)%r(r)%c(c)
!     end do
!     ! Hankel of the first kind
!     ! n = 0,1,2
!     ! tomando de la continuaci´on anal´itica modificamos el argumento
!     z = w_c(c,reg) * radios(r) * exp(-pi*UI) 
!     call hankels(z,H20,H21)
!     H22 = -H20 + 2./z * H21
!     Bess(reg)%JYH1H2(3,-1)%r(r)%c(c) = -H21
!     Bess(reg)%JYH1H2(3,0)%r(r)%c(c) = H20
!     Bess(reg)%JYH1H2(3,1)%r(r)%c(c) = H21
!     Bess(reg)%JYH1H2(3,2)%r(r)%c(c) = H22
!     do n=2,nmax
!     Bess(reg)%JYH1H2(3,n+1)%r(r)%c(c) = & 
!     2 * n / z * Bess(reg)%JYH1H2(3,n)%r(r)%c(c) - Bess(reg)%JYH1H2(3,n-1)%r(r)%c(c)
!     end do
!     ! continuación analítica
!     do n=-1,nmax+1
!     Bess(reg)%JYH1H2(3,n)%r(r)%c(c) = exp(-real(n)*pi*UI) * Bess(reg)%JYH1H2(3,n)%r(r)%c(c) ! (-) ?
!     end do
!     
!     z = w_c(c,reg) * radios(r) ! alfa(c) * radios(r)
!     call BESJC(z,nmax,Bess(reg)%JYH1H2(1,0:nmax+1)%r(r)%c(c))
!     Bess(reg)%JYH1H2(1,-1)%r(r)%c(c) = - Bess(reg)%JYH1H2(1,1)%r(r)%c(c) 
!     z = w_c(c,reg) * radios(r) ! alfa(c) * radios(r)
!     call BESYC(z,nmax,Bess(reg)%JYH1H2(2,0:nmax+1)%r(r)%c(c))
!     Bess(reg)%JYH1H2(2,-1)%r(r)%c(c) = - Bess(reg)%JYH1H2(2,1)%r(r)%c(c) 
!     
      if (imprimo .eqv. .true.) then
      do imec = 1,4
        write(tx,'(a,a,i0,a,i0,a,i0,a)') nom(imec),"_",reg,"_",r,"_",c,".pdf"
        print*,trim(tx)," | z=",radios(r) * w_c(c,reg)
        call showMNmatrixZ(NMAX,1,Bess(reg)%JYH1H2(imec,0:NMAX)%r(r)%c(c),"-----",6)
          
          x = (/(real(i*dns,4),i=0,nmax-1)/)
      CALL METAFL('PDF')
      call filmod('DELETE')
      CALL SETPAG('DA4P')
      CALL SETFIL(trim(adjustl(tx)))
      
          pt_RES = real(Bess(reg)%JYH1H2(imec,0:NMAX-1)%r(r)%c(c),4)
      call qplcrv(real(x,4),pt_RES, nmax,'FIRST')
      call color('RED')
      call qplcrv(real(x,4),pt_RES, nmax,'NEXT')
      call color('BLUE')
      pt_RES = real(aimag (Bess(reg)%JYH1H2(imec,0:NMAX-1)%r(r)%c(c)),4)
      call qplcrv(real(x,4),pt_RES, nmax,'LAST')
      
          
      end do
      end if
      end do !c
      end do !r
      end do !reg
      end subroutine makeBessels
      end module hank

        module specfun
        contains
        SUBROUTINE CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
!
!       ===========================================================
!       Purpose: Compute Bessel functions Jv(z), Yv(z) and their 
!                derivatives for a complex argument
!       Input :  z --- Complex argument
!                v --- Order of Jv(z) and Yv(z)
!                      ( v = n+v0, n = 0,1,2,..., 0 Û v0 < 1 )
!       Output:  CBJ(n) --- Jn+v0(z)
!                CDJ(n) --- Jn+v0'(z)
!                CBY(n) --- Yn+v0(z)
!                CDY(n) --- Yn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ===========================================================
!
        IMPLICIT DOUBLE PRECISION (A,B,G,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:*),CDJ(0:*),CBY(0:*),CDY(0:*)
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z1=Z
        Z2=Z*Z
        N=INT(V)
        V0=V-N
        PV0=PI*V0
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           IF (V0.EQ.0.0) THEN
              CBJ(0)=(1.0D0,0.0D0)
              CDJ(1)=(0.5D0,0.0D0)
           ELSE
              CDJ(0)=(1.0D+300,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0D0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CJV0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 15 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+V0))
              CJV0=CJV0+CR
              IF (CDABS(CR).LT.CDABS(CJV0)*1.0D-15) GO TO 20
15         CONTINUE
20         VG=1.0D0+V0
           CALL GAMMA(VG,GA)
           CA=(0.5D0*Z1)**V0/GA
           CJV0=CJV0*CA
        ELSE
           K0=11
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           VV=4.0D0*V0*V0
           CPZ=(1.0D0,0.0D0)
           CRP=(1.0D0,0.0D0)
           DO 25 K=1,K0
              CRP=-0.78125D-2*CRP*(VV-(4.0*K-3.0)**2.0)*(VV-&
                          (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*Z2)
25            CPZ=CPZ+CRP
           CQZ=(1.0D0,0.0D0)
           CRQ=(1.0D0,0.0D0)
           DO 30 K=1,K0
              CRQ=-0.78125D-2*CRQ*(VV-(4.0*K-1.0)**2.0)*(VV-&
                          (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*Z2)
30            CQZ=CQZ+CRQ
           CQZ=0.125D0*(VV-1.0)*CQZ/Z1
           ZK=Z1-(0.5D0*V0+0.25D0)*PI
           CA0=CDSQRT(RP2/Z1)
           CCK=CDCOS(ZK)
           CSK=CDSIN(ZK)
           CJV0=CA0*(CPZ*CCK-CQZ*CSK)
           CYV0=CA0*(CPZ*CSK+CQZ*CCK)
        ENDIF
        IF (A0.LE.12.0) THEN
           IF (V0.NE.0.0) THEN
              CJVN=(1.0D0,0.0D0)
              CR=(1.0D0,0.0D0)
              DO 35 K=1,40
                 CR=-0.25D0*CR*Z2/(K*(K-V0))
                 CJVN=CJVN+CR
                 IF (CDABS(CR).LT.CDABS(CJVN)*1.0D-15) GO TO 40
35            CONTINUE
40            VG=1.0D0-V0
              CALL GAMMA(VG,GB)
              CB=(2.0D0/Z1)**V0/GB
              CJU0=CJVN*CB
              CYV0=(CJV0*DCOS(PV0)-CJU0)/DSIN(PV0)
           ELSE
              CEC=CDLOG(Z1/2.0D0)+.5772156649015329D0
              CS0=(0.0D0,0.0D0)
              W0=0.0D0
              CR0=(1.0D0,0.0D0)
              DO 45 K=1,30
                 W0=W0+1.0D0/K
                 CR0=-0.25D0*CR0/(K*K)*Z2
45               CS0=CS0+CR0*W0
              CYV0=RP2*(CEC*CJV0-CS0)
           ENDIF
        ENDIF
        IF (N.EQ.0) N=1
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 50 K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1-CF2
           IF (K.LE.N) CBJ(K)=CF
           CF2=CF1
50         CF1=CF
        CS=CJV0/CF
        DO 55 K=0,N
55         CBJ(K)=CS*CBJ(K)
        IF (REAL(Z).LT.0.0D0) THEN
           CFAC0=CDEXP(PV0*CI)
           IF (DIMAG(Z).LT.0.0D0) THEN
              CYV0=CFAC0*CYV0-2.0D0*CI*DCOS(PV0)*CJV0
           ELSE IF (DIMAG(Z).GT.0.0D0) THEN
              CYV0=CYV0/CFAC0+2.0D0*CI*DCOS(PV0)*CJV0
           ENDIF
           DO 60 K=0,N
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBJ(K)=CDEXP(-PI*(K+V0)*CI)*CBJ(K)
              ELSE IF (DIMAG(Z).GT.0.0D0) THEN
                 CBJ(K)=CDEXP(PI*(K+V0)*CI)*CBJ(K)
              ENDIF
60         CONTINUE
           Z1=Z1
        ENDIF
        CBY(0)=CYV0
        DO 65 K=1,N
           CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           CBY(K)=CYY
65      CONTINUE
        CDJ(0)=V0/Z*CBJ(0)-CBJ(1)
        DO 70 K=1,N
70         CDJ(K)=-(K+V0)/Z*CBJ(K)+CBJ(K-1)
        CDY(0)=V0/Z*CBY(0)-CBY(1)
        DO 75 K=1,N
75         CDY(K)=CBY(K-1)-(K+V0)/Z*CBY(K)
        VM=N+V0
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
!
!       ==================================================
!       Purpose: Compute gamma function ‚(x)
!       Input :  x  --- Argument of ‚(x)
!                       ( x is not equal to 0,-1,-2,˙˙˙)
!       Output:  GA --- ‚(x)
!       ==================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,&
              -0.6558780715202538D0, -0.420026350340952D-1,&
               0.1665386113822915D0,-.421977345555443D-1,&
               -.96219715278770D-2, .72189432466630D-2,&
               -.11651675918591D-2, -.2152416741149D-3,&
               .1280502823882D-3, -.201348547807D-4,&
               -.12504934821D-5, .11330272320D-5,&
               -.2056338417D-6, .61160950D-8,&          
                .50020075D-8, -.11812746D-8,&          
                .1043427D-9, .77823D-11,&          
               -.36968D-11, .51D-12,&          
               -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward  
!                recurrence such that the magnitude of    
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point   
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
        end module specfun
      PROGRAM tunelRecub
! Solución analítica de la difración de una onda P incidente 
! en una cavidad cilíndrica circular con recubrimiento
      use datos; use vars_func_of_w; use Hank
      use debug; use RES; use plotter
      implicit none
      integer :: J,n,et,info,i
      real*8, pointer :: r,th
      integer, pointer :: reg
      complex*16, dimension(6,6) :: M
      complex*16, dimension(6,2) :: B !P , SV
      integer, dimension(6) :: ipiv 
      complex*16,dimension(nplot) :: pt_RES
      ! terms defined in the appendix (functions)
      complex*16 :: e11,e12,e21,e22,e41,e42,e71,e72,e81,e82,sig0
      real*8,dimension(2) :: BEALF
      call set_radios !a,b
      call setup_resu !puntos receptores
!     OMEI = - periodicdamper*PI/TW
      eta = radios(2)/radios(1) !b/a
      BEALF(1:2)=SQRT((0.5-NU(1:2))/(1.0-NU(1:2))) !IF POISSON RATIO IS GIVEN
      alf(1:2) = bet(1:2)/BEALF(1:2)
      print*,"alf", alf(1),alf(2)
      print*,"bet", bet(1),bet(2)
      print*,"muR=",(RHO(1) * bet(1)**2.)/(RHO(2) * bet(2)**2.)
      print*,"gamm=",alf(1) / alf(2)
      print*,"nu1=",NU(1)
      print*,"nu2=",NU(2)
      print*,"eta= b/a =",eta 
      print*,"DFREC=",DFREC
      print*,"Fmax= [Hz]", DFREC* NFREC
      write(6,'(A)', ADVANCE = "NO") " "
      write(6,'(A)', ADVANCE = "NO") repeat("_",58)
      write(6,'(A)', ADVANCE = "NO") " "
      print*,""
      write(6,'(A)', ADVANCE = "NO") repeat(" ",60)
      do J=1,NFREC !*********************************************
      FREC=DFREC*real(J); if (J .eq. 1)  FREC = 0.5_8 * DFREC ! Hz
      OME=2.0*PI*FREC !rad/s
!     COME = CMPLX(OME, OMEI,8)
!     VEL(1,1:2) = cmplx(alf(1:2)*& 
!                  (1.+1./pi/Qq*log(ome/2./pi)),-1./2./Qq,8)
!     VEL(2,1:2) = cmplx(bet(1:2)*& 
!                  (1.+1./pi/Qq*log(ome/2./pi)),-1./2./Qq,8)
      COME = OME * CMPLX(1., -1./(2.*Qq),8)
      VEL(1,1:2) = cmplx(alf(1:2),(/0.,0./),8)
      VEL(2,1:2) = cmplx(bet(1:2),(/0.,0./),8)
      cp(1:2) => VEL(1,1:2) ! dilatación
      cs(1:2) => VEL(2,1:2) ! corte
      w_c(1:2,1:2) = cOME / VEL(1:2,1:2) 
      beta(1:2) => w_c(2,1:2) !shear wave number
      alfa(1:2) => w_c(1,1:2) !compressional wave number
      p2(1:2) = (alfa(1:2))**2.0 !compressional wave number (square)
      s2(1:2) = (beta(1:2))**2.0 !shear wave number (square)
      aMU(1:2) = RHO(1:2) * cs(1:2)**2.
      Lambda(1:2) = RHO(1:2)* cp(1:2)**2. - real(2.)* aMU(1:2)
      muR = amu(1)/amu(2) !shear moduli ratio
      gammaP = z0! spacing variable  2.5D
      gammaS = z0! spacing variable
      call makeBessels(.false.) ! n = -1,0,1,2,...,nmax+1 (imprimir)
      do n=0,nmax*nfracs ! ensamblar matriz 4.26 y terminos independientes
      M = z0; B = z0; iPIV = 6
      ! sigma_{rr}1 = sigma_{rr}2   @ r = b
      M(1,1) = - muR * e11(3,1,1,2,n)
      M(1,2) = - muR * e12(3,2,1,2,n)
      M(1,3) = e11(3,1,2,2,n)
      M(1,4) = e11(4,1,2,2,n)
      M(1,5) = e12(3,2,2,2,n)
      M(1,6) = e12(4,2,2,2,n) 
      B(1,1) = (1.0) * et(n) * UI**n * muR * e11(1,1,1,2,n)
      ! sigma_{r th}1 = sigma_{r th}2   @ r = b
      M(2,1) = - muR * e41(3,1,1,2,n)
      M(2,2) = - muR * e42(3,2,1,2,n)
      M(2,3) = e41(3,1,2,2,n)
      M(2,4) = e41(4,1,2,2,n)
      M(2,5) = e42(3,2,2,2,n)
      M(2,6) = e42(4,2,2,2,n)
      B(2,1) = (1.0) * et(n) * UI**n * muR * e41(1,1,1,2,n)

      ! u_{r}1 = u_{r}2   @ r = b
      M(3,1) = - e71(3,1,1,2,n) 
      M(3,2) = - e72(3,2,1,2,n) 
      M(3,3) = e71(3,1,2,2,n)
      M(3,4) = e71(4,1,2,2,n) 
      M(3,5) = e72(3,2,2,2,n)
      M(3,6) = e72(4,2,2,2,n)
      B(3,1) = (1.0) * et(n) * UI**n * e71(1,1,1,2,n)
      ! u_{th}1 = u_{th}2   @ r = b
      M(4,1) = - e81(3,1,1,2,n) 
      M(4,2) = - e82(3,2,1,2,n) 
      M(4,3) = e81(3,1,2,2,n) 
      M(4,4) = e81(4,1,2,2,n) 
      M(4,5) = e82(3,2,2,2,n)
      M(4,6) = e82(4,2,2,2,n)
      B(4,1) = (1.0) * et(n) * UI**n * e81(1,1,1,2,n)
      ! sigma_{rr}2 = 0   @ r = a
      M(5,1) = z0
      M(5,2) = z0
      M(5,3) = e11(3,1,2,1,n)
      M(5,4) = e11(4,1,2,1,n)
      M(5,5) = e12(3,2,2,1,n)
      M(5,6) = e12(4,2,2,1,n)
      B(5,1) = z0
      ! sigma_{r th}2 = 0   @ r = a
      M(6,1) = z0
      M(6,2) = z0
      M(6,3) = e41(3,1,2,1,n)
      M(6,4) = e41(4,1,2,1,n)
      M(6,5) = e42(3,2,2,1,n)
      M(6,6) = e42(4,2,2,1,n)
      B(6,1) = z0

      call zgesv(6,1,M(1:6,1:6),6,IPIV,B(:,1),6,info)
      if(info .ne. 0) then
        write(6,'(A,I0,a,I0)', ADVANCE = "NO") &
        "se corta la suma en ",n, "system info = ",info
        exit
      else if (abs(B(1,1)) .lt. 0.00000001) then  !NaN
       write(6,'(A,I0,a)', ADVANCE = "NO") &
      "se corta la suma en ",n, " por chiquito"
       exit
      end if
!     call showMNmatrixZ(6,1,B,"  A  ",6)
      !elementos mecanicos
      do i = 1,nRes
       r => Rw(i)%r
       th => Rw(i)%th
       reg => Rw(i)%reg
        if (reg .eq. 1) then
          Rw(i)%s_rr_1(J) = Rw(i)%s_rr_1(J) + &
          (1. * et(n) * UI**n * e11(1,1,1,2+i,n) + &
                       B(1,1) * e11(3,1,1,2+i,n) + &
                       B(2,1) * e12(3,2,1,2+i,n)) * (cos(n * th))
          Rw(i)%s_tt_1(J) = Rw(i)%s_tt_1(J) + &
          (1. * et(n) * UI**n * e21(1,1,1,2+i,n) + & 
                       B(1,1) * e21(3,1,1,2+i,n) + & 
                       B(2,1) * e22(3,2,1,2+i,n)) * (cos(n * th))
          Rw(i)%s_rt_1(J) = Rw(i)%s_rt_1(J) + &
          (1. * et(n) * UI**n * e41(1,1,1,2+i,n) + & 
                       B(1,1) * e41(3,1,1,2+i,n) + & 
                       B(2,1) * e42(3,2,1,2+i,n)) * (sin(n * th))
          Rw(i)%u_r_1(J) = Rw(i)%u_r_1(J) + &
          (1. * et(n) * UI**n * e71(1,1,1,2+i,n) + & 
                       B(1,1) * e71(3,1,1,2+i,n) + & 
                       B(2,1) * e72(3,2,1,2+i,n)) * (cos(n * th))
          Rw(i)%u_t_1(J) = Rw(i)%u_t_1(J) + &
          (1. * et(n) * UI**n * e81(1,1,1,2+i,n) + & 
                       B(1,1) * e81(3,1,1,2+i,n) + & 
                       B(2,1) * e82(3,2,1,2+i,n)) * (sin(n * th))
        else if (reg .eq. 2) then
          Rw(i)%s_rr_2(J) = Rw(i)%s_rr_2(J) + &
                      (B(3,1) * e11(3,1,2,2+i,n) + &
                       B(4,1) * e11(4,1,2,2+i,n) + &
                       B(5,1) * e12(3,2,2,2+i,n) + &
                       B(6,1) * e12(4,2,2,2+i,n)) * (cos(n * th))
          Rw(i)%s_tt_2(J) = Rw(i)%s_tt_2(J) + &
                      (B(3,1) * e21(3,1,2,2+i,n) + &
                       B(4,1) * e21(4,1,2,2+i,n) + &
                       B(5,1) * e22(3,2,2,2+i,n) + &
                       B(6,1) * e22(4,2,2,2+i,n)) * (cos(n * th))
          Rw(i)%s_rt_2(J) = Rw(i)%s_rt_2(J) + &
                      (B(3,1) * e41(3,1,2,2+i,n) + &
                       B(4,1) * e41(4,1,2,2+i,n) + &
                       B(5,1) * e42(3,2,2,2+i,n) + &
                       B(6,1) * e42(4,2,2,2+i,n)) * (sin(n * th))
          Rw(i)%u_r_2(J) = Rw(i)%u_r_2(J) + &
                      (B(3,1) * e71(3,1,2,2+i,n) + &
                       B(4,1) * e71(4,1,2,2+i,n) + &
                       B(5,1) * e72(3,2,2,2+i,n) + &
                       B(6,1) * e72(4,2,2,2+i,n)) * (cos(n * th))
          Rw(i)%u_t_2(J) = Rw(i)%u_t_2(J) + &
                      (B(3,1) * e81(3,1,2,2+i,n) + &
                       B(4,1) * e81(4,1,2,2+i,n) + &
                       B(5,1) * e82(3,2,2,2+i,n) + &
                       B(6,1) * e82(4,2,2,2+i,n)) * (sin(n * th))
        end if
      end do! i:nRes
      end do !n
      
      ! términos fuera de la suma:
      do i = 1,nRes
       r => Rw(i)%r
       th => Rw(i)%th
       reg => Rw(i)%reg
       sig0 = amu(1) * beta(1)**2. !eq 3.15 
        if (reg .eq. 1) then
          Rw(i)%s_rr_1(J) = Rw(i)%s_rr_1(J) * 2. *amu(1) / r**2. 
          Rw(i)%s_tt_1(J) = Rw(i)%s_tt_1(J) * 2. *amu(1) / r**2. / sig0
          Rw(i)%s_rt_1(J) = Rw(i)%s_rt_1(J) * 2. *amu(1) / r**2. 
          Rw(i)%u_r_1(J) = Rw(i)%u_r_1(J) / r 
          Rw(i)%u_t_1(J) = Rw(i)%u_t_1(J) / r 
        else if (reg .eq. 2) then
          Rw(i)%s_rr_2(J) = Rw(i)%s_rr_2(J) * 2. *amu(2) / r**2. 
          Rw(i)%s_tt_2(J) = Rw(i)%s_tt_2(J) * 2. *amu(2) / r**2. / sig0
          Rw(i)%s_rt_2(J) = Rw(i)%s_rt_2(J) * 2. *amu(2) / r**2. 
          Rw(i)%u_r_2(J) = Rw(i)%u_r_2(J) / r 
          Rw(i)%u_t_2(J) = Rw(i)%u_t_2(J) / r 
        end if
      end do! i:nRes

!       write(6,'(A)', ADVANCE = "NO") repeat(char(8),60)
        print*,""
        write(6,'(A)', ADVANCE = "NO") "["
        write(6,'(A,A)', ADVANCE = "NO") & 
        repeat("X",58-int((58.0/NFREC)*(NFREC+1-J))), &
        repeat("_",int((58.0/NFREC)*(NFREC+1-J)))
        write(6,'(A)', ADVANCE = "NO") "]"
      abscisa(J) = dfrec * 2 * pi * J * radios(1) / cp(1)
      end do !J   
      ! plot curves
      ! end program tunelRecub
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
       do i = 1,nRes
        reg => Rw(i)%reg 
        if (reg .eq. 1) then
!         pt_RES = Rw(i)%s_rr_1(1:nplot)
!         call plot(i,1,reg, abscisa ,pt_RES,nplot)
          pt_RES = Rw(i)%s_tt_1(1:nplot)
          call plot(i,2,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)%s_rt_1(1:nplot)
!         call plot(i,3,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)% u_r_1(1:nplot)
!         call plot(i,4,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)% u_t_1(1:nplot)
!         call plot(i,5,reg, abscisa ,pt_RES,nplot)
        else if (reg .eq. 2) then
!         pt_RES = Rw(i)%s_rr_2(1:nplot)
!         call plot(i,1,reg, abscisa ,pt_RES,nplot)
          pt_RES = Rw(i)%s_tt_2(1:nplot)
          call plot(i,2,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)%s_rt_2(1:nplot)
!         call plot(i,3,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)% u_r_2(1:nplot)
!         call plot(i,4,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)% u_t_2(1:nplot)
!         call plot(i,5,reg, abscisa ,pt_RES,nplot)
        end if
      end do! i:nRes      
      end program tunelRecub
      function et(n)
      integer,intent(in) :: n
      integer :: et
      et = 2
      if (n .eq. 0) et = 1
      end function et
      
      function e11(i,c,reg,r,n)
      use vars_func_of_w, only : alfa,s2,gammaP
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e11
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(c)
      e11 = (n**2 + n - s2(reg) * radios(r)**2 / 2. + & 
      gammaP**2 * radios(r)**2) * Bess_n - alfa(reg)* radios(r) * Bess_n_1
      end function e11
      function e12(i,c,reg,r,n)
      use vars_func_of_w, only : beta
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e12
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(c)
      e12 = n * (-(n+1)* Bess_n + beta(reg) * radios(r) * Bess_n_1)
      end function e12
      function e21(i,c,reg,r,n)
      use vars_func_of_w, only : alfa,s2,p2
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e21
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(c)
      e21 = - (n**2 + n + s2(reg)*radios(r)**2 / 2. - p2(reg)*radios(r)**2) * &
             Bess_n + alfa(reg) * radios(r) * Bess_n_1
      end function e21
      function e22(i,c,reg,r,n)
      use vars_func_of_w, only : beta
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e22
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(c)
      e22 = n * ((n+1)*Bess_n - beta(reg)*radios(r)*Bess_n_1)
      end function e22
      function e41(i,c,reg,r,n)
      use vars_func_of_w, only : alfa
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e41
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(c)
      e41 = - n * (-(n+1) * Bess_n + alfa(reg)*radios(r)*Bess_n_1)
      end function e41
      function e42(i,c,reg,r,n)
      use vars_func_of_w, only : beta
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e42
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(c)
      e42 = - (n**2 + n - beta(reg)**2 * radios(r)**2 / 2.) * Bess_n + &
            beta(reg) * radios(r) * Bess_n_1
      end function e42
      function e71(i,c,reg,r,n)
      use vars_func_of_w, only : alfa
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e71
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(c)
      e71 = alfa(reg) * radios(r) * Bess_n_1 - n * Bess_n
      end function e71
      function e72(i,c,reg,r,n)
      use Hank, only : Bess
      implicit none
      complex*16 :: e72
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      e72 = n * Bess_n
      end function e72
      function e81(i,c,reg,r,n)
      use Hank, only : Bess
      implicit none
      complex*16 :: e81
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      e81 = - n * Bess_n
      end function e81
      function e82(i,c,reg,r,n)
      use vars_func_of_w, only : beta
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e82
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(c)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(c)
      e82 = - (beta(reg) * radios(r) * Bess_n_1 & 
              - n * Bess_n)
      end function e82
