Subroutine PotShape(l, rb, vb, curv_out, rmin, vmin)
  Use ccfull_initialization_mod
  
  Implicit None
  Integer :: n, i
  Integer, Intent(In) :: l
  Real(8), Intent(Out) :: rb, vb, curv_out, rmin, vmin
  Real(8) :: r, u0, u1, ra, rb_tmp, tolk, u, ddv00
  Real(8) :: dV_dr, V
  external dV_dr, V
  r = 50.5D0
  u0 = dV_dr(r, l)

  Do
    r = r - 1.0D0
    If (r < 0.0D0) Then
      rb = -5.0D0
      vb = -5.0D0
      curv_out = -1.0D0
      rmin = -1.0D0
      vmin = -1.0D0
      Return
    End If
    u1 = dV_dr(r, l)
    If (u0 * u1 <= 0.0D0) Exit
    u0 = u1
  End Do

  ra = r + 1.0D0
  rb_tmp = r
  tolk = 1.0D-6
  n = Int(Log10(Abs(rb_tmp - ra) / tolk) / Log10(2.0D0) + 0.5D0)

  Do i = 1, n
    r = (ra + rb_tmp) / 2.0D0
    u = dV_dr(r, l)
    If (u0 * u < 0.0D0) Then
      rb_tmp = r
    Else
      ra = r
    End If
  End Do

  rb = r
  ddv00 = (dV_dr(rb + 1.0D-5, l) - dV_dr(rb - 1.0D-5, l)) / (2.0D-5)

  If (ddv00 > 0.0D0) Then
    Print *, "Error: Second derivative positive at barrier top."
    Stop
  End If

  curv_out = Hbar * Sqrt(Abs(ddv00) / ReduceMass)
  vb = V(rb, l)

  ra = rb - 0.5D0
  u0 = dV_dr(ra, l)
  rb_tmp = 0.5D0
  u1 = dV_dr(rb_tmp, l)
  n = Int(Log10(Abs(rb - ra) / tolk) / Log10(2.0D0) + 0.5D0)

  Do i = 1, n
    r = (ra + rb_tmp) / 2.0D0
    u = dV_dr(r, l)
    If (u0 * u < 0.0D0) Then
      rb_tmp = r
    Else
      ra = r
    End If
  End Do

  rmin = r
  vmin = V(rmin, l)

End Subroutine PotShape

Real(8) Function Vn(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Real(8) :: R12
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  Vn = -V0 / (1.0D0 + Exp((r - R12) / A0))
End Function Vn

!Real(8) Function W(r)
!  Use ccfull_initialization_mod
!  Real(8), Intent(In) :: r
!  Real(8) :: rw
!  rw = 1.0D0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
!  W = 50.0D0 / (1.0D0 + Exp((r - rw) / 0.3D0))
!End Function W

Real(8) Function dVn_dr(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Real(8) :: R12, exp_term
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  exp_term = Exp((r - R12) / A0)
  dVn_dr = V0 / A0 * exp_term / (1.0D0 + exp_term)**2
End Function dVn_dr

Real(8) Function dVcent_dr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  dVcent_dr = -2.0D0 * l * (l + 1.0D0) * Hbar**2 / (2.0D0 * ReduceMass * r**3)
End Function dVcent_dr

Real(8) Function dVc_dr(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  dVc_dr = -Zpro * Ztar * Hbar / 137.0D0 / r**2
End Function dVc_dr

Real(8) Function dV_dr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Real(8)             :: dVc_dr, dVn_dr, dVcent_dr
  external dVc_dr, dVn_dr, dVcent_dr
  dV_dr = dVn_dr(r) + dVc_dr(r) + dVcent_dr(r, l)
End Function dV_dr

Real(8) Function Vc(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Vc = Zpro * Ztar * Hbar / 137.0D0 / r
End Function Vc

Real(8) Function Vr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Vr = h2m * l * (l + 1.0D0) / r**2
End Function Vr

Real(8) Function V(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Real(8)             :: Vn, Vc, Vr
  external Vn, Vc, Vr
  V = Vn(r) + Vc(r) + Vr(r, l)
End Function V

Real(8) Function VnCC(r, Xt)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r, Xt
  Real(8) :: R12
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  VnCC = -V0 / (1.0D0 + Exp((r - R12 - Xt) / A0))
End Function VnCC

Real(8) Function Fct(r)
  Use ccfull_initialization_mod
  IMPLICIT NONE
  Real(8), Intent(In) :: r
  Real(8) ::  result, Lambda

  Lambda = LambdaT
  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  result = result * BetaT / Sqrt(4.D0 * Pi)

  Fct = result
  Return
End Function Fct

Real(8) Function Fct2(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 2

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  Fct2 = result
  Return
End Function Fct2

Real(8) Function Fct3(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 3

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  result = result * (Beta4T + 7.D0 * Beta2T**2 / 7.D0 / Sqrt(Pi))

  Fct3 = result
  Return
End Function Fct3

Real(8) Function Fct4(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 4

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  result = result * (Beta4T + 9.D0 * Beta2T**2 / 7.D0 / Sqrt(Pi))

  Fct4 = result
  Return
End Function Fct4

Real(8) Function Fcp(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result, Lambda

  Lambda = LambdaP

  If (r > Rpro) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro) ** Lambda
  End If

  result = result * BetaP / SQRT(4.D0 * Pi)

  Fcp = result
  Return
End Function Fcp

Real(8) Function Fcp2(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  Fcp2 = result
  Return
End Function Fcp2

Real(8) Function Fcp3(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 3
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  result = result * (Beta4P + 7.D0 * Beta2P**2 / 7.D0 / Sqrt(Pi))

  Fcp3 = result
  Return
End Function Fcp3

Real(8) Function Fcp4(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 4
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  result = result * (Beta4P + 9.D0 * Beta2P**2 / 7.D0 / Sqrt(pi))

  Fcp4 = result
  Return
End Function Fcp4

Real(8) Function Fctt(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = LambdaT2
  result = 0.D0

  If (r > Rtar) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  result = result * BetaT2 / SQRT(4.D0 * Pi)

  Fctt = result
  Return
End Function Fctt

Real(8) Function Fct2v(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2
  result = 0.D0

  If (r > Rtar) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  Fct2v = result
  Return
End Function Fct2v

Real(8) Function Fcp2v(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  Fcp2v = result
  Return
End Function Fcp2v

Real(8) Function Ftrans(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result, dVn_dr
  external dVn_dr

  result = 0.D0
  result = Ftr * dVn_dr(r)

  Ftrans = result
  Return
End Function Ftrans

Real(8) Function CG(j1, m1, j2, m2, j3, m3)
    ! Clebsch-Gordan coefficient <j1 m1 j2 m2 | j3 m3>
    implicit real*8 (a-h, o-z)
    external fact

    if (m1 + m2 .ne. m3) then
        cg = 0.d0
        return
    end if

    if (j3 .lt. abs(j1 - j2)) then
        cg = 0.d0
        return
    end if

    if (j3 .gt. j1 + j2) then
        cg = 0.d0
        return
    end if

    ka = j1 + j2 - j3
    kb = j3 + j1 - j2
    kc = j2 + j3 - j1
    kd = j1 + j2 + j3 + 1

    del = sqrt(fact(ka) * fact(kb) * fact(kc) / fact(kd))

    cg = 0.d0
    do n = 0, max(j1 + j2 - j3, j1 - m1, j2 + m2)
        ka1 = j1 + j2 - j3 - n
        if (ka1 .lt. 0.d0) cycle

        ka2 = j3 - j2 + m1 + n
        if (ka2 .lt. 0.d0) cycle

        ka3 = j3 - j1 - m2 + n
        if (ka3 .lt. 0.d0) cycle

        ka4 = j1 - m1 - n
        if (ka4 .lt. 0.d0) cycle

        ka5 = j2 + m2 - n
        if (ka5 .lt. 0.d0) cycle

        cg = cg + (-1.d0)**n / &
                 (fact(n) * fact(ka1) * fact(ka2) * fact(ka3) * fact(ka4) * fact(ka5))
    end do

    cg = cg * sqrt(fact(j1 + m1) * fact(j1 - m1))
    cg = cg * sqrt(fact(j2 + m2) * fact(j2 - m2))
    cg = cg * sqrt(fact(j3 + m3) * fact(j3 - m3))

    cg = cg * sqrt(2.d0 * j3 + 1.d0) * del

    return
End Function CG

Subroutine mdiag(a,n,d,v)

	implicit real*8(a-h,o-z)
	parameter(nmax=100)
  parameter (nlevelmax=30)
  dimension a(nlevelmax,nlevelmax)
	dimension b(nmax),z(nmax),d(nlevelmax),v(nlevelmax,nlevelmax)

	do 12 ip=1,n
	do 11 iq=1,n
	v(ip,iq)=0.d0
 11	continue
	v(ip,ip)=1.d0
 12	continue

	do 13 ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.d0
 13	continue

	nrot=0

	do 24 i=1,50
	sm=0.d0

	do 15 ip=1,n-1
	do 14 iq=ip+1,n
	sm=sm+abs(a(ip,iq))
 14	continue
 15     continue

	if(sm.eq.0.d0) return

	if(i.lt.4) then
	tresh=0.2d0*sm/n**2
	else
	tresh=0.d0
	endif

	do 22 ip=1,n-1
	do 21 iq=ip+1,n
	g=100.d0*abs(a(ip,iq))
	
	if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq)))) then
	   a(ip,iq)=0.d0
	elseif(abs(a(ip,iq)).gt.tresh) then
	   h=d(iq)-d(ip)

	   if(abs(h)+g.eq.abs(h)) then
	     t=a(ip,iq)/h
	   else
	     theta=0.5d0*h/a(ip,iq)
	     t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
	     if(theta.lt.0.d0) t=-t
	   endif

	c=1.d0/sqrt(1.d0+t**2)
	s=t*c
	tau=s/(1.d0+c)
	h=t*a(ip,iq)
	z(ip)=z(ip)-h
	z(iq)=z(iq)+h
	d(ip)=d(ip)-h
	d(iq)=d(iq)+h
	a(ip,iq)=0.d0

	do 16 j=1,ip-1
	g=a(j,ip)
	h=a(j,iq)
	a(j,ip)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 16	continue

	do 17 j=ip+1,iq-1
	g=a(ip,j)
	h=a(j,iq)
	a(ip,j)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 17	continue

	do 18 j=iq+1,n
	g=a(ip,j)
	h=a(iq,j)
	a(ip,j)=g-s*(h+g*tau)
	a(iq,j)=h+s*(g-h*tau)
 18	continue

	do 19 j=1,n
	g=v(j,ip)
	h=v(j,iq)
	v(j,ip)=g-s*(h+g*tau)
	v(j,iq)=h+s*(g-h*tau)
 19	continue
	
	nrot=nrot+1
	endif

 21	continue
 22	continue

	do 23 ip=1,n
	b(ip)=b(ip)+z(ip)
	d(ip)=b(ip)
	z(ip)=0.d0

 23	continue
 24	continue

	return

End Subroutine  mdiag

SUBROUTINE matinv(nmax, c, d)
  IMPLICIT NONE
  
  ! Argument declarations
  INTEGER, INTENT(IN) :: nmax
  COMPLEX(8), DIMENSION(30,30), INTENT(INOUT) :: c, d
  
  ! Parameter and local variable declarations
  COMPLEX(8) :: u(30), v(30), t, a, b
  REAL(8) :: deter
  INTEGER :: m, n, j, k, l
  
  ! Initialize
  deter = 1.0_8
  d = (0.0_8, 0.0_8)
  DO m = 1, nmax
      d(m, m) = (1.0_8, 0.0_8)
  END DO
  
  ! Main inversion loop
  DO n = 1, nmax
      t = c(n, n)
      
      ! Check for zero pivot element
      IF (ABS(t) < 1.0E-10_8) THEN
          j = n
          DO
              j = j + 1
              IF (j > nmax) THEN
                  PRINT *, 'matrix not invertible'
                  RETURN
              END IF
              
              t = c(n, j)
              deter = -deter
              IF (ABS(t) > 1.0E-10_8) EXIT
          END DO
          
          ! Swap rows
          DO k = 1, nmax
              u(k) = c(n, k)
              v(k) = d(n, k)
              c(n, k) = c(j, k)
              d(n, k) = d(j, k)
              c(j, k) = u(k)
              d(j, k) = v(k)
          END DO
      END IF
      
      ! Eliminate column
      DO k = 1, nmax
          IF (k == n) CYCLE
          
          a = c(k, n) / c(n, n)
          DO l = 1, nmax
              c(k, l) = c(k, l) - a * c(n, l)
              d(k, l) = d(k, l) - a * d(n, l)
          END DO
      END DO
      
      ! Normalize row
      b = c(n, n)
      deter = b * deter
      DO m = 1, nmax
          c(n, m) = c(n, m) / b
          d(n, m) = d(n, m) / b
      END DO
  END DO
  
END SUBROUTINE matinv


Real(8) Function fact(n)
    Implicit none
    Integer n, i

    If (n .lt. 0) Then
        fact = 0.d0
        Return
    End If

    If (n .eq. 0) Then
        fact = 1.d0
        Return
    End If

    fact = 1.d0
    Do i = 1, n
        fact = fact * i * 1.d0
    End Do

    Return
End Function fact

Subroutine grotation()
    Use ccfull_initialization_mod
    Implicit None
    Integer :: it, jt, ip, jp, i
    Character(len=1) :: ans
    Real(8) :: erott0, erotp0

    Do it = 0, Ntar
      Do jt = 0, Ntar
        bett(it,jt) = 0.0D0
        If (it == jt - 1 .Or. it == jt + 1 .Or. it == jt) bett(it,jt) = Beta2T
      End Do
      erott(it) = 2.0D0 * it * (2.0D0 * it + 1.0D0) / 6.0D0 * E2T
    End Do

    If (IVIBROTT == 1 .And. Ntar > 1) Then
      Print *, ''
      Print *, 'Generalised E2 couplings for the target rotor (n/y)?'
      Read(*,'(A)') ans
      If (ans == 'y' .Or. ans == 'Y') Then
        Print *, 'Beta2_targ=', Beta2T
        Do it = 1, Ntar - 1
          jt = it + 1
          Print *, 'Modify beta2 for transition from', 2*it, '+ to', 2*jt, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) bett(it,jt)
            bett(jt,it) = bett(it,jt)
            Print *, 'Beta2 set to: ', bett(it,jt)
          End If
        End Do

        Print *, 'Reorientation terms:'
        Do it = 1, Ntar
          jt = it
          Print *, 'Modify beta2 for transition from', 2*it, '+ to', 2*jt, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) bett(it,jt)
            bett(jt,it) = bett(it,jt)
          End If
        End Do

        Print *, 'Excitation energies:'
        Do i = 2, Ntar
          Print *, 'Energy of the', 2*i, '+ state for pure rotor =', erott(i)
          Print *, 'Modify this energy? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            erott0 = erott(i)
            Print *, 'Energy=?'
            Read(*,*) erott(i)
            Print *, 'Modified from', erott0, ' to ', erott(i)
          End If
        End Do
      End If
    End If

    Do ip = 0, Npro
      Do jp = 0, Npro
        betp(ip,jp) = 0.0D0
        If (ip == jp - 1 .Or. ip == jp + 1 .Or. ip == jp) betp(ip,jp) = Beta2P
      End Do
      erotp(ip) = 2.0D0 * ip * (2.0D0 * ip + 1.0D0) / 6.0D0 * E2P
    End Do

    If (IVIBROTP == 1 .And. Npro > 1) Then
      Print *, 'Generalised E2 couplings for the projectile rotor (n/y)?'
      Read(*,'(A)') ans
      If (ans == 'y' .Or. ans == 'Y') Then
        Print *, 'Beta2_proj=', Beta2P
        Do ip = 1, Npro - 1
          jp = ip + 1
          Print *, 'Modify beta2 for transition from', 2*ip, '+ to', 2*jp, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) betp(ip,jp)
            betp(jp,ip) = betp(ip,jp)
          End If
        End Do

        Print *, 'Reorientation terms:'
        Do ip = 1, Npro
          jp = ip
          Print *, 'Modify beta2 for transition from', 2*ip, '+ to', 2*jp, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) betp(ip,jp)
            betp(jp,ip) = betp(ip,jp)
          End If
        End Do

        Print *, 'Excitation energies:'
        Do i = 2, Npro
          Print *, 'Energy of the', 2*i, '+ state for pure rotor =', erotp(i)
          Print *, 'Modify this energy? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            erotp0 = erotp(i)
            Print *, 'Energy=?'
            Read(*,*) erotp(i)
            Print *, 'Modified from', erotp0, ' to ', erotp(i)
          End If
        End Do
      End If
    End If

End Subroutine grotation

Subroutine Mutual()
    Use ccfull_initialization_mod
    Implicit None
    Integer :: i, j
    Character(len=1) :: ans, sign1, sign2
    Integer :: imut

    Do i = 0, Nlevelmax
      Do j = 0, Nlevelmax
        imutual(i,j) = 0
      End Do
    End Do

    Print *, ''
    Print *, 'Mutual excitations in the *target* nucleus'
    Print *, 'Include the mutual excitations (y/n)?'
    Read(*,'(A)') ans

    If (ans == 'n' .Or. ans == 'N') Then
      imut = 0
      Nlevel = (Ntar + NphononT2 + 1) * (Npro + 1)
      Do i = 0, Ntar
        Do j = 0, NphononT2
          If (i == 0 .Or. j == 0) imutual(i,j) = 1
        End Do
      End Do
      Return
    End If

    Print *, 'All the possible mutual excitation channels (n/y)?'
    Read(*,'(A)') ans

    If (ans == 'y' .Or. ans == 'Y') Then
      imut = 1
      Nlevel = (Ntar + 1) * (NphononT2 + 1) * (Npro + 1)
      Do i = 0, Ntar
        Do j = 0, NphononT2
          imutual(i,j) = 1
        End Do
      End Do
      Return
    End If

    imut = 2
    If (Mod(LambdaT, 2) == 0) Then
      sign1 = '+'
    Else
      sign1 = '-'
    End If
    If (Mod(LambdaT2, 2) == 0) Then
      sign2 = '+'
    Else
      sign2 = '-'
    End If

    Do i = 0, Ntar
      Do j = 0, NphononT2
        If (i == 0 .Or. j == 0) Then
          imutual(i,j) = 1
        Else
          Print *, 'Include (', LambdaT, sign1, '^', i, ',', LambdaT2, sign2, '^', j, ') state ? (y/n)'
          Read(*,'(A)') ans
          If (ans == 'n' .Or. ans == 'N') Then
            imutual(i,j) = 0
          Else
            imutual(i,j) = 1
          End If
        End If
      End Do
    End Do

    Print *, 'Excited states in the target to be included:'
    Nlevel = 0

    Do i = 0, Ntar
      Do j = 0, NphononT2
        If (imutual(i,j) == 0) Cycle
        Print *, '(', LambdaT, sign1, '^', i, ',', LambdaT2, sign2, '^', j, ') state'
        Nlevel = Nlevel + 1
      End Do
    End Do

    Nlevel = Nlevel * (Npro + 1)

End Subroutine Mutual

Subroutine Anharmonicity()
  Use ccfull_initialization_mod
  Implicit None

  Character(1) :: ans, sign
  Integer :: it, jt, i

  ! Target: AHV Couplings (1st mode)
  Do it = 0, Ntar
    Do jt = 0, Ntar
      betnahv(it, jt) = 0.0
      betcahv(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahv(it, jt) = BetaTn * Sqrt(Real(jt))
        betcahv(it, jt) = BetaT  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahv(it, jt) = BetaTn * Sqrt(Real(it))
        betcahv(it, jt) = BetaT  * Sqrt(Real(it))
      End If
    End Do
    omeahv(it) = it * OmegaT
  End Do

  If (IVIBROTT == 0 .And. Ntar > 1) Then
    Write(*, '(A)') 'AHV couplings for the first mode in the target phonon (n/y)?'
    Read(*, '(A1)') ans
    If (ans == 'Y' .Or. ans == 'y') Then
      Write(*, '(A)') '**** AHV Couplings in the target (the 1st mode)'
      If ((-1.)**LambdaT ==  1) sign = '+'
      If ((-1.)**LambdaT == -1) sign = '-'

      Do it = 0, Ntar
        Do jt = 0, Ntar
          If (it > jt)Cycle
          If ((it == 0 .And. jt <= 1)) Cycle
          Write(*, *) ' '
          Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaT, sign, '^', it, ' to the', LambdaT, sign, '^', jt, 'state:'
          Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahv(it, jt), betcahv(it, jt)
          Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   beta_N and beta_C =?'
            Read(*, *) betnahv(it, jt), betcahv(it, jt)
            Write(*, '(G15.5, G15.5)') betnahv(it, jt), betcahv(it, jt)
          End If
        End Do
      End Do

      Write(*, *) 'Excitation energy for the first mode:'
      Do i = 2, Ntar
        Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaT
        Write(*, *) '   Modify this energy(n/y)?'
        Read(*, '(A1)') ans
        If (ans == 'Y' .Or. ans == 'y') Then
          Write(*, *) '   Energy=?'
          Read(*, *) omeahv(i)
          Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahv(i), '(HO: ', i * OmegaT, ')'
        End If
      End Do
    End If
  End If


  ! Target: AHV Couplings (2st mode)
  Do it = 0, NphononT2
    Do jt = 0, NphononT2
      betnahv2(it, jt) = 0.0
      betcahv2(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahv2(it, jt) = BetaT2n * Sqrt(Real(jt))
        betcahv2(it, jt) = BetaT2  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahv(it, jt) = BetaT2n * Sqrt(Real(it))
        betcahv(it, jt) = BetaT2  * Sqrt(Real(it))
      End If
    End Do
    omeahv2(it) = it * OmegaT2
  End Do

  If (NphononT2 > 1) Then
    Write(*, '(A)') 'AHV couplings for the second mode in the target phonon (n/y)?'
    Read(*, '(A1)') ans
    If (ans == 'Y' .Or. ans == 'y') Then
      Write(*, '(A)') '**** AHV Couplings in the target (the 1st mode)'
      If ((-1.)**LambdaT2 ==  1) sign = '+'
      If ((-1.)**LambdaT2 == -1) sign = '-'

      Do it = 0, NphononT2
        Do jt = 0, NphononT2
          If (it > jt)Cycle
          If ((it == 0 .And. jt <= 1)) Cycle
          Write(*, *) ' '
          Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaT2, sign, '^', it, ' to the', LambdaT2, sign, '^', jt, 'state:'
          Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahv2(it, jt), betcahv2(it, jt)
          Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   beta_N and beta_C =?'
            Read(*, *) betnahv2(it, jt), betcahv2(it, jt)
            Write(*, '(G15.5, G15.5)') betnahv2(it, jt), betcahv2(it, jt)
          End If
        End Do
      End Do

      Write(*, *) 'Excitation energy for the second mode:'
      Do i = 2, Ntar
        Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaT2
        Write(*, *) '   Modify this energy(n/y)?'
        Read(*, '(A1)') ans
        If (ans == 'Y' .Or. ans == 'y') Then
          Write(*, *) '   Energy=?'
          Read(*, *) omeahv2(i)
          Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahv2(i), '(HO: ', i * OmegaT2, ')'
        End If
      End Do
    End If
  End If


  Do it = 0, Npro
    Do jt = 0, Npro
      betnahvp(it, jt) = 0.0
      betcahvp(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahvp(it, jt) = BetaPn * Sqrt(Real(jt))
        betcahvp(it, jt) = BetaP  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahvp(it, jt) = BetaPn * Sqrt(Real(it))
        betcahvp(it, jt) = BetaP  * Sqrt(Real(it))
      End If
    End Do
    omeahvp(it) = it * OmegaP
  End Do

  If (IVIBROTP == 0 .And. Npro > 1) Then
      Write(*, '(A)') 'AHV couplings for the first mode in the projectile phonon (n/y)?'
      Read(*, '(A1)') ans
      If (ans == 'Y' .Or. ans == 'y') Then
        Write(*, '(A)') '**** AHV Couplings in the projectile'
        If ((-1.)**LambdaP ==  1) sign = '+'
        If ((-1.)**LambdaP == -1) sign = '-'

        Do it = 0, Npro
          Do jt = 0, Npro
            If (it > jt)Cycle
            If ((it == 0 .And. jt <= 1)) Cycle
            Write(*, *) ' '
            Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaP, sign, '^', it, ' to the', LambdaP, sign, '^', jt, 'state:'
            Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahvp(it, jt), betcahvp(it, jt)
            Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
            Read(*, '(A1)') ans
            If (ans == 'Y' .Or. ans == 'y') Then
              Write(*, *) '   beta_N and beta_C =?'
              Read(*, *) betnahvp(it, jt), betcahvp(it, jt)
              Write(*, '(G15.5, G15.5)') betnahvp(it, jt), betcahvp(it, jt)
            End If
          End Do
        End Do

        Write(*, *) 'Excitation energy for the first mode:'
        Do i = 2, Npro
          Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaP
          Write(*, *) '   Modify this energy(n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   Energy=?'
            Read(*, *) omeahvp(i)
            Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahvp(i), '(HO: ', i * OmegaP, ')'
          End If
        End Do
      End If
    End If

End Subroutine Anharmonicity

Subroutine coupled_matrix0()
  Use ccfull_initialization_mod
  implicit none
  real(8) :: A(Nlevelmax, Nlevelmax)
  Real(8) :: C, CG
  integer :: i, j, ip, it, it2, jp, jt, jt2
  external CG

  !initialization of A matrix
  Do i = 1, Nlevelmax
    Do j = 1, Nlevelmax
     A(i,j) = 0.0d0
    End Do
  End Do

  If((Nlevel - Ntrans) == 1)Return
  i = 0
  Do ip = 0, Npro
    Do it = 0, Ntar
      Do it2 = 0, NphononT2

        if (NphononT2 .Ne. 0)Then
          if(imutual(it, it2) == 0)Cycle
        End if

        i = i + 1 
        j = 0
        Do jp = 0, Npro
          Do jt = 0, Ntar
            Do jt2 = 0, NphononT2

              if (NphononT2 .Ne. 0)Then
                if(imutual(jt, jt2) == 0)Cycle
              End if

              j = j + 1 

              If (i > j)Then
                A(i,j) = A(j,i)
                Cycle
              End If

              C = 0.0d0

              if (ip == jp .and. it2 == jt2) then
                  if (IVIBROTT == 0) then
                      C = Rtar * betnahv(it, jt) / sqrt(4.0d0 * pi)

                  else
                      C = Rtar * bett(it, jt) * sqrt((2 * 2 * it + 1) * 5.0d0 * (2 * 2 * jt + 1) / (4.0d0 * pi)) &
                          * CG(2 * it, 0, 2, 0, 2 * jt, 0)**2 / (2.0d0 * 2 * jt + 1)

                      C = C + Rtar * Beta4T * sqrt((2 * 2 * it + 1) * 9.0d0 * (2 * 2 * jt + 1) / (4.0d0 * pi)) &
                          * CG(2 * it, 0, 4, 0, 2 * jt, 0)**2 / (2.0d0 * 2 * jt + 1)
                  end if
              end if

              if (ip == jp .and. it == jt) then
                  C = C + Rtar * betnahv2(it2, jt2) / sqrt(4.0d0 * pi)
              end if

              if (it == jt .and. it2 == jt2) then
                  if (IVIBROTP == 0) then
                      C = C + Rpro * betnahvp(ip, jp) / sqrt(4.0d0 * pi)
                  else
                      C = C + Rpro * betp(ip, jp) * sqrt((2 * 2 * ip + 1) * 5.0d0 * (2 * 2 * jp + 1) / (4.0d0 * pi)) &
                          * CG(2 * ip, 0, 2, 0, 2 * jp, 0)**2 / (2.0d0 * 2 * jp + 1)

                      C = C + Rpro * Beta4P * sqrt((2 * 2 * ip + 1) * 9.0d0 * (2 * 2 * jp + 1) / (4.0d0 * pi)) &
                          * CG(2 * ip, 0, 4, 0, 2 * jp, 0)**2 / (2.0d0 * 2 * jp + 1)
                  end if
              end if

              A(i,j) = C


            End Do
          End Do
        End Do

      End Do
    End Do
  End Do

  Call mdiag(A, Nlevel - Ntrans, Ev, Evec)

  Return

End Subroutine coupled_matrix0

Subroutine coupled_matrix(r, cpot_matrix)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  real(8), intent(out) :: cpot_matrix(Nlevelmax, Nlevelmax)
  Real(8) :: C, A, C0
  Real(8) :: Fct, Fct2, Fct3, Fct4, Fct2v, Fctt
  Real(8) :: Ftrans, Fcp, Fcp2, Fcp3, Fcp4, Fcp2v
  Real(8) :: CG, Vn, VnCC
  integer :: i, j, k, ip, it, it2, jp, jt, jt2
  external CG, Fct, Fct2, Fct3, Fct4, Fct2v, Fctt
  external Ftrans, Fcp, Fcp2, Fcp3, Fcp4, Fcp2v, Vn, VnCC
  
  Do i = 1, nlevelmax
    eps(i) = 0.
  End Do

  Do i = 1, Nlevelmax
    Do j = 1, Nlevelmax
     cpot_matrix(i,j) = 0.0d0
    End Do
  End Do

  If((Nlevel - Ntrans) == 1)Return

  i = 0
  Do ip = 0, Npro
    Do it = 0, Ntar
      Do it2 = 0, NphononT2

        if (NphononT2 .Ne. 0)Then
          if(imutual(it, it2) == 0.)Cycle
        End if

        i = i + 1 
        j = 0
        Do jp = 0, Npro
          Do jt = 0, Ntar
            Do jt2 = 0, NphononT2

            if (NphononT2 .Ne. 0)Then
              if(imutual(jt, jt2) == 0.)Cycle
            End if

            j = j + 1 

            If (i > j)Then

              cpot_matrix(i,j) = cpot_matrix(j,i)
              Cycle
            End If

            C = 0.0d0
            ! nuclear coupling
            Do k = 1, Nlevel-Ntrans
              C = C + VnCC(r, Ev(k))*Evec(i, k)*Evec(j, k)

            End Do

            ! ---- Target contribution ----
            if (ip == jp .and. it2 == jt2) then
                if (IVIBROTT == 0) then
                    if (it /= jt) then
                        A = betcahv(it, jt) * Fct(r)
                        if (BetaT /= 0.d0) A = A / BetaT
                        
                    else
                        A = betcahv(it, jt) * Fct2v(r) / sqrt(4.d0 * pi)

                    end if
                    C = C + A
                else
                    C = C + sqrt((2 * 2 * it + 1) * 5.d0 * (2 * 2 * jt + 1) / (4.d0 * pi)) &
                          * CG(2 * it, 0, 2, 0, 2 * jt, 0)**2 / (2.d0 * 2 * jt + 1) * Fct2(r) &
                          * (bett(it, jt) + 2.d0 * sqrt(5.d0 / pi) * bett(it, jt)**2 / 7.d0)

                    C = C + sqrt((2 * 2 * it + 1) * 9.d0 * (2 * 2 * jt + 1) / (4.d0 * pi)) &
                          * CG(2 * it, 0, 4, 0, 2 * jt, 0)**2 / (2.d0 * 2 * jt + 1) * Fct4(r)
                end if
            end if

            ! ---- Projectile contribution ----
            if (it == jt .and. it2 == jt2) then
                if (IVIBROTP == 0) then
                    if (ip /= jp) then
                        A = betcahvp(ip, jp) * Fcp(r)
                        if (BetaP /= 0.d0) A = A / BetaP
                    else
                        A = betcahvp(ip, jp) * Fcp2v(r) / sqrt(4.d0 * pi)
                    end if
                    C = C + A
                else
                    C = C + sqrt((2 * 2 * ip + 1) * 5.d0 * (2 * 2 * jp + 1) / (4.d0 * pi)) &
                          * CG(2 * ip, 0, 2, 0, 2 * jp, 0)**2 / (2.d0 * 2 * jp + 1) * Fcp2(r) &
                          * (betp(ip, jp) + 2.d0 * sqrt(5.d0 / pi) * betp(ip, jp)**2 / 7.d0)

                    C = C + sqrt((2 * 2 * ip + 1) * 9.d0 * (2 * 2 * jp + 1) / (4.d0 * pi)) &
                          * CG(2 * ip, 0, 4, 0, 2 * jp, 0)**2 / (2.d0 * 2 * jp + 1) * Fcp4(r)
                end if
            end if

            ! ---- BetaT2 part (vibrational excitation of target) ----
            if (ip == jp .and. it == jt) then
                if (it2 /= jt2) then
                    A = betcahv2(it2, jt2) * Fctt(r)
                    if (BetaT2 /= 0.d0) A = A / BetaT2
                else
                    A = betcahv2(it2, jt2) * Fct2v(r) / sqrt(4.d0 * pi)
                end if
                C = C + A
            end if

            ! ---- Excitation energy shift ----
            If (it == jt .and. ip == jp .and. it2 == jt2) then
                if (IVIBROTT == 0) then
                    C = C + omeahv(it)
                    eps(i) = eps(i) + omeahv(it)
                else
                    C = C + erott(it)
                    eps(i) = eps(i) + erott(it)
                end if

                if (IVIBROTP == 0) then
                    C = C + omeahvp(ip)
                    eps(i) = eps(i) + omeahvp(ip)
                else
                    C = C + erotp(ip)
                    eps(i) = eps(i) + erotp(ip)
                end if

                C = C + omeahv2(it2)
                eps(i) = eps(i) + omeahv2(it2)
            End if

            cpot_matrix(i, j) = C


            End Do
          End Do
        End Do

      End Do
    End Do
  End Do

  C0 = cpot_matrix(1,1)
  Do i = 1, Nlevel-Ntrans
    cpot_matrix(i,i) = cpot_matrix(i,i) - Vn(r)
    !cpot_matrix(i,i) = cpot_matrix(i,i) - C0
  End Do

  ! transfer coupling
  If (Ntrans == 1) Then
      cpot_matrix(1, Nlevel) = Ftrans(r)
      cpot_matrix(Nlevel, 1) = Ftrans(r)
      cpot_matrix(Nlevel, Nlevel) = cpot_matrix(1, 1) - Qtrans + (Zpro + iq) * (Ztar - iq) / r * Hbar / 137.0 - Zpro * Ztar / r * Hbar / 137.0
  End If

  Return

End Subroutine coupled_matrix