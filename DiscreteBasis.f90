Subroutine DiscreteBasis(P)
  Use ccfull_initialization_mod

  Implicit None
  real(8), intent(out) :: P
  Real(8) ::  k, kk
  Integer ::  i, j
  complex*16 Hpsi(Nx1, Nx1)
  complex*16 Hpsi0(Nx1, 1)
  complex*16 cc(Nx1)

  Call Prameterize_H(Hamilton, Ec_psi, Ec_psi0)

  Call MATMUL_CPLX(Hamilton, Ec_psi,  Hpsi, Nx1,  Nx2, Nx1)
  Call MATMUL_CPLX(Hamilton, Ec_psi0, Hpsi0, Nx1, Nx2, 1)
  call CSOLVE_EQUATION(Nx1,Hpsi,Hpsi0, cc)

  P = 0.
  Do i = 1, Nlevel
    if(Echannel(i) .gt. 0)then
      !k  = acos(1. - Echannel(i)/2. / t)/dr
      !kk = acos(1. - E/2. / t)/dr
      k=sqrt((2.d0*ReduceMass/Hbar**2*Echannel(i)))
      kk=sqrt(2.d0*ReduceMass/hbar**2*E)
      P = P+(abs(cc(i)))**2*k/kk
    end if 
  end do

  Do i = 1, Nx2
    wf(i) = Ec_psi0(i)
  End Do

  Do i = 1, Nx1
    Do j = 1, Nx2
      wf(j) = wf(j) + cc(i) * Ec_psi(j, i)
    End Do
  End Do
  !wf = wf / wf(R_iterat)
End Subroutine


Subroutine Prameterize_H(H2, dbpsi, dbpsi0)
  Use ccfull_initialization_mod
  implicit none
  !H2, dbpsi, dbpsi0
  complex*16 H1(Nx_cc,Nx_cc)
  complex*16 H2(Nx1, Nx2)
  complex*16 dbpsi(Nx2, Nx1)
  complex*16 dbpsi0(Nx2)
  complex*16 ai

  Integer ::  i, j, ir, ie, il, idx, jdx
  Real(8) ::  r,rh, V, eta, rho, k, ec
  Real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  Real(8), dimension(0:200) :: sigmad
  Integer, dimension(0:200) :: iexp
  external V
  H1 = (0.d0, 0.d0)
  H2 = (0.d0, 0.d0)
  ai = (0.d0,1.d0)

  !BetaT = para
  ! contruct the cpot
  Call coupled_matrix0

  rh = R_min + dr / 2.0
  Call coupled_matrix(rh, CPOTH)

  Do ir = 0, R_iterat + 1
    r = R_min + dr * ir

    Call coupled_matrix(r, CPOT0)
    Do  i = 1, Nlevel
      Do  j = 1, Nlevel
        CPOT(i, j, ir) = CPOT0(i, j)

        ! contruct the H, potential term
        If(ir .lt. R_iterat)Then

          idx = Nlevel * ir + i
          jdx = Nlevel * ir + j
          H1(idx, jdx) = H1(idx, jdx) + CPOT0(i, j)
          
          
          If(i == j)then
            H1(idx, jdx) = H1(idx, jdx) + 2*t - E
          End If

        End If 

      End Do
    End Do

    If(ir .lt. R_iterat)Then
      Do i = 1, Nlevel
          idx = i + ir* Nlevel
          H1(idx, idx) = H1(idx, idx) + V(r, L_i) 

      End Do
    End If
  End Do

  ! contruct the H, knetic term 
  Do i = 1, Nlevel*(R_iterat-1)
        H1(i, i+nlevel) = -t
        H1(i+nlevel, i) = -t
  End Do

  Do i = 1, Nlevel
    Echannel(i) = E - V(R_min, L_i) - CPOT(i,i,0)
    
    If(Echannel(i) > 0)Then
     
      !k = acos(1. - Echannel(i)/2. / t)/dr
      k = sqrt(2. * ReduceMass / Hbar**2 * Echannel(i))

      H1(i, i) = H1(i, i) + (-t * exp(ai * k * dr))

    Else
      !k = acos(1. - abs(Echannel(i))/2. / t)/dr
      k = sqrt(2. * ReduceMass / Hbar**2 * abs(Echannel(i)))
      H1(i, i) = H1(i, i) + (-t * exp(k * dr))
    End If 

  End do

  Do i = 1, Nlevel*R_iterat
      Do j = 1, Nlevel*R_iterat
        H2(i,j) = H1(i,j)
      End Do
  End Do

  Do i = 1, Nlevel
    H2(Nlevel*R_iterat+i, Nlevel*R_iterat+i) = 2.0d0*t - E + eps(i) + V(R_max,L_i) 
    H2(Nlevel*R_iterat + i - nlevel, Nlevel*R_iterat + i) = -t
    H2(Nlevel*R_iterat + i, Nlevel*R_iterat + i  - nlevel) = -t
    H2(Nlevel*R_iterat + i, Nlevel*R_iterat + i  + nlevel) = -t
  End Do

  Do i = 1, Nlevel*R_iterat
      dbpsi(i, i) = (1.0d0, 0.0d0)
  End Do


  Do i = 1, Nlevel
    ec = E - CPOT(i,i,R_iterat)
    rho = acos(1. - ec/2. / t)/dr * (R_max)
    rho = sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max)

    eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
    Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
    dbpsi(Nlevel*R_iterat + i, Nlevel*R_iterat + i) = gcw(L_i) + ai*fcw(L_i)

    ec = E - CPOT(i,i,R_iterat+1)
    !rho = acos(1. - ec/2. / t)/dr  * (R_max + dr)
    rho = sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max + dr)
    eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
    Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
    dbpsi(Nlevel*R_iterat + i + Nlevel, Nlevel*R_iterat + i) = gcw(L_i) + ai*fcw(L_i)
  End Do

  do i = 1, Nx2
      dbpsi0(i) = dconjg(dbpsi(i, Nx1 - Nlevel + 1))
  end do

End Subroutine Prameterize_H


SUBROUTINE MATMUL_CPLX(H2, PSI, HPSI, M, N, P)
    IMPLICIT NONE
    INTEGER M, N, P
    COMPLEX*16 H2(M, N), PSI(N, P), HPSI(M, P)

    COMPLEX*16 ALPHA, BETA
    CHARACTER*1 TRANSA, TRANSB

    ALPHA = (1.0D0, 0.0D0)
    BETA  = (0.0D0, 0.0D0)
    TRANSA = 'N'
    TRANSB = 'N'
    CALL ZGEMM(TRANSA, TRANSB, M, P, N, ALPHA,H2,M,PSI, N,BETA, HPSI, M)
    RETURN
END

SUBROUTINE CSOLVE_EQUATION(N, HPSI, HPSI0, CC)
      IMPLICIT NONE
      INTEGER N, INFO
      COMPLEX*16 HPSI(N,N), HPSI0(N), CC(N)

      INTEGER IPIV(N)
      COMPLEX*16 A(N,N), B(N)
      INTEGER I, J

      DO 20 J = 1, N
          DO 10 I = 1, N
              A(I,J) = HPSI(I,J)
   10     CONTINUE
   20 CONTINUE

      DO 30 I = 1, N
          B(I) = -HPSI0(I)
   30 CONTINUE

      CALL ZGESV(N, 1, A, N, IPIV, B, N, INFO)

      DO 40 I = 1, N
          CC(I) = B(I)
   40 CONTINUE

      RETURN
END