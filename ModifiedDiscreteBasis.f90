Subroutine ModifiedDiscreteBasis(P)
  Use ccfull_initialization_mod

  Implicit None
  real(8), intent(out) :: P
  Real(8) ::  k, kk
  Integer ::  i, j
  complex*16 Hpsi(Nx1, Nx1)
  complex*16 Hpsi0(Nx1, 1)
  complex*16 cc(Nx1)

  Call Prameterize_H_numerove(Hamilton, Ec_psi, Ec_psi0)

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


Subroutine Prameterize_H_numerove(H2, dbpsi,dbpsi0)
  Use ccfull_initialization_mod
  implicit none
  !H2, dbpsi, dbpsi0
  complex*16, dimension(Nlevelmax, Nlevelmax) :: dd0, dd1
  complex*16, dimension(Nlevelmax, Nlevelmax) :: B, B2, C, cc, tt, ee, Hc
  complex*16 H1(Nx_cc,Nx_cc)
  complex*16 H2(Nx1, Nx2)
  complex*16 dbpsi(Nx2, Nx1)
  complex*16 dbpsi0(Nx2)
  complex*16 ai
  real(8), dimension(Nlevelmax, Nlevelmax, -1:R_iterat + 1) :: A

  Integer ::  i, j, ir, ie, il, ik, ic, i0, idx, jdx
  Real(8) ::  r,rh, V, eta, rho, k, ec, fac,temp
  !real(8), intent(in) :: para
  Real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  Real(8), dimension(0:200) :: sigmad
  Integer, dimension(0:200) :: iexp
  external V
  H1 = (0.d0, 0.d0)
  H2 = (0.d0, 0.d0)
  ai = (0.d0, 1.d0)
  fac = dr**2*(2.d0*ReduceMass/Hbar**2) 
  !BetaT = para
  ! contruct the cpot
  Call coupled_matrix0

  rh = R_min - dr 
  Call coupled_matrix(rh, CPOTH)
  Do  i = 1, Nlevel
    Do  j = 1, Nlevel
      A(i, j, -1) = 2*ReduceMass/Hbar**2 * CPOTH(i, j)
      if(i == j)Then
        A(i, j, -1) = A(i, j, -1) + 2*ReduceMass/Hbar**2 * (V(rh, L_i) - E)
      End If
    End Do
  End Do

  rh = R_min + dr / 2.0
  Call coupled_matrix(rh, CPOTH)


  Do ir = 0, R_iterat + 1
    r = R_min + dr * ir
    Call coupled_matrix(r, CPOT0)
    Do  i = 1, Nlevel
      Do  j = 1, Nlevel
        CPOT(i, j, ir) = CPOT0(i, j)
        A(i, j, ir) = 2*ReduceMass/Hbar**2 * CPOT0(i, j)
        if(i == j)Then
          A(i, j, ir) = A(i, j, ir) + 2*ReduceMass/Hbar**2 * (V(r, L_i) - E)
        End If
      End Do
    End Do
  End Do
  
  Do ir = 0, R_iterat + 1
    r = R_min + dr * ir
    
    Do  i = 1, Nlevel
      Do  j = 1, Nlevel
        B(i,j) = 0.0D0
        B(i,j) = dr**2 /SQRT(12.) * A(i, j, ir)
        if(i == j)B(i,j) = B(i,j) + sqrt(3.d0)
      End Do
    End Do

    Do i = 1, Nlevel
      Do j = 1, Nlevel
        B2(i,j) = 0.0D0
        Do ik = 1, Nlevel
          B2(i,j) = B2(i,j) + B(i,ik) * B(ik,j)
        End Do
        if(i == j)B2(i,j) = B2(i,j) - 1
      End Do
    End Do

    Do i = 1, Nlevel
      Do j = 1, Nlevel
          cc(i, j) = 0.0D0
          cc(i, j) = - dr**2 / 12.d0 * A(i, j, ir)
          If (i == j)cc(i, j) = cc(i, j) + 1.d0
      End Do
    End Do

    Do i = 1, Nlevel
      Do j = 1, Nlevel
        C(i,j) = 0.0D0
        Do ik = 1, Nlevel
          C(i,j) = C(i,j) + B2(i,ik) * cc(ik,j)
        End Do
      End Do
    End Do

    If(ir .lt. R_iterat)Then
      Do  i = 1, Nlevel
        Do  j = 1, Nlevel
          idx = Nlevel * ir + i
          jdx = Nlevel * ir + j
          H1(idx, jdx) = H1(idx, jdx) - C(i, j)
        End Do
      End Do
    End If

    If(ir == R_iterat)Then
      Do i = 1, Nlevel
        Do j = 1, Nlevel
          Hc(i,j) = -C(i,j)
        End Do
      End Do
    End If

  End Do


  Do ir = 0, (R_iterat-2)
    Do i = 1, Nlevel
      Do j = 1, Nlevel
          H1(ir * nlevel + i, ir * nlevel + nlevel + j) =  - dr**2 /12. * A(i, j, ir+1)
          H1(ir * nlevel + nlevel + i, ir * nlevel + j) =  - dr**2 /12. * A(i, j, ir)
          if(i == j)Then
            H1(ir * nlevel + i, ir * nlevel + nlevel + j) = H1(ir * nlevel + i, ir * nlevel + nlevel + j) + 1
            H1(ir * nlevel + nlevel + i, ir * nlevel + j) = H1(ir * nlevel + nlevel + i, ir * nlevel + j) + 1
          End If
      End Do
    End Do
  End Do
  
  
  Do i = 1, Nlevel
    Do j = 1, Nlevel
    Echannel(i) = E - V(R_min, L_i) - CPOT(i,i,0)
    ec = E - V(R_min, L_i) - CPOT(i,j,0)

    If(ec > 0)Then
      !k = acos(1. - ec/2. / t)/dr
      k = sqrt(2. * ReduceMass / Hbar**2 * ec)
      if(i == j)Then
        H1(i, j) = H1(i, j) + (1 - dr**2 /12. * A(i, j, -1)) * exp(ai * k * dr)
      else
        H1(i, j) = H1(i, j) + (  - dr**2 /12. * A(i, j, -1)) * exp(ai * k * dr)
      End If

    Else
      k = acos(1. - abs(ec)/2. / t)/dr
      H1(i, j) = H1(i, j) + ((- dr**2 /12. * A(i, j, -1)) * exp(k * dr))
      if(i == j)H1(i, j) = H1(i, j) + 1 * exp(k * dr)

    End If 

    End do
  End do

  Do i = 1, Nlevel*R_iterat
      Do j = 1, Nlevel*R_iterat
        H2(i,j) = H1(i,j)
      End Do
  End Do

  Do i = 1, Nlevel
    Do j = 1, Nlevel

      H2(Nlevel*R_iterat + i, Nlevel*R_iterat + j) = Hc(i, j)
      H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) = (- dr**2 /12. * A(i, j, R_iterat))
      H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) = (- dr**2 /12. * A(i, j, R_iterat-1))
      H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) = (- dr**2 /12. * A(i, j, R_iterat+1))

      If(i == j)Then
      H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) = H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) + 1
      H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) = H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) + 1
      H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) = H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) + 1
      End If

    End Do
  End Do

  Do i = 1, Nlevel*R_iterat
      dbpsi(i, i) = (1.0d0, 0.0d0)
  End Do


  Do i = 1, Nlevel
    ec = E - CPOT(i,i,R_iterat)
    !rho = acos(1. - ec/2. / t)/dr * (R_max)
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

End Subroutine Prameterize_H_numerove