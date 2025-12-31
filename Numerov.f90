Subroutine rkutta00(psi0, phi0, psi1)
  Use ccfull_initialization_mod
  Implicit None
  complex*16, Intent(In)  ::  psi0(Nlevelmax), phi0(Nlevelmax)
  complex*16, intent(out) ::  psi1(Nlevelmax)
  complex*16 :: ai
  complex*16, dimension(Nlevelmax) :: ak1, ak2, ak3, ak4
  complex*16, dimension(Nlevelmax) :: bk1, bk2, bk3, bk4

  Real(8) ::  fac, r, rh, rpp, V
  Integer ::  j1, i0, ic, is
  external V

  ai = (0.d0, 1.d0)
  fac = dr * (2.d0 * ReduceMass / Hbar**2)

  r  = R_min
  rh = R_min + dr/2.d0
  rpp= R_min + dr
  Do j1 = 1, Nlevel
      ak1(j1) = 0.0D0
      ak2(j1) = 0.0D0
      ak3(j1) = 0.0D0
      ak4(j1) = 0.0D0
      bk1(j1) = 0.0D0
      bk2(j1) = 0.0D0
      bk3(j1) = 0.0D0
      bk4(j1) = 0.0D0
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak1(i0) = ak1(i0) + fac * CPOT(i0, ic, 0) * psi0(ic)
      End Do
      ak1(i0) = ak1(i0) - fac * (E - V(r, L_i)) * psi0(i0)
      bk1(i0) = dr * phi0(i0)
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak2(i0) = ak2(i0) + fac * CPOTH(i0, ic) * (psi0(ic) + 1.0D0 / 2.0D0 * bk1(ic))
      End Do
      ak2(i0) = ak2(i0) - fac * (E - V(rh, L_i)) * (psi0(i0) + 1.0D0 / 2.0D0 * bk1(i0))
      bk2(i0) = dr * (phi0(i0) + 1.0D0 / 2.0D0 * ak1(i0))
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak3(i0) = ak3(i0) + fac * CPOTH(i0, ic) * (psi0(ic) + 1.0D0 / 2.0D0 * bk2(ic))
      End Do
      ak3(i0) = ak3(i0) - fac * (E - V(rh, L_i)) * (psi0(i0) + 1.0D0 / 2.0D0 * bk2(i0))
      bk3(i0) = dr * (phi0(i0) + 1.0D0 / 2.0D0 * ak2(i0))
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak4(i0) = ak4(i0) + fac * CPOT(i0, ic, 1) * (psi0(ic) + bk3(ic))
      End Do
      ak4(i0) = ak4(i0) - fac * (E - V(rpp, L_i)) * (psi0(i0) + bk3(i0))
      bk4(i0) = dr * (phi0(i0) + ak3(i0))
  End Do

  Do is = 1, nlevel
      psi1(is) = psi0(is) + (1.0D0 / 6.0D0) * (bk1(is) + 2.0D0 * bk2(is) + 2.0D0 * bk3(is) + bk4(is))
  End Do

  Return

End Subroutine rkutta00

Subroutine stabilize(xi1, xi, aa, ir1)
  Use ccfull_initialization_mod
  Implicit None
  Integer, Intent(In)  ::ir1
  complex*16, Intent(Inout)  ::  aa(Nlevelmax, Nlevelmax)
  complex*16, Intent(Inout)  ::  xi(Nlevelmax, Nlevelmax)
  complex*16, Intent(Inout)  ::  xi1(Nlevelmax, Nlevelmax)
  complex*16 :: ai
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: psi
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: cc, cin
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: aa0
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: xid, xid1
  Real(8) ::  r, V, fac
  Integer ::  i, j, k, i0, ic, ich
  external  V

  Do i = 1, nlevel
      Do j = 1, nlevel
        aa0(i, j)  = aa(i, j)
        xid(i, j)  = xi(i, j)
        xid1(i, j) = xi1(i, j)
      End Do
  End Do

  ai = (0.d0, 1.d0)
  fac = dr**2 * (2.d0 * ReduceMass / Hbar**2)
  r = R_min + ir1*dr

  Do i0 = 1, Nlevel
    Do ic = 1, Nlevel
        cc(i0, ic) = -fac / 12.0D0 * CPOT(i0, ic, ir1)
        If (i0 == ic) Then
          cc(i0, ic) = cc(i0, ic) - fac / 12.0D0 * (V(r, L_i) - E) + 1.0D0
        End If
    End Do
  End Do

  Call matinv(Nlevel,cc,cin)

  Do ich = 1, nlevel
      Do i0 = 1, nlevel
        psi(i0, ich) = 0.0D0
        Do ic = 1, nlevel
            psi(i0, ich) = psi(i0, ich) + cin(i0, ic) * xi(ic, ich)
        End Do
      End Do
  End Do

  Do i = 1, nlevel
      Do j = 1, nlevel
        aa(i, j) = 0.0D0
        Do k = 1, nlevel
            aa(i, j) = aa(i, j) + psi(i, k) * aa0(k, j)
        End Do
      End Do
  End Do

  call matinv(Nlevel,psi,cin)

  Do i = 1, Nlevel
      Do j = 1, Nlevel
        xi(i, j)  = 0.0D0
        xi1(i, j) = 0.0D0
        Do k = 1, Nlevel
            xi(i, j)  = xi(i, j) + xid(i, k)  * cin(k, j)
            xi1(i, j) = xi1(i, j) + xid1(i, k) * cin(k, j)
        End Do
      End Do
  End Do

  Return
End Subroutine stabilize

Subroutine Numerov(P)
  Use ccfull_initialization_mod
  Implicit None
  real(8), intent(out) :: P
  complex*16, dimension(Nlevelmax) :: psi, psi0, psi1
  complex*16, dimension(Nlevelmax, Nlevelmax) :: xi, xi0, xi1
  complex*16, dimension(Nlevelmax) :: phi0
  complex*16, dimension(Nlevelmax, Nlevelmax) :: bb, bin
  complex*16, dimension(Nlevelmax, Nlevelmax) :: bb2
  complex*16, dimension(Nlevelmax, Nlevelmax) :: cc, cin
  complex*16, dimension(Nlevelmax, Nlevelmax) :: dd0, dd1
  complex*16, dimension(Nlevelmax) :: dd
  real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  real(8), dimension(0:200) :: sigmad
  integer, dimension(0:200) :: iexp
  real(8), dimension(Nlevelmax) :: ech, ech2
  complex*16 :: k, kk, k2
  complex*16 :: ai
  complex*16 :: cwup0, cwdown0, cwup1, cwdown1
  complex*16 :: dummy, dummy2
  complex*16, dimension(Nlevelmax, Nlevelmax) :: xi1d, xi0d
  complex*16, dimension(Nlevelmax, Nlevelmax) :: aa, bb0, bb20
  Integer ::  i, j, lc, io, io2, i0, ir, ii, ik, ic, ich, ibarrier
  Integer ::  j1, j2
  Real(8) ::  V, fac, r, r1, r2
  Real(8) ::  ak, ec, eta, rho
  external V

  ai = (0.d0, 1.d0)
  fac = dr**2 * (2.d0 * ReduceMass / Hbar**2)
  ibarrier = (R_barrier - R_min)/dr
  Do lc = 0, 200
    fcw(lc)   = 0.0d0
    gcw(lc)   = 0.0d0
    fpcw(lc)  = 0.0d0
    gpcw(lc)  = 0.0d0
    sigmad(lc) = 0.0d0
    iexp(lc)   = 0
  End Do

  Do i = 1, Nlevelmax
      Do j = 1, Nlevelmax
          bb(i,j)  = 0.0d0
          bin(i,j) = 0.0d0
          cc(i,j)  = 0.0d0
          cin(i,j) = 0.0d0
          aa(i,j)  = 0.0d0
      End Do
      dd(i) = 0.0d0
  End Do

  Do i = 1, Nlevel
      aa(i,i) = 1.0d0
  End Do


  Do io = 1, Nlevel

    Do j1 = 1, Nlevel
        psi(j1)  = 0.0d0
        psi0(j1) = 0.0d0
        psi1(j1) = 0.0d0
        phi0(j1) = 0.0d0
    End Do

    If (io == 1) Then
        Do io2 = 1, Nlevel
            ech(io2) = E - V(R_min, L_i) - CPOT(io2, io2, 0)
            ech2(io2) = E - CPOT(io2, io2, R_iterat)

        End Do
    End If

    If (ech(io) > 0.d0) Then
        k = sqrt(2.d0 * ReduceMass / Hbar**2 * ech(io))
        psi0(io) = exp(-ai * k * R_min)
        phi0(io) = -ai * k * psi0(io)

    Else
        k = sqrt(2.d0 * ReduceMass / hbar**2 * abs(ech(io)))
        psi0(io) = exp(k * R_min)
        phi0(io) = k * psi0(io)
    End If


    call rkutta00(psi0,phi0,psi1)

    do i0 = 1, Nlevel
        xi0(i0, io) = (1.d0 - fac/12.d0 * (V(R_min, L_i) - e)) * psi0(i0)
        xi1(i0, io) = (1.d0 - fac/12.d0 * (V(R_min+dr, L_i) - e)) * psi1(i0)
        do ic = 1, Nlevel
            xi0(i0, io) = xi0(i0, io) - fac/12.d0 * CPOT(i0, ic, 0) * psi0(ic)
            xi1(i0, io) = xi1(i0, io) - fac/12.d0 * CPOT(i0, ic, 1) * psi1(ic)
        end do
    end do
  End Do

  !----------------------------------------------  iterations start
  Do ir = 2, R_iterat + 1
    r = R_min + dr * ir
    !r0 = R_min + dr * (ir - 2)
    r1 = R_min + dr * (ir - 1)

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd0(i0, ic) = fac / sqrt(12.d0) * CPOT(i0, ic, ir - 1)
            If (i0 == ic) Then
                dd0(i0, ic) = dd0(i0, ic) + fac / sqrt(12.d0) * (V(r1, L_i) - E) + sqrt(3.d0)
            End If
        End Do
    End Do

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd1(i0, ic) = 0.d0
            If (i0 == ic) Then
                dd1(i0, ic) = dd1(i0, ic) - 1.d0
            End If
            Do ik = 1, Nlevel
                dd1(i0, ic) = dd1(i0, ic) + dd0(i0, ik) * dd0(ik, ic)
            End Do
        End Do
    End Do

    Do ich = 1, Nlevel
        Do i0 = 1, Nlevel
            xi(i0, ich) = -xi0(i0, ich)
            Do ic = 1, nlevel
                xi(i0, ich) = xi(i0, ich) + dd1(i0, ic) * xi1(ic, ich)
            End Do
        End Do
    End Do

    If(ir == R_iterat+1)Exit
    If(ir == ibarrier) Call stabilize(xi1,xi,aa,ir)
    

    Do ich = 1, nlevel
        Do i0 = 1, nlevel
            xi0(i0, ich) = xi1(i0, ich)
            xi1(i0, ich) = xi(i0, ich)
        End Do
    End Do


  End Do


  !--------------------------------------------------------------
  !  matching to the coulomb wave function at rmax

  Do io = 1, Nlevel

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, R_iterat - 1)
            If (i0 == ic) Then
                cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_max - dr, L_i) - E) + 1.d0
            End If
        End Do
    End Do


    Call matinv(Nlevel,cc,cin)
    Do i0 = 1, Nlevel
        psi0(i0) = 0.d0
        Do ic = 1, Nlevel
            psi0(i0) = psi0(i0) + cin(i0, ic) * xi0(ic, io)
        End Do
    End Do


    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, R_iterat + 1)
            If (i0 == ic) Then
                cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_max + dr, L_i) - E) + 1.d0
            End If
        End Do
    End Do

    Call matinv(Nlevel,cc,cin)

    Do i0 = 1, Nlevel
        psi(i0) = 0.d0
        Do ic = 1, Nlevel
            psi(i0) = psi(i0) + cin(i0, ic) * xi(ic, io)
        End Do
    End Do


    Do ii = 1, Nlevel
        ! coulomb wave function
        ec = E - CPOT(ii, ii, R_iterat - 1)

        If (ec < 0.d0) Then
            r1 = R_max - dr
            r2 = R_max + dr
            ak = sqrt(2.d0 * ReduceMass * abs(ec) / Hbar**2)
            
            bb0(ii, io) = (exp(-ak * r2) * psi0(ii) - exp(-ak * r1) * psi(ii)) &
                          / (exp(ak * (r1 - r2)) - exp(-ak * (r1 - r2)))
            bb20(ii, io) = -(exp(ak * r2) * psi0(ii) - exp(ak * r1) * psi(ii)) &
                          / (exp(ak * (r1 - r2)) - exp(-ak * (r1 - r2)))
        Else
            rho = sqrt(2.d0 * ReduceMass * ec) / hbar * (R_max - dr)
            eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))


            !Call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, L_i, iexp)
            Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
  
 
            cwup0 = (gcw(L_i) + ai * fcw(L_i))
            cwdown0 = gcw(L_i) - ai * fcw(L_i)

            ec = E - CPOT(ii, ii, R_iterat + 1)
            rho = Sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max + dr)
            eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
            !Call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, L_i, iexp)
            Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
            cwup1 = (gcw(L_i) + ai * fcw(L_i))
            cwdown1 = gcw(L_i) - ai * fcw(L_i)

            bb0(ii, io) = (cwup0 * psi(ii) - cwup1 * psi0(ii)) &
                          / (cwup0 * cwdown1 - cwup1 * cwdown0)
            bb20(ii, io) = (cwdown1 * psi0(ii) - cwdown0 * psi(ii)) &
                          / (cwup0 * cwdown1 - cwup1 * cwdown0)
        End If
    End Do

  End Do
  !===============================================================
  !                                        penetration probability

  Do i = 1, Nlevel
      Do j = 1, Nlevel
        bb(i,j) = 0.d0
        bb2(i,j) = 0.d0
        Do j2 = 1, nlevel
            bb(i,j) = bb(i,j) + bb0(i,j2) * aa(j2,j)
            bb2(i,j) = bb2(i,j) + bb20(i,j2) * aa(j2,j)
        End Do
      End Do
  End Do

  Call matinv(Nlevel,bb,bin)

  P = 0.d0

  Do io = 1, Nlevel
      If (ech(io) .lt. 0.d0) Cycle
      k = sqrt((2.d0 * ReduceMass / Hbar**2 * ech(io)))
      kk = sqrt(2.d0 * ReduceMass / Hbar**2 * e)
      p = p + (abs(bin(io,1)))**2 * k / kk
  End Do

  Return

End Subroutine Numerov