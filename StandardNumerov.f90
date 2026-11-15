! StandardNumerov.f90
!
! The coupled-channels equations integrated with the standard Numerov step.
! With xi = (1 - T) u and T = dr**2 (2 mu / hbar**2) (V - E) / 12, V including the coupling matrix,
!     xi(n+1) = [12 (1 - T(n))**(-1) - 10] xi(n) - xi(n-1).
! Numerov.f90 uses the modified step of the original CCFULL, (sqrt(3) + sqrt(12) T)**2 - 1
! = 2 + 12 T + 12 T**2, which truncates (1 - T)**(-1) after T**2.  The standard step is linear
! in the potential.
! The solutions are stabilized at the barrier and every 0.5 fm (Numerov.f90: at the barrier only).
! The potentials, the boundary conditions, and the Coulomb matching are those of Numerov.f90;
! rkutta00 is the one in Numerov.f90 and matinv the one in Interaction.f90.

Subroutine StabilizeStandard(xi1, xi, aa, ir1, psi_out)
  Use ccfull_initialization_mod
  Implicit None
  Integer, Intent(In)  ::ir1
  complex*16, Intent(Out) :: psi_out(Nlevelmax, Nlevelmax)
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
  psi_out = (0.d0, 0.d0)
  psi_out(1:Nlevel, 1:Nlevel) = psi(1:Nlevel, 1:Nlevel)

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
End Subroutine StabilizeStandard

Subroutine StandardNumerov(P)
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
  ! Stabilization at the barrier and every 0.5 fm keeps the propagated solutions
  ! independent.  With the barrier alone, the highest partial waves lose it near
  ! some energies.
  complex*16, dimension(Nlevelmax, Nlevelmax) :: psi_stab, b0in, tmp
  complex*16, allocatable :: stab_psi(:,:,:), coef(:,:)
  integer :: nstab, nseg, kseg, maxseg
  Integer ::  i, j, lc, io, io2, i0, ir, ii, ik, ic, ich, ibarrier
  Integer ::  j1, j2
  Real(8) ::  V, fac, r, r1, r2
  Real(8) ::  ak, ec, eta, rho
  external V

  nstab = Max(1, Nint(0.5d0/dr))
  maxseg = R_iterat/nstab + 3
  Allocate(stab_psi(Nlevelmax,Nlevelmax,maxseg), coef(Nlevelmax,0:maxseg))
  nseg = 0

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

    ! Standard Numerov step in the xi = (1 - T) u representation, T = dr**2 Q / 12:
    !   xi(n+1) = [12 (1 - T(n))**(-1) - 10] xi(n) - xi(n-1).
    ! Numerov.f90 uses the modified step (sqrt(3) + dr**2 Q / sqrt(12))**2 - 1,
    ! which truncates (1 - T)**(-1) after T**2.  The standard step is linear in the potential.
    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd0(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, ir - 1)
            If (i0 == ic) Then
                dd0(i0, ic) = dd0(i0, ic) + 1.d0 - fac / 12.d0 * (V(r1, L_i) - E)
            End If
        End Do
    End Do

    Call matinv(Nlevel, dd0, dd1)
    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd1(i0, ic) = 12.d0 * dd1(i0, ic)
            If (i0 == ic) Then
                dd1(i0, ic) = dd1(i0, ic) - 10.d0
            End If
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
    If (ir == ibarrier .Or. Mod(ir, nstab) == 0) Then
      Call StabilizeStandard(xi1, xi, aa, ir, psi_stab)
      nseg = nseg + 1
      stab_psi(:, :, nseg) = psi_stab
    End If
    

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

  ! c(last) = bb0^(-1) e_1 in the last segment's basis; c(k-1) = psi_k^(-1) c(k).
  ! c(0) are the r_min-basis coefficients that the matrix bb (= bb0 aa) gave.
  tmp = bb0
  Call matinv(Nlevel, tmp, b0in)
  coef = (0.d0, 0.d0)
  coef(1:Nlevel, nseg) = b0in(1:Nlevel, 1)
  Do kseg = nseg, 1, -1
    tmp = stab_psi(:, :, kseg)
    Call matinv(Nlevel, tmp, cin)
    Do i = 1, Nlevel
      Do j = 1, Nlevel
        coef(i, kseg-1) = coef(i, kseg-1) + cin(i, j) * coef(j, kseg)
      End Do
    End Do
  End Do
  Do i = 1, Nlevel
    bin(i, 1) = coef(i, 0)
  End Do

  P = 0.d0

  Do io = 1, Nlevel
      If (ech(io) .lt. 0.d0) Cycle
      k = sqrt((2.d0 * ReduceMass / Hbar**2 * ech(io)))
      kk = sqrt(2.d0 * ReduceMass / Hbar**2 * e)
      p = p + (abs(bin(io,1)))**2 * k / kk
  End Do


  Deallocate(stab_psi, coef)

  Return

End Subroutine StandardNumerov
