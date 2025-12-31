Subroutine Sub_R_matrix(R_matrix)
  Use ccfull_initialization_mod
  Implicit None
  Complex*16, Intent(Out) :: R_matrix(Nlevel, Nlevel)
  real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  real(8), dimension(0:200) :: sigmad
  integer, dimension(0:200) :: iexp
  Integer ::  n, m, ch_i, ch_j, row, col, Ntot, i
  Real(8) ::  a, phi_i, phi_j, x_n, x_m
  Real(8) ::  V, r_point
  Real(8) ::  Delta_r

  complex*16 :: ai
  complex*16, Allocatable :: c(:,:), cinv(:,:), k_abs(:)
  complex*16 :: total_sum

  Complex*16 :: k_loc
  Real(8) :: E_local
  Real(8) :: phi_min_n, phi_min_m


  external V
  

  ai=(0.d0,1.d0)
  a = R_max
  Delta_r = R_max - R_min
  Ntot = Nlevel * Nbase
  Allocate(c(Ntot, Ntot))
  Allocate(cinv(Ntot, Ntot))
  Allocate(k_abs(Nlevel))
  c = 0.d0

  CPOT0 = 0.d0
  Call coupled_matrix(R_min, CPOT0)


  Do n = 1, Nbase
    x_n = xgauss(n)
    r_point = R_min + Delta_r * x_n
    CPOT0 = 0.d0
    Call coupled_matrix(r_point, CPOT0)
    Do m = 1, Nbase
      x_m = xgauss(m)

      phi_min_n = (-1.d0)**(n+1) * sqrt(xgauss(n)*(1.d0-xgauss(n))*Delta_r) / (xgauss(n)*Delta_r)
      phi_min_m = (-1.d0)**(m+1) * sqrt(xgauss(m)*(1.d0-xgauss(m))*Delta_r) / (xgauss(m)*Delta_r)

      Do ch_i = 1, Nlevel
        Do ch_j = 1, Nlevel

          ! the idex of the matrix element
          row = (ch_i - 1) * Nbase + n
          col = (ch_j - 1) * Nbase + m
          
          !Diagonal Blocks
          If (ch_i == ch_j) Then
            If (n == m) Then
              !kinetic energy 
              c(row,col) = Hbar**2 / (2.d0*ReduceMass) &
                         * (dble(Nbase)**2 + dble(Nbase) + 6.d0 - 2.d0/(x_n*(1.d0 - x_n))) &
                         / (3.d0 * Delta_r**2 * x_n * (1.d0 - x_n))

              c(row, col) = c(row, col) - E + V(R_min + Delta_r * x_n, L_i) ! - ai*W(R_max*x_n)
            Else
              !kinetic energy
              c(row,col) = Hbar**2 / (2.d0*ReduceMass)  &
                         * (-1.d0)**(n+m) / Delta_r**2 * Sqrt(x_m*(1.d0-x_m)/x_n/(1.d0-x_n))&
                         * (2.d0*x_n*x_m + 3*x_n - x_m - 4*x_n**2) &
                         / (x_n*(1.d0 - x_n)*(x_m - x_n)**2)

            End If
            !Bloch operator
            c(row, col) = c(row, col) + Hbar**2 / (2.d0*ReduceMass)* (-1.d0)**(n+m) / Delta_r**2 &
                        * Sqrt(x_n*x_m/(1.d0-x_n)/(1.d0-x_m)) &
                        * (dble(Nbase)*(dble(Nbase) + 1) - 1.d0/(1.d0-x_m))

            c(row, col) = c(row, col) - Hbar**2 / (2.d0*ReduceMass) * (-1.d0)**(n+m) / Delta_r**2 &
                        * Sqrt((1.d0-x_n)*(1.d0-x_m)/x_n/x_m) &
                        * (-dble(Nbase)*(dble(Nbase) + 1) + 1.d0/x_m)

          If (ch_i == ch_j) Then
             ! Energy on the R_min 
             E_local = E - V(R_min, L_i) - CPOT0(ch_i, ch_i)
             
             ! wave number on the R_min
             k_loc = Sqrt(DCMPLX(2.d0 * ReduceMass * E_local, 0.d0)) / Hbar
             
             ! -i * (hbar^2/2u) * k * phi_n * phi_m
             c(row, col) = c(row, col) - ai * (Hbar**2 / (2.d0*ReduceMass)) * k_loc * phi_min_n * phi_min_m
          End If



          End If

          !Off-Diagonal Blocks
          !(Coupling Potential)
          If (n == m) Then
          c(row, col) = c(row, col) + CPOT0(ch_i, ch_j)
          End If
          

        End Do
      End Do

    End Do
  End Do
    
    !Call matinv_R(Ntot, Ntot, c, cinv)
    Call matinv_lapack(Ntot, c, cinv)

  R_matrix = 0.d0
  Do ch_i = 1, Nlevel
    Do ch_j = 1, Nlevel
      
      Do n = 1, Nbase
        Do m = 1, Nbase
           row = (ch_i - 1) * Nbase + n
           col = (ch_j - 1) * Nbase + m
           
           phi_i = (-1)**(Nbase + n)* sqrt(xgauss(n)*(1-xgauss(n))*Delta_r)/(a-xgauss(n)*Delta_r - R_min)
           phi_j = (-1)**(Nbase + m)* sqrt(xgauss(m)*(1-xgauss(m))*Delta_r)/(a-xgauss(m)*Delta_r - R_min)

           
           R_matrix(ch_i, ch_j) = R_matrix(ch_i, ch_j) + &
                                      (Hbar**2 / (2.d0 * ReduceMass * a)) * &
                                      phi_i * cinv(row, col) * phi_j
        End Do
      End Do
      
    End Do
  End Do

  Deallocate(c, cinv)

End Subroutine Sub_R_matrix


Subroutine S_matrix(P)
  Use ccfull_initialization_mod
  
  Implicit None
  Real(8), Intent(Out) :: P  
  real(8) :: fcw_val, gcw_val, fpcw_val, gpcw_val
  integer :: iexp_val
  Integer :: i, j, ch, io
  Real(8) :: a, ak_ch, rho, eta_ch, E_ch, k, kk
  Real(8) :: qk_vec(Nlevel), V

  Complex*16 :: R_mat(Nlevel, Nlevel)
  Complex*16 :: S_mat(Nlevel, Nlevel)
  Complex*16 :: Z_in(Nlevel, Nlevel), Z_out(Nlevel, Nlevel)
  Complex*16 :: Z_out_inv(Nlevel, Nlevel)
  
  Complex*16 :: H_plus_vec(Nlevel), dH_plus_vec(Nlevel)   ! Outgoing
  Complex*16 :: H_minus_vec(Nlevel), dH_minus_vec(Nlevel) ! Incoming
  
  Complex*16 :: ai
  external V
  ai = (0.d0, 1.d0)
  a = R_max

  ! 1. R_matrix (Nlevel x Nlevel)
  Call Sub_R_matrix(R_mat)


  Call coupled_matrix(R_max, CPOT0)
  ! 2. Coulomb wavfunction (external function)
  Do ch = 1, Nlevel

    E_ch = E - CPOT0(ch,ch)
    !E_ch = E -  V(R_min, L_i) - CPOT0(ch,ch)
    ak_ch = sqrt(2.d0 * ReduceMass * E_ch / hbar**2)
    qk_vec(ch) = ak_ch
    rho = ak_ch * a
    eta_ch = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * E_ch))

    Call myCOULFG(L_i*1.d0, eta_ch, rho, &
                  fcw_val, fpcw_val, gcw_val, gpcw_val, iexp_val)
    
    ! H+ (Outgoing) and H- (Incoming)
    ! H = G + iF
    H_plus_vec(ch)  = gcw_val + ai * fcw_val
    dH_plus_vec(ch) = ak_ch * (gpcw_val + ai * fpcw_val) ! å¯¼æ•° d/dr
    
    H_minus_vec(ch)  = gcw_val - ai * fcw_val
    dH_minus_vec(ch) = ak_ch * (gpcw_val - ai * fpcw_val)
      

  End Do
  ! S = Z_out^{-1 * Z_in
  ! Z_ij = H_i * delta_ij - R_ij * a * dH_j
  Z_out = (0.d0, 0.d0)
  Z_in  = (0.d0, 0.d0)
  S_mat = (0.d0, 0.d0)

  Do i = 1, Nlevel
     Do j = 1, Nlevel
        
        Z_out(i, j) = - R_mat(i, j) * a * dH_plus_vec(j)
        Z_in(i, j)  = - R_mat(i, j) * a * dH_minus_vec(j)
        
        If (i == j) Then
           Z_out(i, j) = Z_out(i, j) + H_plus_vec(i)
           Z_in(i, j)  = Z_in(i, j) + H_minus_vec(i)
        End If
        
     End Do
  End Do

  ! 4. matrix inversion
  !Call matinv_R(Nlevel, Nlevel, Z_out, Z_out_inv)
  Call matinv_lapack(Nlevel, Z_out, Z_out_inv)
  ! S_mat = Z_out_inv * Z_in
  S_mat = Matmul(Z_out_inv, Z_in)

  ! 5. (Symmetric S-matrix)
  ! S_phys(i,j) = S(i,j) * sqrt(k_i / k_j)
  Do i = 1, Nlevel
     Do j = 1, Nlevel
        If (qk_vec(i) > 0.d0 .and. qk_vec(j) > 0.d0) Then
           S_mat(i, j) = S_mat(i, j) * Sqrt(qk_vec(i) / qk_vec(j))
        Else
           S_mat(i, j) = (0.d0, 0.d0)
        End If
     End Do
  End Do
  
  ! 6. 
  ! P = 1 - Sum(|S_1j|^2)
  P = 1.d0
  Do io = 1, Nlevel
      If (E - eps(io) .lt. 0.d0) Cycle
      k = sqrt((2.d0 * ReduceMass / Hbar**2 * (E - eps(io))))
      kk = sqrt(2.d0 * ReduceMass / Hbar**2 * E)
      P = P - Abs(S_mat(1, io))**2 

  End Do

  If (P < 0.d0) P = 0.d0

End Subroutine S_matrix




SUBROUTINE gaussh(n, x, w)
  !---------------------------------------------------------
  ! Purpose : Compute the zeros of Legendre polynomial Pn(x)
  !           in the interval [0,1] (mapped from [-1,1]),
  !           and the corresponding weighting coefficients
  !           for Gauss-Legendre integration.
  !
  ! Input :   n    --- Order of the Legendre polynomial (Number of points)
  ! Output:   x(n) --- Zeros of the Legendre polynomial (Mapped to [0,1])
  !           w(n) --- Corresponding weighting coefficients
  !
  ! Algorithm: Newton-Raphson Method
  ! [cite_start]Source: Adapted from rmatrix.f [cite: 31-36]
  !---------------------------------------------------------
  Implicit None

  ! Argument declarations
  Integer, Intent(In)  :: n
  Real(8), Intent(Out) :: x(n), w(n)

  ! Local variables
  Real(8) :: z, z0, p, f0, f1, pf, pd, fd, q, wp, gd
  Real(8), Parameter :: pi = 3.1415926535898D0
  Real(8), Parameter :: one = 1.0D0
  Integer :: n0, nr, i, k, j

  ! Determine the number of roots to compute (due to symmetry)
  n0 = (n + 1) / 2

  Do nr = 1, n0
      ! Initial guess using high-precision trigonometric approximation
      z = Cos(pi * (Dble(nr) - 0.25D0) / (Dble(n) + 0.5D0))

      ! Newton-Raphson iteration loop
      Do
          z0 = z
          p = 1.0D0
          
          ! Deflation: Remove roots that have already been found
          Do i = 1, nr - 1
              p = p * (z - x(i))
          End Do

          f0 = 1.0D0
          ! Handle the center root for odd N
          If (nr == n0 .And. Mod(n, 2) /= 0) Then
              z = 0.0D0
          End If
          f1 = z

          ! Compute Pn(z) and its derivative using recurrence relations
          Do k = 2, n
              pf = (2.0D0 - one / Dble(k)) * z * f1 - (1.0D0 - one / Dble(k)) * f0
              pd = Dble(k) * (f1 - z * pf) / (1.0D0 - z * z)
              f0 = f1
              f1 = pf
          End Do

          ! Exit if the root is exactly zero to avoid division errors
          If (z == 0.0D0) Exit

          fd = pf / p
          q = 0.0D0
          
          ! Compute correction terms based on deflation
          Do i = 1, nr - 1
              wp = 1.0D0
              Do j = 1, nr - 1
                  If (j /= i) Then
                      wp = wp * (z - x(j))
                  End If
              End Do
              q = q + wp
          End Do

          gd = (pd - q * fd) / p
          
          ! Update Z
          z = z - fd / gd

          ! Convergence check (using machine precision tolerance)
          If (Abs(z - z0) <= Abs(z) * 1.0D-15) Exit
      End Do

      ! Store roots and weights (utilizing symmetry about 0)
      x(nr) = z
      x(n + 1 - nr) = -z
      w(nr) = 2.0D0 / ((1.0D0 - z * z) * pd * pd)
      w(n + 1 - nr) = w(nr)

      !x(n0 + 1 - nr) = z
      !x(n + 1 - nr) = -z
      !w(n0 + 1 - nr) = 2.0D0 / ((1.0D0 - z * z) * pd * pd)
      !w(n + 1 - nr) = w(n0 + 1 - nr)
  End Do

  !---------------------------------------------------------
  ! Interval Mapping: [-1, 1] -> [0, 1]
  ! This step is specific to the Lagrange Mesh Method used in
  ![cite_start]! the R-matrix code[cite: 36].
  !---------------------------------------------------------
  x(1:n) = (1.0D0 + x(n:1:-1)) / 2.0D0
  w(1:n) = w(n:1:-1) / 2.0D0

END SUBROUTINE gaussh


subroutine matinv_lapack(n, A, Ainv)
  implicit none
  integer, intent(in) :: n
  complex*16, intent(in) :: A(n,n)
  complex*16, intent(out) :: Ainv(n,n)

  integer :: info, lwork
  integer, allocatable :: ipiv(:)
  complex*16, allocatable :: work(:)
  complex*16 :: tmp(n,n)

  tmp = A
  Ainv = (0.0d0, 0.0d0)

  allocate(ipiv(n))

  ! LU 
  call zgetrf(n, n, tmp, n, ipiv, info)
  if (info /= 0) then
    print *, "Error in zgetrf: info = ", info
    stop
  end if


  lwork = n * n
  allocate(work(lwork))


  call zgetri(n, tmp, n, ipiv, work, lwork, info)
  if (info /= 0) then
    print *, "Error in zgetri: info = ", info
    stop
  end if

  Ainv = tmp

  deallocate(ipiv)
  deallocate(work)
end subroutine matinv_lapack

subroutine matinv_R(nlmax, nmax, cc, d)
    implicit none
    integer, intent(in) :: nlmax, nmax
    complex*16, intent(in) :: cc(nlmax, nlmax)
    complex*16, intent(out)   :: d(nlmax, nlmax)

    complex*16 :: u(nlmax), v(nlmax)
    complex*16 :: a, b, t
    real*8     :: eps1, eps2
    integer :: n, m, k, l, j
    complex*16 :: c(nlmax,nlmax)
    ! thresholds
    eps1 = 1.0d-20
    eps2 = 1.0d-30

    ! --- Initialize d as identity matrix ---
    c = cc
    do m = 1, nmax
        do n = 1, nmax
            if (n == m) then
                d(n, m) = (1.0d0, 0.0d0)
            else
                d(n, m) = (0.0d0, 0.0d0)
            end if
        end do
    end do

    ! --- Gauss-Jordan elimination ---
    do n = 1, nmax
        t = c(n, n)

        ! If pivot is too small, try row swapping
        if (abs(t) < eps1) then
            j = n
            do
                j = j + 1
                if (j > nmax) then
                    print *, "matrix not invertible"
                    return
                end if

                t = c(n, j)

                if (abs(t) > eps2) then
                    ! swap row n and j
                    do k = 1, nmax
                        u(k) = c(n, k)
                        v(k) = d(n, k)

                        c(n, k) = c(j, k)
                        d(n, k) = d(j, k)

                        c(j, k) = u(k)
                        d(j, k) = v(k)
                    end do
                    exit
                end if
            end do
        end if

        ! --- Eliminate other rows ---
        do k = 1, nmax
            if (k == n) cycle
            a = c(k, n) / c(n, n)

            do l = 1, nmax
                c(k, l) = c(k, l) - a * c(n, l)
                d(k, l) = d(k, l) - a * d(n, l)
            end do
        end do

        ! --- Normalize the pivot row ---
        b = c(n, n)
        do m = 1, nmax
            c(n, m) = c(n, m) / b
            d(n, m) = d(n, m) / b
        end do

    end do

end subroutine matinv_R