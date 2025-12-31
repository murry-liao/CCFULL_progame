!'''
!@File    :   CCFULL solved by different methods 
!@Time    :   2026/01/01 
!@Author  :   Liao ZeHong, Kouichi Hagino
!@Version :   1.0
!@Contact :   liaozh26@mail2.sysu.edu.cn
!

Program CCFull_Main
  Use ccfull_initialization_mod
  Implicit None
  Integer ::  ie, il
  Integer ::  Iestep
  Real(8) ::  sigma, spin, s0, P_0, P, V
  Real(8) :: P_1, P_2, P_3, P_4
  Real(8) :: P0_1, P0_2, P0_3, P0_4
  Real(8) :: sigma_1 = 0.0D0, sigma_2 = 0.0D0, sigma_3 = 0.0D0, sigma_4 = 0.0D0
  Real(8) :: spin_1 = 0.0D0, spin_2 = 0.0D0, spin_3 = 0.0D0, spin_4 = 0.0D0

  external V
  Call Initialize_CCFull()
  
  Iestep = int((Emax - Emin)/dE)
  Do ie = 0, Iestep

    E = Emin + dE * ie
    sigma = 0.
    spin  = 0.
    P_0 = 0.

    sigma_1 = 0.0D0; sigma_2 = 0.0D0; sigma_3 = 0.0D0; sigma_4 = 0.0D0
    spin_1 = 0.0D0; spin_2 = 0.0D0; spin_3 = 0.0D0; spin_4 = 0.0D0

    Do il = 0, 200
      L_i = il
      !s0 = sigma
      s0 = sigma_1 ! use numerov method as a standard to cutoff

      Call PotShape(L_i, R_barrier_l, V_barrier_l, curv_l, R_bottom_l, V_bottom_l)
      If( V(R_bottom_l, L_i) > E .or. R_bottom_l < 0 )Then
          P = 0.
          exit  
      End If

      Call Numerov(P_1)
      Call DiscreteBasis(P_2)
      Call ModifiedDiscreteBasis(P_3)
      Call S_matrix(P_4)

      sigma_1 = sigma_1 + (2.0D0 * L_i + 1.0D0) * P_1 * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0
      sigma_2 = sigma_2 + (2.0D0 * L_i + 1.0D0) * P_2 * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0
      sigma_3 = sigma_3 + (2.0D0 * L_i + 1.0D0) * P_3 * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0
      sigma_4 = sigma_4 + (2.0D0 * L_i + 1.0D0) * P_4 * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0

      spin_1 = spin_1 + (2.0D0 * L_i + 1.0D0) * P_1 * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0 * L_i
      spin_2 = spin_2 + (2.0D0 * L_i + 1.0D0) * P_2 * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0 * L_i
      spin_3 = spin_3 + (2.0D0 * L_i + 1.0D0) * P_3 * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0 * L_i
      spin_4 = spin_4 + (2.0D0 * L_i + 1.0D0) * P_4 * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0 * L_i

      !sigma = sigma + (2.0D0 * L_i + 1.0D0) * P * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0
      !spin  = spin  + (2.0D0 * L_i + 1.0D0) * P * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0 * L_i

      !If (L_i == 0)P_0 = P
      If(L_i == 0)Then
        P0_1 = P_1; P0_2 = P_2
        P0_3 = P_3; P0_4 = P_4
      End If
      If(sigma_1 - s0 < s0*1.d-4)Exit 


      !Write(*,*)E, L_i, P

    End Do

    ! output the fusion cross section as a function of energy
    If (sigma_1 < 0.01D0)then
        write(*,'(A,F8.3,A,4(A,E15.5,A))') ' E = ', E, ' MeV,  ', &
                       'Num Sigma =', sigma_1, ' mb  ', 'Dbm Sigma =', sigma_2, ' mb  ', &
                       'Mdbm Sigma =', sigma_3, ' mb  ', 'Rmatrix Sigma =', sigma_4, ' mb  '
                 
        write(20,'(F15.5,4E15.5,4F15.5,4E15.5)') E, sigma_1, sigma_2, sigma_3, sigma_4, &
                        spin_1/sigma_1, spin_1/sigma_1, spin_1/sigma_1, spin_1/sigma_1, &
                        P0_1, P0_2, P0_3, P0_4
     
          
    else
        write(*,'(A,F8.3,A,4(A,F15.5,A))') ' E = ', E, ' MeV,  ', &
                       'Num Sigma =', sigma_1, ' mb  ', 'Dbm Sigma =', sigma_2, ' mb  ', &
                       'Mdbm Sigma =', sigma_3, ' mb  ', 'Rmatrix Sigma =', sigma_4, ' mb  '
                 
        write(20,'(F15.5,4F15.5,4F15.5,4E15.5)') E, sigma_1, sigma_2, sigma_3, sigma_4, &
                        spin_1/sigma_1, spin_1/sigma_1, spin_1/sigma_1, spin_1/sigma_1, &
                        P0_1, P0_2, P0_3, P0_4
    end if



    !if (sigma < 0.01D0) then
    !    write(*,'(A,F8.3,A,E15.5,A,F8.5,A,E15.5)') ' E = ', E, ' MeV,   Sigma = ', sigma, ' mb, <L> = ', spin / sigma, ' hbar, <P_0> = ', P_0
    !    write(20,'(F15.5,E15.5,F15.5,E15.5)') E,  sigma,  spin / sigma,  P_0
    !else
    !    write(*,'(A,F8.3,A,F15.5,A,F8.5,A,E15.5)') ' E = ', E, ' MeV,   Sigma = ', sigma, ' mb, <L> = ', spin / sigma, ' hbar, <P_0> = ', P_0
    !    write(20,'(4F15.5)') E,  sigma,  spin / sigma,  P_0
    !end if

    ! output the wave function


  End Do

  Print *, "Program completed successfully."

End Program CCFull_Main