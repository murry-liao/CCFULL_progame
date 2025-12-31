Module ccfull_initialization_mod
    Implicit None   
  !  """Principal global variables
  !  P                   - penetrability 
  !  Sigma               - fusion cross section, unit mb
  !  Spin                - mean angular momentum
  !    
  !  Apro, Zpro, Rpro    - atomic #, proton # and radius of the projectile
  !  Atar, Ztar, Rtar    - those of the target
  !  ReduceMass          - reduced mass, unit MeV/c**2
  !  Amu                 - nucleon mass, unit MeV/c**2
  !  E                   - bombarding energy in the center of mass frame, unit MeV
  !    
  !  V0,R0,A0            - depth, range, and dissuseness parameters of uncoupled 
  !                        nuclear potential, which is assumed to be a Woods-Saxon form
  !    
  !  IVIBROTT (IVIBROTP) - option for intrinsic degree of freedom 
  !                        = -1; no excitation (inert)
  !                        =  0; vibrational coupling 
  !                        =  1; rotational coupling
  !
  !  Ntar (Npro)         - the number of levels to be included    
  !
  !  BetaTn (BetaPn)     - i am not very clear about it (Liao)
  !  BetaT (BetaP)       - defeormation parameter
  !  OmegaT (OmegaP)     - excitation energy of the oscillator
  !  LambdaT (LambdaP)   - multipolarity
  !  NphononT (NphononP) - the number of phonon to be included
  !  Beta2T (Beta2P)     - static quadrupole deformation parameter
  !  Beta4T (Beta4P)     - static hexadecapole deformation parameter
  !  BetaT2n             - is the same as BetaTn but for the second mode
  !  BetaT2              - is the same as BETAT but for the second mode
  !  OmegaT2             - is the same as OmegaT but for the second mode
  !  LambdaT2            - is the same as LambdaT but for the second mode
  !  NphononT2           - is the same as NphononT but for the second mode
  !
  !  E2T (E2P)           - excitation energy of 2+ state in a rotational band
  !  NrotT (NrotP)       - the number of levels in the rotational band to be included 
  !                        (up to I^pi=2*NROT+ states are included)
  !
  !  L                   - angular momentum of the relative motion
  !  CPOT                - coupling matrix
  !  """
  Real(8), parameter :: Hbar = 197.329d0
  Real(8), parameter :: Pi   = 3.141592653d0
  Real(8), parameter :: Amu  = 938.0
  Real(8), parameter :: e2   = 1.44

  integer, parameter :: Nlevelmax = 30
  integer, parameter :: Nbase = 120
  integer, Public ::  Apro, Zpro, Atar, Ztar, L_i
  Real(8), Public ::  R0P, R0T, Rpro, Rtar, ReduceMass, h2m
  Real(8), Public ::  OmegaT, BetaT, OmegaT2, BetaT2, E2T, Beta2T, Beta4T
  Real(8), Public ::  OmegaP, BetaP, E2P, Beta2P, Beta4P
  Real(8), Public ::  V0, R0, A0, Emin, Emax, dE, R_max, R_min, dr
  Real(8), Public ::  BetaTn, BetaPn, BetaT2n, BetaP2n
  Real(8), Public ::  R_barrier, V_barrier, curv, R_bottom, V_bottom
  Real(8), Public ::  R_barrier_l, V_barrier_l, curv_l, R_bottom_l, V_bottom_l
  Real(8), Public ::  t
  Real,    Public ::  Ftr, Qtrans
  Real,    Public ::  E
  Integer, Public ::  IVIBROTP, IVIBROTT, LambdaT, NphononT, LambdaT2, NphononT2
  Integer, Public ::  NrotT, Ntar, LambdaP, NphononP, NrotP, Npro
  Integer, Public ::  Ntrans, iq
  Integer, Public ::  Nlevel, R_iterat
  Real(8), Allocatable, Public :: Pot_para1(:), Pot_para2(:), Pot_para3(:)
  Real(8), Allocatable, Public :: Ev(:), eps(:)
  Real(8), Allocatable, Public :: omeahv(:), omeahv2(:), omeahvp(:), erott(:), erotp(:)
  Real(8), Allocatable, Public :: betnahv(:,:), betcahv(:,:), betnahv2(:,:), betcahv2(:,:)
  Real(8), Allocatable, Public :: betnahvp(:,:), betcahvp(:,:), bett(:,:), betp(:,:)
  Real(8), Allocatable, Public :: Evec(:,:), CPOT0(:,:), CPOTH(:,:), CPOT(:,:,:)
  Real(8), Allocatable, Public :: H(:,:), unit(:,:)
  Real(8), Allocatable, Public :: wgauss(:), xgauss(:)

  Integer, Allocatable, Public :: imutual(:,:)


  Integer, Public :: Nx, Nx_cc, Nx1, Nx2
  Real(8),    Allocatable, Public :: Echannel(:)
  complex*16, Allocatable, Public :: wf(:)
  complex*16, Allocatable, Public :: wf_base(:,:), wf_base2(:,:), wf_base3_dagger(:,:)
  complex*16, Allocatable, Public :: Hamilton(:,:)
  complex*16, Allocatable, Public :: Ec_psi(:,:), Ec_psi0(:), Ec_psi0_2(:)
  complex*16, Allocatable, Public :: cc_ec(:,:), cc_ec_inv(:,:),dd_ec(:), dd2_ec(:), a_tra(:)

  character(len=2), dimension(0:109) :: Element

  data Element / &
       'O ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
       'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ', &
       'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
       'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In', &
       'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm', &
       'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re', &
       'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra', &
       'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md', &
       'No','Lr','XX','X1','X2','X3','X4','04' /

Contains

  Subroutine Read_Input()

    Open(10, File='ccfull.inp', Status = 'Unknown')
    Open(20, File='sigma.dat',  Status = 'Unknown')
    write(20,'(13A15)') 'E',  'sigma_num', 'sigma_dbm', 'sigma_mdbm', 'sigma_rmatrix',  &
                      '<L>(hbar)',  '<L>(hbar)', '<L>(hbar)', '<L>(hbar)','P (L = 0)', 'P (L = 0)', 'P (L = 0)', 'P (L = 0)'
    Open(30, File='Output.dat', Status = 'Unknown')
    Read(10,*)Apro, Zpro, Atar, Ztar
    Read(10,*)R0P, IVIBROTP, R0T, IVIBROTT

    Rpro = R0P * Apro ** (1.d0 / 3.d0)
    Rtar = R0T * Atar ** (1.d0 / 3.d0)
    ReduceMass = Apro * Atar * 1.d0 / (Apro + Atar) * Amu

    h2m = Hbar**2 / (2 * ReduceMass)
    
    If (IVIBROTT == 0) Then

      Read(10,*)OmegaT, BetaT, LambdaT, NphononT
      Read(10,*)OmegaT2, BetaT2, LambdaT2, NphononT2
      Ntar = NphononT

    Else If (IVIBROTT == 1) Then

      Read(10,*)E2T, Beta2T, Beta4T, NrotT 
      Read(10,*)
      Ntar = NrotT
    Else
      Read(10,*)
      Read(10,*)
      Ntar = 0
    End If

    If (IVIBROTP == 0) Then

      Read(10,*)OmegaP, BetaP, LambdaP, NphononP
      Npro = NphononP

    Else If (IVIBROTP == 1) Then

      Read(10,*)E2P, Beta2P, Beta4P, NrotP
      Npro = NrotP
    Else
      Read(10,*)
      Npro = 0
    End If
    Nlevel=(Npro + 1)*(Ntar + 1)

    Read(10,*)Qtrans, Ftr
    Read(10,*)V0, R0, A0
    Read(10,*)Emin, Emax, dE
    Read(10,*)R_max, dr
  End Subroutine Read_Input

  Subroutine Output_information()
      Character(len=1) :: ans


      print '(A, I0 ,A, A,I0 ,A)', 'System: ', Int(Apro), Element(Zpro),'+', Int(Atar), Element(Ztar)

      print '(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)', &
            "Simulation range E from ", Emin, " MeV To ", Emax, " MeV, dE = ", dE, " MeV"

      print '(A, F8.4, A, F8.4, A, F8.3, A, F8.3, A)', &
            "Simulation range R from ", R_min, " fm To ", R_max, " fm, dR = ", dr, " fm"

      print '(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)', &
            "Potential parameters: V0= ", V0, "(MeV), a= ", A0, "(fm), r0= ", R0, "(fm), R= ", R0*(Apro**(1.0/3.0) + Atar**(1.0/3.0)), "(fm)"

      print '(A, F8.4, A)', "Coulomb barrier position :", R_barrier, " fm"

      print '(A, F8.4, A)', "Coulomb barrier energy   :", V_barrier, " MeV"

      print '(A, F8.4, A)', "Coulomb barrier curv     :", curv, " MeV"

      print '(A, F8.4, A)', "Coulomb bottom position  :", R_bottom, " fm"

      print '(A, F10.4, A)', "Coulomb bottom energy    :", V_bottom, " MeV"

      write(30,'(A, I0 ,A, A,I0 ,A)') 'System: ', Int(Apro), Element(Zpro),'+', Int(Atar), Element(Ztar)

      write(30,'(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)') &
            "Simulation range E from ", Emin, " MeV To ", Emax, " MeV, dE = ", dE, " MeV"

      write(30,'(A, F8.4, A, F8.4, A, F8.3, A, F8.3, A)') &
            "Simulation range R from ", R_min, " fm To ", R_max, " fm, dR = ", dr, " fm"

      write(30,'(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)') &
            "Potential parameters: V0= ", V0, "(MeV), a= ", A0, "(fm), r0= ", R0, "(fm), R= ", R0*(Apro**(1.0/3.0) + Atar**(1.0/3.0)), "(fm)"

      write(30,'(A, F8.4, A)') "Coulomb barrier position :", R_barrier, " fm"

      write(30,'(A, F8.4, A)') "Coulomb barrier energy   :", V_barrier, " MeV"

      write(30,'(A, F8.4, A)') "Coulomb barrier curv     :", curv, " MeV"

      write(30,'(A, F8.4, A)') "Coulomb bottom position  :", R_bottom, " fm"

      write(30,'(A, F10.4, A)') "Coulomb bottom energy    :", V_bottom, " MeV"

      print *, "------------------------------"
      print *, "Mode of excitation for taget"

      if (Ntar /= 0) then
          if (IVIBROTT == 0) then
              write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the targ.: beta=", BetaT, ", omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
              BetaTn = BetaT
              Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
              Read(*,*)ans
              IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
                Write(*,*)'beta_N=?'
                Read(*,*)BetaTn
              End If
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaTn, ", beta_C=", BetaT, ", r0=", R0T, " (fm),"
              write(*,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaTn, ", beta_C=", BetaT, ", r0=", R0T, " (fm),"
              write(30,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
          else if (IVIBROTT == 1) then
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the targ.: beta2=", Beta2T, ", beta4=", Beta4T, ", r0=", R0T, " (fm),"
              write(*,'(A, F6.3, A, I0)') "                                   E2=", E2T, " (MeV), Nrot=", NrotT
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the targ.: beta2=", Beta2T, ", beta4=", Beta4T, ", r0=", R0T, " (fm),"
              write(30,'(A, F6.3, A, I0)') "                                   E2=", E2T, " (MeV), Nrot=", NrotT
          end if
      end if

      if (NphononT2 /= 0) then
          write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the targ.: beta=", BetaT2, ", omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          BetaT2n = BetaT2
          Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
          Read(*,*)ans
          IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
            Write(*,*)'beta_N=?'
            Read(*,*)BetaT2n
          End If
          write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaT2n, ", beta_C=", BetaT2, ", r0=", R0T, " (fm),"
          write(*,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaT2n, ", beta_C=", BetaT2, ", r0=", R0T, " (fm),"
          write(30,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          Call Mutual
      End if

      print *, "------------------------------"
      print *, "Mode of excitation for projectile"

      if (Npro /= 0) then
          if (IVIBROTP == 0) then
              write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the proj.: beta=", BetaP, ", omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
              BetaPn = BetaP
              Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
              Read(*,*)ans
              IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
                Write(*,*)'beta_N=?'
                Read(*,*)BetaPn
              End If
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the proj.: beta_N=", BetaPn, ", beta_C=", BetaP, ", r0=", R0P, " (fm),"
              write(*,'(A, F6.3, A, F6.3, A, I0)') "                              omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the proj.: beta_N=", BetaPn, ", beta_C=", BetaP, ", r0=", R0P, " (fm),"
              write(30,'(A, F6.3, A, F6.3, A, I0)') "                              omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
          else if (IVIBROTP == 1) then
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the proj.: beta2=", Beta2P, ", beta4=", Beta4P, ", r0=", R0P, " (fm),"
              write(*,'(A, F6.3, A, I0)') "                                   E2=", E2P, " (MeV), Nrot=", NrotP
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the proj.: beta2=", Beta2P, ", beta4=", Beta4P, ", r0=", R0P, " (fm),"
              write(30,'(A, F6.3, A, I0)') "                                   E2=", E2P, " (MeV), Nrot=", NrotP
          end if
      end if
      Call Anharmonicity
      Call grotation

      print *, "------------------------------"
      print *, "Transfer coupled mode:"
      Write(30,*)"------------------------------"
      Write(30,*)"Transfer coupled mode:"
      If (Ntrans == 0)then
        Write(*,*)'No transfer coupled mode'
        Write(30,*)'No transfer coupled mode'
      End if 

      if (Ntrans /= 0) then
          write(*,'(A, F6.3, A, F6.3, A)') "Transfer channel: Strength= ", Ftr, ", Q = ", Qtrans, " MeV"
          write(30,'(A, F6.3, A, F6.3, A)') "Transfer channel: Strength= ", Ftr, ", Q = ", Qtrans, " MeV"
      end if
      print *, "------------------------------"
      Write(*,*)'Nlevle = ', Nlevel
      Write(30,*)'Nlevle = ', Nlevel
  End Subroutine Output_information


  Subroutine Initialize_CCFull()

    Implicit None
    Integer ::  i, j, ir
    Real(8) ::  r, rh

    Allocate(imutual(0:Nlevelmax,0:Nlevelmax))
    Allocate(Ev(Nlevelmax))
    Allocate(Evec(Nlevelmax,Nlevelmax))
    Allocate(betnahv(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahv(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahv(0:Nlevelmax))
    Allocate(betnahv2(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahv2(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahv2(0:Nlevelmax+1))
    Allocate(betnahvp(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahvp(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahvp(0:Nlevelmax))
    Allocate(bett(0:Nlevelmax,0:Nlevelmax))
    Allocate(erott(0:Nlevelmax+1))
    Allocate(betp(0:Nlevelmax,0:Nlevelmax))
    Allocate(erotp(0:Nlevelmax))
    Allocate(eps(0:Nlevelmax))
    Allocate(xgauss(Nbase), wgauss(Nbase))

    Call Read_Input()

    Call PotShape(0, R_barrier, V_barrier, curv, R_bottom, V_bottom)
    R_iterat = int((R_max - R_bottom) / dr)
    R_min = R_bottom
    R_max = R_bottom + dr * R_iterat
    t = Hbar**2 / (2 * ReduceMass * dr**2)
    
    Allocate(CPOT(Nlevelmax,Nlevelmax,0:R_iterat+2))
    Allocate(CPOT0(Nlevelmax,Nlevelmax))
    Allocate(CPOTH(Nlevelmax,Nlevelmax))

    ! The variable for discrete basis method
    Nx = R_iterat
    Nx_cc = Nlevel * Nx
    Nx1 = Nx_cc + Nlevel
    Nx2 = Nx_cc + 2 * Nlevel
    Allocate(Hamilton(Nx1, Nx2))
    Allocate(Ec_psi(Nx2, Nx1))
    Allocate(Ec_psi0(Nx2))
    Allocate(Ec_psi0_2(Nx1))
    Allocate(Echannel(1:Nlevelmax))
    Allocate(a_tra(Nlevel))
    Allocate(wf(Nx2))

    ! The variable for R matrix
    Call gaussh(Nbase, xgauss, wgauss)

    Call Output_information()

    Call coupled_matrix0

    Do ir = 0, R_iterat + 1
      r = R_min + dr * ir
      Call coupled_matrix(r, CPOT0)
      Do  i = 1, Nlevel
        Do  j = 1, Nlevel
          CPOT(i, j, ir) = CPOT0(i, j)
        End Do
      End Do
    End Do

    rh = R_min + dr / 2.0
    Call coupled_matrix(rh, CPOTH)


  End Subroutine Initialize_CCFull


End Module ccfull_initialization_mod