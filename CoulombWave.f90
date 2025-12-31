Subroutine myCOULFG(angL,eta,rho,myf,mydf,myg,mydg,iv)
  double precision  angL,eta, rho, myf, mydf,myg,mydg
  integer iv,lmax 
  double precision XLMIN, XLMAX
  double precision XX, ETA1
  integer  MODE1, KFN, IFAIL
  double precision , allocatable, dimension(:):: FC,GC,FCP,GCP
  intent(in) angL
  intent(in) eta
  intent(in) rho
  intent(out) myf
  intent(out) mydf
  intent(out) myg
  intent(out) mydg
  intent(out) iv     
  XLMIN=angL
  XLMAX=angL
  lmax=IDINT(XLMAX)+1
  allocate (FC(lmax),GC(lmax),FCP(lmax),GCP(lmax))
  MODE1=1
  KFN=0 
  XX=rho
  ETA1=eta
  call COULFG(XX,ETA1,XLMIN,XLMAX,FC,GC,FCP,GCP,MODE1,KFN,IFAIL)
  myf=FC(lmax)
  mydf=FCP(lmax)
  myg=GC(lmax)
  mydg=GCP(lmax) 
  iv=IFAIL
  return 
End Subroutine

SUBROUTINE COULFG(XX, ETA1, XLMIN, XLMAX, FC, GC, FCP, GCP, MODE1, KFN, IFAIL)
    IMPLICIT NONE
    ! Argument type declarations
    DOUBLE PRECISION, INTENT(IN) :: XX, ETA1, XLMIN, XLMAX
    DOUBLE PRECISION, DIMENSION(*), INTENT(OUT) :: FC, GC, FCP, GCP
    INTEGER, INTENT(IN) :: MODE1, KFN
    INTEGER, INTENT(OUT) :: IFAIL

    ! Local variable declarations
    DOUBLE PRECISION :: ACCUR, ACC, ACC4, ACCH, ETA, GJWKB, PACCQ, X, XLM, E2MM1, DELL, XLL, XI, FCL, PK, PX, F, D, C
    DOUBLE PRECISION :: PK1, EK, RK2, TK, DF, FPL, XL, RL, EL, SL, FCL1, TA, WI, P, Q, AR, AI, BR, BI, DR, DI, DP, DQ
    DOUBLE PRECISION :: A, B, GAM, W, ALPHA, BETA, FCM, GCL, GPL, FJWKB, GCL1
    INTEGER :: MODE, NFP, NPQ, IEXP, M1, LXTRA, L1, LP, MAXL, L, LP1
    LOGICAL :: ETANE0, XLTURN

    ! Common block
    COMMON /STEE/ PACCQ, NFP, NPQ, IEXP, M1

    ! Constants
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, TEN2 = 1.0D2, ABORT = 2.0D4
    DOUBLE PRECISION, PARAMETER :: HALF = 0.5D0, TM30 = 1.0D-30, BIG = 1.0D+100
    DOUBLE PRECISION, PARAMETER :: RT2DPI = 0.7978845608286535587989211986876373D0

    ! Initialize accuracy and other variables
    ACCUR = 1.0D-16
    MODE = 1
    IF (MODE1 == 2 .OR. MODE1 == 3) MODE = MODE1
    IFAIL = 0
    IEXP = 1
    NPQ = 0
    ETA = ETA1
    GJWKB = ZERO
    PACCQ = ONE
    IF (KFN /= 0) ETA = ZERO
    ETANE0 = (ETA /= ZERO)
    ACC = ACCUR * 10.0D0
    ACC4 = ACC * TEN2 * TEN2
    ACCH = SQRT(ACC)

    ! Check for small XX
    IF (XX <= ACCH) THEN
        IFAIL = -1
        RETURN
    END IF

    X = XX
    XLM = XLMIN
    IF (KFN == 2) XLM = XLM - HALF
    IF (XLM <= -ONE .OR. XLMAX < XLMIN) THEN
        IFAIL = -2
        RETURN
    END IF

    E2MM1 = ETA*ETA + XLM*XLM + XLM
    XLTURN = (X*(X - TWO*ETA) < XLM*XLM + XLM)
    DELL = XLMAX - XLMIN + ACC
    LXTRA = INT(DELL)
    XLL = XLM + DBLE(LXTRA)

    M1 = MAX(INT(XLMIN + ACC), 0) + 1
    L1 = M1 + LXTRA

    XI = ONE / X
    FCL = ONE
    PK = XLL + ONE
    PX = PK + ABORT
    F = ETA/PK + PK*XI
    IF (ABS(F) < TM30) F = TM30
    D = ZERO
    C = F

    DO
        PK1 = PK + ONE
        EK = ETA / PK
        RK2 = ONE + EK*EK
        TK = (PK + PK1)*(XI + EK/PK1)
        D = TK - RK2 * D
        C = TK - RK2 / C
        IF (ABS(C) < TM30) C = TM30
        IF (ABS(D) < TM30) D = TM30
        D = ONE / D
        DF = D * C
        F = F * DF
        IF (D < ZERO) FCL = -FCL
        PK = PK1
        IF (PK > PX) THEN
            IFAIL = -3
            RETURN
        END IF
        IF (ABS(DF-ONE) < ACC) EXIT
    END DO

    NFP = PK - XLL - 1
    IF (LXTRA /= 0) THEN
        FCL = FCL / BIG
        FPL = FCL * F
        IF (MODE == 1) FCP(L1) = FPL
        FC(L1) = FCL
        XL = XLL
        RL = ONE
        EL = ZERO
        DO LP = 1, LXTRA
            IF (ETANE0) EL = ETA / XL
            IF (ETANE0) RL = SQRT(ONE + EL*EL)
            SL = EL + XL*XI
            L = L1 - LP
            FCL1 = (FCL * SL + FPL) / RL
            FPL = FCL1 * SL - FCL * RL
            FCL = FCL1
            FC(L) = FCL
            IF (MODE == 1) FCP(L) = FPL
            IF (MODE /= 3 .AND. ETANE0) GC(L+1) = RL
            IF (ABS(FCL) > BIG) THEN
                DO LP1 = L, M1+LXTRA
                    IF (MODE == 1) FCP(LP1) = FCP(LP1) * 1.0D-20
                    FC(LP1) = FC(LP1) * 1.0D-20
                END DO
                FCL = FC(L)
                FPL = FPL * 1.0D-20
            END IF
            XL = XL - ONE
        END DO
        IF (FCL == ZERO) FCL = ACC
        F = FPL / FCL
    END IF

    IF (XLTURN) CALL JWKB(X, ETA, MAX(XLM, ZERO), FJWKB, GJWKB, IEXP)
    IF (IEXP > 1 .OR. GJWKB > ONE/(ACCH*TEN2)) THEN
        W = FJWKB
        GAM = GJWKB * W
        P = F
        Q = ONE
    ELSE
        XLTURN = .FALSE.
        TA = TWO * ABORT
        PK = ZERO
        WI = ETA + ETA
        P = ZERO
        Q = ONE - ETA*XI
        AR = -E2MM1
        AI = ETA
        BR = TWO*(X - ETA)
        BI = TWO
        DR = BR/(BR*BR + BI*BI)
        DI = -BI/(BR*BR + BI*BI)
        DP = -XI*(AR*DI + AI*DR)
        DQ = XI*(AR*DR - AI*DI)
        DO
            P = P + DP
            Q = Q + DQ
            PK = PK + TWO
            AR = AR + PK
            AI = AI + WI
            BI = BI + TWO
            D = AR*DR - AI*DI + BR
            DI = AI*DR + AR*DI + BI
            C = ONE/(D*D + DI*DI)
            DR = C * D
            DI = -C * DI
            A = BR*DR - BI*DI - ONE
            B = BI*DR + BR*DI
            C = DP*A - DQ*B
            DQ = DP*B + DQ*A
            DP = C
            IF (PK > TA) THEN
                IFAIL = -4
                RETURN
            END IF
            IF (ABS(DP) + ABS(DQ) < (ABS(P) + ABS(Q)) * ACC) EXIT
        END DO
        NPQ = INT(PK / TWO)
        PACCQ = HALF * ACC / MIN(ABS(Q), ONE)
        IF (ABS(P) > ABS(Q)) PACCQ = PACCQ * ABS(P)
        GAM = (F - P) / Q
        IF (Q <= ACC4 * ABS(P)) THEN
            IFAIL = -5
            RETURN
        END IF
        W = ONE / SQRT((F - P)*GAM + Q)
    END IF

    ALPHA = ZERO
    IF (KFN == 1) ALPHA = XI
    IF (KFN == 2) ALPHA = XI * HALF
    BETA = ONE
    IF (KFN == 1) BETA = XI
    IF (KFN == 2) BETA = SQRT(XI) * RT2DPI
    FCM = SIGN(W, FCL) * BETA
    FC(M1) = FCM
    IF (MODE /= 3) THEN
        IF (.NOT. XLTURN) THEN
            GCL = FCM * GAM
        ELSE
            GCL = GJWKB * BETA
        END IF
        IF (KFN /= 0) GCL = -GCL
        GC(M1) = GCL
        GPL = GCL*(P - Q/GAM) - ALPHA*GCL
        IF (MODE /= 2) THEN
            GCP(M1) = GPL
            FCP(M1) = FCM*(F - ALPHA)
        END IF
    END IF

    IF (LXTRA == 0) RETURN

    W = BETA * W / ABS(FCL)
    MAXL = L1 - 1
    DO L = M1, MAXL
        IF (MODE /= 3) THEN
            XL = XL + ONE
            IF (ETANE0) EL = ETA / XL
            IF (ETANE0) RL = GC(L+1)
            SL = EL + XL*XI
            GCL1 = ((SL - ALPHA)*GCL - GPL) / RL
            IF (ABS(GCL1) > BIG) THEN
                IFAIL = L1 - L
                RETURN
            END IF
            GPL = RL*GCL - (SL + ALPHA)*GCL1
            GCL = GCL1
            GC(L+1) = GCL1
            IF (MODE /= 2) THEN
                GCP(L+1) = GPL
                FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
            END IF
        END IF
        FC(L+1) = W * FC(L+1)
    END DO

END SUBROUTINE COULFG

subroutine JWKB(XX, ETA1, XL, FJWKB, GJWKB, IEXP)
    implicit none
    ! ===== 输入输出参数 =====
    real(8), intent(in)  :: XX, ETA1, XL
    real(8), intent(out) :: FJWKB, GJWKB
    integer, intent(out) :: IEXP

    ! ===== 内部变量 =====
    real(8) :: ZERO, HALF, ONE, SIX, TEN
    real(8) :: DZ, RL35
    real(8) :: aloge, X, ETA, GH2, XLL1, HLL, HL, SL, RL2, GH
    real(8) :: PHI, PHI10

    ! ===== 常量赋值 =====
    ZERO = 0.0d0
    HALF = 0.5d0
    ONE  = 1.0d0
    SIX  = 6.0d0
    TEN  = 10.0d0
    DZ   = 0.0d0
    RL35 = 35.0d0

    aloge = log(TEN)

    ! ===== 核心计算 =====
    X   = XX
    ETA = ETA1

    GH2  = X * (2.0d0*ETA - X)
    XLL1 = max(XL*XL + XL, DZ)

    if (GH2 + XLL1 <= ZERO) then
        FJWKB = 0.0d0
        GJWKB = 0.0d0
        IEXP  = 0
        return
    end if

    HLL = XLL1 + SIX / RL35
    HL  = sqrt(HLL)
    SL  = ETA/HL + HL/X
    RL2 = ONE + ETA*ETA/HLL
    GH  = sqrt(GH2 + HLL) / X

    PHI = X*GH - HALF * ( HL*log((GH+SL)**2 / RL2) - log(GH) )
    if (ETA /= ZERO) PHI = PHI - ETA*atan2(X*GH, X - ETA)

    PHI10 = -PHI * aloge
    IEXP  = int(PHI10)

    if (IEXP > 70) then
        GJWKB = TEN**(PHI10 - dble(IEXP))
    else
        GJWKB = exp(-PHI)
        IEXP  = 0
    end if

    FJWKB = HALF / (GH * GJWKB)

end subroutine JWKB

Subroutine WHIT(HETA, R, XK, E, LL, F, FD, IE)
    implicit none
    ! ====== 输入输出参数 ======
    integer, intent(in) :: LL
    integer, intent(inout) :: IE
    real(8), intent(in) :: HETA, R, XK, E
    real(8), intent(out) :: F(LL+1), FD(LL+1)

    ! ====== 内部变量 ======
    integer :: L, LP1, LM, LMP1, IS, H
    integer :: I, M, M1, M2, JS, IFEQL
    real(8) :: fpmax, fpminl
    real(8) :: EE, AK, ETA, RHO, RHOA, PJE, A, B, C, D
    real(8) :: T(12), S(7)

    fpmax = 1.0d290
    L = LL + 1

    EE = -1.0d0
    AK = XK
    ETA = HETA
    LP1 = L + 1
    RHO = AK * R

    S(:) = 0.0d0

    ! ====== 设置 LM ======
    if (L <= 50) then
        LM = 60
    else
        LM = L + 10
    end if
    LMP1 = LM + 1
    IS = 7

    PJE = 30.0d0 * RHO + 1.0d0
    H = max(int(PJE), 4)
    H = RHO / H

    RHOA = 10.0d0 * (ETA + 1.0d0)
    if (RHOA <= RHO) then
        IFEQL = 1
        RHOA = RHO
    else
        IFEQL = 0
    end if

    PJE = RHOA / H + 0.5d0
    RHOA = H * int(PJE)

    if (IFEQL == 1) then
        if (RHOA <= RHO + 1.5d0 * H) then
            RHOA = RHO + 2.0d0 * H
        end if
    end if

    if (EE < 0.0d0) then
        ! ----------- 级数展开部分 ----------
        C = 1.0d0 / RHOA
        A = 1.0d0
        B = 1.0d0 - C * ETA
        F(1) = A
        FD(1) = B
        do M = 1, 26
            D = 0.5d0 * (ETA + dble(M - 1)) * (ETA + dble(M)) * C / dble(M)
            A = -A * D
            B = -B * D - A * C
            F(1) = F(1) + A
            FD(1) = FD(1) + B
        end do

        A = -ETA * log(2.0d0 * RHOA) - RHOA
        fpminl = -log(fpmax)
        if (IE == 0 .and. A < fpminl) IE = int(fpminl - A)
        A = exp(A + dble(IE))

        F(1) = A * F(1)
        FD(1) = A * FD(1) * (-1.0d0 - 2.0d0 * ETA / RHOA)

        if (IFEQL == 1) then
            S(IS) = F(1)
            if (IS == 7) then
                IS = 6
                RHOA = RHOA + H
                ! 再做一次展开
                C = 1.0d0 / RHOA
                A = 1.0d0
                B = 1.0d0 - C * ETA
                F(1) = A
                FD(1) = B
                do M = 1, 26
                    D = 0.5d0 * (ETA + dble(M - 1)) * (ETA + dble(M)) * C / dble(M)
                    A = -A * D
                    B = -B * D - A * C
                    F(1) = F(1) + A
                    FD(1) = FD(1) + B
                end do
                A = -ETA * log(2.0d0 * RHOA) - RHOA
                fpminl = -log(fpmax)
                if (IE == 0 .and. A < fpminl) IE = int(fpminl - A)
                A = exp(A + dble(IE))
                F(1) = A * F(1)
                FD(1) = A * FD(1) * (-1.0d0 - 2.0d0 * ETA / RHOA)
            end if
        else
            F(1) = T(1)
            FD(1) = T(2)
        end if

    else
        stop "WHIT1: EE >= 0 not implemented"
    end if

    ! ----------- 递推关系 ----------
    C = 1.0d0 / RHO
    do M = 1, L-1
        A = ETA / dble(M)
        B = A + C * dble(M)
        F(M+1) = (B * F(M) - FD(M)) / (A + 1.0d0)
        FD(M+1) = (A - 1.0d0) * F(M) - B * F(M+1)
    end do

    do M = 1, L
        FD(M) = AK * FD(M)
    end do

    return
End Subroutine