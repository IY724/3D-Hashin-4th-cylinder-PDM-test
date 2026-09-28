C=======================================================================
C B（能量线性软化） UMAT: local +/-theta plies, 3D Hashin damage.
C Material slots are read and checked by WCM_CHECK below.
C Elastic: E11/E22/E33=1:3, PR12/PR13/PR23=4:6,
C          G12/G13/G23=7:9 (MPa; Poisson ratios unitless).
C Strength: XT/XC/YT/YC=10:13, S12/S13/S23=14:16 (MPa).
C B: GFT/GFC/GMT/GMC=17:20 (fracture energy, N/mm).
C    ALPHA=21 (reserved), SMT/SMC=22:23 (shear factors).
C    ETA_F/ETA_M=24:25 (viscosity in step-time units).
C    THETA/WPLUS=26:27 (deg, fraction).
C    CAPF/CAPM=28:29 (stiffness loss limits, 0..1).
C    UPDATE=30 (0: next increment, 1: current feedback).
C=======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER
      INTEGER KSPT,JSTEP(4),KINC,J,IERR,IADV,MODE,IREASON,NREQ
      PARAMETER (MODE=2)
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
      DOUBLE PRECISION DDSDDE(NTENS,NTENS),DDSDDT(NTENS)
      DOUBLE PRECISION DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS)
      DOUBLE PRECISION TIME(2),PREDEF(*),DPRED(*),PROPS(NPROPS)
      DOUBLE PRECISION COORDS(3),DROT(3,3),DFGRD0(3,3)
      DOUBLE PRECISION DFGRD1(3,3),SSE,SPD,SCD,RPL,DRPLDT
      DOUBLE PRECISION DTIME,TEMP,DTEMP,PNEWDT,CELENT
      DOUBLE PRECISION NEW(84),OLD(84),S(6),C(6,6)
      DOUBLE PRECISION SOLD(6),UOLD,WORK
C--- Abaqus outputs and interface check.
      RPL=0D0
      DRPLDT=0D0
      DDSDDT=0D0
      DRPLDE=0D0
      SCD=0D0
      NREQ=20
      IF(MODE.EQ.2) NREQ=62
      IF(NDI.NE.3.OR.NSHR.NE.3.OR.NTENS.NE.6.OR.
     1   NSTATV.NE.NREQ) THEN
          WRITE(7,*) 'WCM_ERROR dimensions: ',NDI,NSHR,NTENS,
     1                NSTATV
          CALL XIT
          RETURN
      END IF
      SOLD=STRESS
      UOLD=SSE
      IADV=1
      IF(KINC.EQ.0.OR.DTIME.EQ.0D0) IADV=0
C--- Restore the two independent angle-family histories.
      CALL WCM_UNPACK(MODE,STATEV,NSTATV,OLD)
C--- Stress, failure criteria, damage, and DDSDDE.
      CALL WCM_UPDATE(MODE,PROPS,NPROPS,OLD,STRAN,DSTRAN,
     1                DTIME,CELENT,IADV,NEW,S,C,IERR)
      IF(IERR.NE.0) THEN
C 未收敛的Newton试应变也可能越界；先请求切步，不提交坏历史。
C 真正的材料卡错误直接退出；持续不满足能量条件则由最小步长
C 终止，不能通过无限切步或擅自增大G来修复其物理不相容性。
          WRITE(7,*) 'WCM_ERROR code/mode/element/point: ',
     1               IERR,MODE,NOEL,NPT
          WRITE(7,*) 'WCM_ERROR material/length: ',CMNAME,CELENT
          IREASON=MOD(IERR,1000)
          IF(IERR.GE.1000.AND.IREASON.GE.21.AND.
     1       IREASON.LE.34.AND.DTIME.GT.0D0) THEN
              PNEWDT=MIN(PNEWDT,0.5D0)
              CALL WCM_UPDATE(MODE,PROPS,NPROPS,OLD,STRAN,
     1             DSTRAN,DTIME,CELENT,0,NEW,S,C,IERR)
              IF(IERR.EQ.0) THEN
                  DDSDDE=C
                  STRESS=S
                  RETURN
              END IF
          END IF
          CALL XIT
          RETURN
      END IF
C--- Return stress and the material Jacobian.
      STRESS=S
      DDSDDE=C
C NLGEOM中的体积项来自J*sigma的线性化；不对sigma重复旋转。
C 此处仍采用原项目随动小材料应变框架，非有限应变超弹性。
      IF(JSTEP(3).EQ.1) THEN
          DO J=1,3
              DDSDDE(:,J)=DDSDDE(:,J)+STRESS
          END DO
      END IF
C--- Save damage histories and numbered SDV results.
      IF(IADV.EQ.1) THEN
          CALL WCM_PACK(MODE,PROPS,NPROPS,OLD,NEW,STATEV,NSTATV)
          SSE=0.5D0*DOT_PRODUCT(STRESS,STRAN+DSTRAN)
C SPD为离散应力功减弹性能变化的有符号诊断，含算法/黏性影响。
C 保留旧本构与分步时序时不能将它当作精确断裂耗散。
C 不强行截掉负值来掩盖能量误差；基准测试检查步长敏感性。
          WORK=0.5D0*DOT_PRODUCT(SOLD+STRESS,DSTRAN)
          SPD=SPD+WORK-(SSE-UOLD)
      END IF
      RETURN
      END

C Persistent history is separate for the two weighted families.
      SUBROUTINE WCM_UNPACK(MODE,STATE,NS,OLD)
      IMPLICIT NONE
      INTEGER MODE,NS
      DOUBLE PRECISION STATE(NS),OLD(84)
      OLD=0D0
      IF(MODE.EQ.1) THEN
          OLD(7:10)=STATE(9:12)
          OLD(47:50)=STATE(13:16)
      ELSE
          OLD(1:24)=STATE(11:34)
          OLD(41:64)=STATE(35:58)
      END IF
      OLD(81:84)=STATE(NS-3:NS)
      END

C Combine modes within a family before taking family maxima.
      SUBROUTINE WCM_PACK(MODE,P,NP,OLD,NEW,STATE,NS)
      IMPLICIT NONE
      INTEGER MODE,NP,NS,K,IB,IT,IUP,IUSED
      DOUBLE PRECISION P(NP),OLD(84),NEW(84),STATE(NS),W
      DOUBLE PRECISION D(4),V(4),U(4),DF,DM,UF,UM
      STATE=0D0
      STATE(1:4)=-1D100
      IT=23
      IUP=25
      IUSED=7
      IF(MODE.EQ.1) THEN
          STATE(9:12)=NEW(7:10)
          STATE(13:16)=NEW(47:50)
      ELSE
          IT=26
          IUP=30
          IUSED=9
          STATE(11:34)=NEW(1:24)
          STATE(35:58)=NEW(41:64)
      END IF
      STATE(NS-3:NS)=NEW(81:84)
      DO K=1,2
          IB=40*(K-1)
          W=P(IT+1)
          IF(K.EQ.2) W=1D0-W
          IF(W.EQ.0D0) CYCLE
          STATE(1:4)=MAX(STATE(1:4),NEW(IB+25:IB+28))
          IF(MODE.EQ.1) THEN
              D=NEW(IB+7:IB+10)
              U=OLD(IB+7:IB+10)*P(17:20)
              IF(NINT(P(IUP)).EQ.1) U=D*P(17:20)
          ELSE
              D=NEW(IB+1:IB+4)
              V=NEW(IB+21:IB+24)
              STATE(7)=MAX(STATE(7),
     1                    1D0-(1D0-V(1))*(1D0-V(2)))
              STATE(8)=MAX(STATE(8),
     1                    1D0-(1D0-V(3))*(1D0-V(4)))
              U=OLD(IB+21:IB+24)
              IF(NINT(P(IUP)).EQ.1) U=V
          END IF
          DF=1D0-(1D0-D(1))*(1D0-D(2))
          DM=1D0-(1D0-D(3))*(1D0-D(4))
          UF=1D0-(1D0-U(1))*(1D0-U(2))
          UM=1D0-(1D0-U(3))*(1D0-U(4))
          IF(MODE.EQ.2) THEN
              UF=MIN(UF,P(28))
              UM=MIN(UM,P(29))
          END IF
          STATE(5)=MAX(STATE(5),DF)
          STATE(6)=MAX(STATE(6),DM)
          STATE(IUSED)=MAX(STATE(IUSED),UF)
          STATE(IUSED+1)=MAX(STATE(IUSED+1),UM)
      END DO
      END

C=======================================================================
C WCM UMAT core: parameter check, local rotation, stress, Hashin, damage.
C Abaqus orientation provides W axes; +/-theta rotate W to fiber L axes.
C A: fixed degradation. B: energy softening and viscous damage.
C Units follow the input card: MPa, mm, N/mm and step time.
C=======================================================================
      SUBROUTINE WCM_UPDATE(MODE,P,NP,OLD,E0,DE,DT,LC,IADV,
     1                     NEW,S,C,IERR)
      IMPLICIT NONE
      INTEGER MODE,NP,IADV,IERR,I
      DOUBLE PRECISION P(NP),OLD(84),NEW(84),E0(6),DE(6)
      DOUBLE PRECISION DT,LC,S(6),C(6,6),C0(6,6),E1(6)
      IERR=0
      NEW=OLD
      S=0D0
      C=0D0

C--- Material card and undamaged stiffness.
      CALL WCM_CHECK(MODE,P,NP,OLD,DT,LC,IERR)
      IF(IERR.NE.0) RETURN
      CALL WCM_ELASTIC(P,C0,IERR)
      IF(IERR.NE.0) RETURN

C--- Advance total strain.
      E1=E0+DE
      DO I=1,6
          IF(E1(I).NE.E1(I).OR.ABS(E1(I)).GT.1D10) THEN
              IERR=11
              RETURN
          END IF
      END DO

C--- Angle conversion, stress trial, criterion and damage.
      CALL WCM_EVAL(MODE,P,NP,C0,OLD,E1,DT,LC,IADV,
     1              NEW,S,C,IERR)
      IF(IERR.NE.0) RETURN

C--- B current feedback: DDSDDE from the same OLD history.
      IF(MODE.EQ.2.AND.IADV.EQ.1) THEN
          IF(NINT(P(30)).EQ.1)
     1      CALL WCM_TANGENT_B(P,NP,C0,OLD,E1,DT,LC,C,IERR)
      END IF
      END

      SUBROUTINE WCM_TANGENT_B(P,NP,C0,OLD,E,DT,LC,C,IERR)
      IMPLICIT NONE
      INTEGER NP,IERR,I,J
      DOUBLE PRECISION P(NP),C0(6,6),OLD(84),E(6),DT,LC,C(6,6)
      DOUBLE PRECISION SP(6),SM(6),EP(6),EM(6),H
      DOUBLE PRECISION TMP(84),CT(6,6)
      DO J=1,6
          H=MAX(1D-9,1D-6*ABS(E(J)))
          EP=E
          EM=E
          EP(J)=EP(J)+H
          EM(J)=EM(J)-H
          CALL WCM_EVAL(2,P,NP,C0,OLD,EP,DT,LC,1,
     1                  TMP,SP,CT,IERR)
          IF(IERR.NE.0) RETURN
          CALL WCM_EVAL(2,P,NP,C0,OLD,EM,DT,LC,1,
     1                  TMP,SM,CT,IERR)
          IF(IERR.NE.0) RETURN
          DO I=1,6
              C(I,J)=(SP(I)-SM(I))/(2D0*H)
          END DO
      END DO
      END

      SUBROUTINE WCM_CHECK(MODE,PROPS,NPROPS,OLD,DT,LC,IERR)
      IMPLICIT NONE
      INTEGER MODE,NPROPS,IERR,I,IB
      DOUBLE PRECISION PROPS(NPROPS),OLD(84),DT,LC
      DOUBLE PRECISION E11,E22,E33,PR12,PR13,PR23
      DOUBLE PRECISION G12,G13,G23,XT,XC,YT,YC,S12,S13,S23
      DOUBLE PRECISION DFT0,DFC0,DMT0,DMC0,SMT,SMC
      DOUBLE PRECISION GFT,GFC,GMT,GMC,ALPHA,ETA_F,ETA_M
      DOUBLE PRECISION THETA,WPLUS,CAPF,CAPM,UPDATE
      IERR=0
      IF ((MODE.EQ.1.AND.NPROPS.NE.25).OR.
     1    (MODE.EQ.2.AND.NPROPS.NE.30).OR.
     2    (MODE.NE.1.AND.MODE.NE.2)) THEN
          IERR=1
          RETURN
      END IF
      DO I=1,NPROPS
          IF(PROPS(I).NE.PROPS(I).OR.
     1       ABS(PROPS(I)).GT.1D100) IERR=2
      END DO
      DO I=1,84
          IF(OLD(I).NE.OLD(I).OR.ABS(OLD(I)).GT.1D100) IERR=3
      END DO
      IF (IERR.NE.0) RETURN
C--- Elastic constants: MPa for modulus, dimensionless Poisson ratio.
      E11=PROPS(1)
      E22=PROPS(2)
      E33=PROPS(3)
      PR12=PROPS(4)
      PR13=PROPS(5)
      PR23=PROPS(6)
      G12=PROPS(7)
      G13=PROPS(8)
      G23=PROPS(9)
C--- Ply strengths in MPa; S13 stays in the legacy slot.
      XT=PROPS(10)
      XC=PROPS(11)
      YT=PROPS(12)
      YC=PROPS(13)
      S12=PROPS(14)
      S13=PROPS(15)
      S23=PROPS(16)
      IF(MIN(E11,E22,E33,G12,G13,G23,
     1       XT,XC,YT,YC,S12,S13,S23).LE.0D0) IERR=6
C--- Variant-specific damage inputs; original slot order is retained.
      IF(MODE.EQ.1) THEN
          DFT0=PROPS(17)
          DFC0=PROPS(18)
          DMT0=PROPS(19)
          DMC0=PROPS(20)
          SMT=PROPS(21)
          SMC=PROPS(22)
          THETA=PROPS(23)
          WPLUS=PROPS(24)
          UPDATE=PROPS(25)
          IF(MIN(DFT0,DFC0,DMT0,DMC0,SMT,SMC).LT.0D0.OR.
     1       MAX(DFT0,DFC0,DMT0,DMC0,SMT,SMC).GT.1D0)
     2       IERR=9
      ELSE
          GFT=PROPS(17)
          GFC=PROPS(18)
          GMT=PROPS(19)
          GMC=PROPS(20)
          ALPHA=PROPS(21)
          SMT=PROPS(22)
          SMC=PROPS(23)
          ETA_F=PROPS(24)
          ETA_M=PROPS(25)
          THETA=PROPS(26)
          WPLUS=PROPS(27)
          CAPF=PROPS(28)
          CAPM=PROPS(29)
          UPDATE=PROPS(30)
          IF(MIN(GFT,GFC,GMT,GMC).LE.0D0.OR.
     1       MIN(SMT,SMC,ETA_F,ETA_M).LT.0D0.OR.
     2       MAX(SMT,SMC).GT.1D0.OR.
     3       MIN(CAPF,CAPM).LE.0D0.OR.
     4       MAX(CAPF,CAPM).GT.1D0) IERR=9
C         ALPHA is a retained legacy input; the current FI is fixed.
          IF(ALPHA.NE.ALPHA) IERR=2
          IF(LC.NE.LC.OR.LC.LE.0D0.OR.LC.GT.1D100) IERR=10
      END IF
      IF(THETA.LT.0D0.OR.THETA.GT.90D0.OR.
     1   WPLUS.LT.0D0.OR.WPLUS.GT.1D0) IERR=7
      IF(UPDATE.NE.0D0.AND.UPDATE.NE.1D0) IERR=8
      IF(DT.NE.DT.OR.DT.LT.0D0.OR.DT.GT.1D100) IERR=10
      IF (OLD(84).NE.0D0) THEN
          IF (OLD(84).NE.202609D0.OR.OLD(83).NE.DBLE(MODE))
     1        IERR=5
          IF(ABS(OLD(81)-THETA).GT.1D-12.OR.
     1       ABS(OLD(82)-WPLUS).GT.1D-12) IERR=5
      ELSE IF (ANY(OLD(1:80).NE.0D0)) THEN
          IERR=5
      END IF
      IF (IERR.NE.0) RETURN
C--- Validate the independent +/-theta damage histories.
      DO IB=0,40,40
          IF(MODE.EQ.1) THEN
              DO I=7,10
                  IF(OLD(IB+I).NE.0D0.AND.OLD(IB+I).NE.1D0)
     1                IERR=13
              END DO
          ELSE
              IF(MINVAL(OLD(IB+1:IB+24)).LT.0D0) IERR=13
              IF(MAXVAL(OLD(IB+1:IB+2)).GT.1D0.OR.
     1           MAXVAL(OLD(IB+3:IB+4)).GT.1D0.OR.
     2           MAXVAL(OLD(IB+21:IB+22)).GT.1D0.OR.
     3           MAXVAL(OLD(IB+23:IB+24)).GT.1D0) IERR=13
              DO I=9,12
                  IF(OLD(IB+I).NE.0D0.AND.OLD(IB+I).NE.1D0)
     1                IERR=13
              END DO
          END IF
      END DO
      END

      SUBROUTINE WCM_ROTATE(THETA,B)
      IMPLICIT NONE
      DOUBLE PRECISION THETA,B(6,6),C,S,T,PI
      PARAMETER (PI=3.1415926535897932384626433832795D0)
      T=THETA*PI/180D0
      C=COS(T)
      S=SIN(T)
      IF (ABS(C).LT.1D-14) C=0D0
      IF (ABS(S).LT.1D-14) S=0D0
      B=0D0
      B(1,1)=C*C
      B(1,2)=S*S
      B(1,4)=C*S
      B(2,1)=S*S
      B(2,2)=C*C
      B(2,4)=-C*S
      B(3,3)=1D0
      B(4,1)=-2D0*C*S
      B(4,2)=2D0*C*S
      B(4,4)=C*C-S*S
      B(5,5)=C
      B(5,6)=S
      B(6,5)=-S
      B(6,6)=C
      END

      SUBROUTINE WCM_ELASTIC(P,C,IERR)
      IMPLICIT NONE
      DOUBLE PRECISION P(*),C(6,6),V21,V31,V32,DEN,FAC
      INTEGER IERR
      C=0D0
      IERR=0
      V21=P(2)*P(4)/P(1)
      V31=P(3)*P(5)/P(1)
      V32=P(3)*P(6)/P(2)
      DEN=1D0-P(4)*V21-P(5)*V31-P(6)*V32
     1    -2D0*V21*V32*P(5)
      IF (1D0-P(4)*V21.LE.1D-12.OR.DEN.LE.1D-12) THEN
          IERR=12
          RETURN
      END IF
      FAC=1D0/DEN
      C(1,1)=P(1)*FAC*(1D0-P(6)*V32)
      C(2,2)=P(2)*FAC*(1D0-P(5)*V31)
      C(3,3)=P(3)*FAC*(1D0-P(4)*V21)
      C(1,2)=P(1)*FAC*(V21+V31*P(6))
      C(1,3)=P(1)*FAC*(V31+V21*V32)
      C(2,3)=P(2)*FAC*(V32+P(4)*V31)
      C(2,1)=C(1,2)
      C(3,1)=C(1,3)
      C(3,2)=C(2,3)
      C(4,4)=P(7)
      C(5,5)=P(8)
      C(6,6)=P(9)
      END

      SUBROUTINE WCM_STIFF(MODE,P,D,C0,C)
      IMPLICIT NONE
      INTEGER MODE,I,J
      DOUBLE PRECISION P(*),D(4),C0(6,6),C(6,6)
      DOUBLE PRECISION DF,DM,RF,RM,RS,RMT,RMC,RN,RMIN
      PARAMETER (RMIN=1D-6)
      DF=1D0-(1D0-D(1))*(1D0-D(2))
      DM=1D0-(1D0-D(3))*(1D0-D(4))
      IF (MODE.EQ.1) THEN
C 严格保留hashin_constant的各系数及各自1e-6数值下限。
          RF=MAX(1D0-DF,1D-6)
          RM=MAX(1D0-DM,1D-6)
          RMT=MAX(1D0-P(21)*D(3),1D-6)
          RMC=MAX(1D0-P(22)*D(4),1D-6)
          RS=RF*RMT*RMC
      ELSE
C P(28:29) cap stiffness degradation, NOT damage history.
          DF=MIN(MAX(DF,0D0),P(28))
          DM=MIN(MAX(DM,0D0),P(29))
          RF=MAX(1D0-DF,RMIN)
          RM=MAX(1D0-DM,RMIN)
          RS=(1D0-P(22)*D(3))*(1D0-P(23)*D(4))
          RS=MIN(MAX(RS,RMIN),1D0)
      END IF
C Floor the combined normal factor, not only each factor.
C This retains a positive elastic stiffness at full damage.
      RN=MAX(RF*RM,RMIN)
      RS=MAX(RS,RMIN)
      C=C0
      DO I=1,3
          DO J=1,3
              C(I,J)=C0(I,J)*RN
          END DO
      END DO
      C(1,1)=C0(1,1)*RF
      C(4,4)=C0(4,4)*RS
      C(5,5)=C0(5,5)*RS
      C(6,6)=C0(6,6)*RS
      END

      SUBROUTINE WCM_HASHIN(S,P,FI)
      IMPLICIT NONE
      DOUBLE PRECISION S(6),P(*),FI(4),Q
C 判据与朱宁-2023式(1.4)-(1.7)一致：FC为平方形式；保持alpha=1的剪切
C 项；tau13用S12(P(14))而非S13(P(15))，与原版唯一差别即此项。
C 原版(基准baseline)此处tau13用S13，本行为有意分歧，勿改回。
      FI=0D0
      IF(S(1).GE.0D0) THEN
          FI(1)=(S(1)/P(10))**2+(S(4)/P(14))**2
     1         +(S(5)/P(14))**2
      ELSE
          FI(2)=(-S(1)/P(11))**2
      END IF
      Q=S(2)+S(3)
      IF(Q.GE.0D0) THEN
          FI(3)=(Q/P(12))**2+(S(6)**2-S(2)*S(3))/P(16)**2
     1         +(S(4)/P(14))**2+(S(5)/P(14))**2
      ELSE
          FI(4)=((P(13)/(2D0*P(16)))**2-1D0)*Q/P(13)
     1         +Q**2/(2D0*P(16))**2
     2         +(S(6)**2-S(2)*S(3))/P(16)**2
     3         +(S(4)/P(14))**2+(S(5)/P(14))**2
      END IF
      END

      SUBROUTINE WCM_EQ(E,S,LC,DEL,SIG)
      IMPLICIT NONE
      DOUBLE PRECISION E(6),S(6),LC,DEL(4),SIG(4),SH
C 原版工程剪应变等效量；不替换为VUMAT的张量剪应变。
      DEL(1)=LC*SQRT(MAX(E(1),0D0)**2+E(4)**2+E(5)**2)
      DEL(2)=LC*MAX(-E(1),0D0)
      SH=E(4)**2+E(5)**2+E(6)**2
      DEL(3)=LC*SQRT(MAX(E(2),0D0)**2+MAX(E(3),0D0)**2+SH)
      DEL(4)=LC*SQRT(MAX(-E(2),0D0)**2+MAX(-E(3),0D0)**2+SH)
      SIG=0D0
      SH=S(4)*E(4)+S(5)*E(5)+S(6)*E(6)
      IF(DEL(1).GT.0D0) SIG(1)=LC*(MAX(S(1),0D0)
     1  *MAX(E(1),0D0)+S(4)*E(4)+S(5)*E(5))/DEL(1)
      IF(DEL(2).GT.0D0) SIG(2)=LC*MAX(-S(1),0D0)
     1  *MAX(-E(1),0D0)/DEL(2)
      IF(DEL(3).GT.0D0) SIG(3)=LC*(MAX(S(2),0D0)
     1  *MAX(E(2),0D0)+MAX(S(3),0D0)*MAX(E(3),0D0)+SH)/DEL(3)
      IF(DEL(4).GT.0D0) SIG(4)=LC*(MAX(-S(2),0D0)
     1  *MAX(-E(2),0D0)+MAX(-S(3),0D0)*MAX(-E(3),0D0)+SH)/DEL(4)
      END

      SUBROUTINE WCM_EVAL(MODE,P,NP,C0,OLD,E,DT,LC,IADV,
     1                    NEW,S,C,IERR)
      IMPLICIT NONE
      INTEGER MODE,NP,IADV,IERR,K,IT,IB
      DOUBLE PRECISION P(NP),C0(6,6),OLD(84),E(6),DT,LC
      DOUBLE PRECISION NEW(84),S(6),C(6,6),B(6,6),EL(6)
      DOUBLE PRECISION SL(6),CL(6,6),LOCAL(40),TH,W,WP
      NEW=OLD
      S=0D0
      C=0D0
      IERR=0
      IT=23
      IF(MODE.EQ.2) IT=26
      WP=P(IT+1)
      DO K=1,2
          IB=40*(K-1)
          TH=P(IT)
          W=WP
          IF(K.EQ.2) THEN
              TH=-TH
              W=1D0-WP
          END IF
C 零占比族跳过，避免单族验证中无权重族触发错误或损伤。
          IF(W.EQ.0D0) CYCLE
C--- W axes to the local +/-theta fiber axes.
          CALL WCM_ROTATE(TH,B)
          EL=MATMUL(B,E)
C--- Local stress, Hashin criterion and damage update.
          CALL WCM_PLY(MODE,P,C0,OLD(IB+1:IB+40),EL,DT,LC,
     1                 IADV,LOCAL,SL,CL,IERR)
          IF(IERR.NE.0) THEN
C 1000/2000标记出错族，余数为实际错误码。
              IERR=IERR+1000*K
              RETURN
          END IF
          NEW(IB+1:IB+40)=LOCAL
C--- Rotate stress and tangent back; sum the two families.
          S=S+W*MATMUL(TRANSPOSE(B),SL)
          C=C+W*MATMUL(TRANSPOSE(B),MATMUL(CL,B))
      END DO
      IF(IADV.EQ.1) THEN
          NEW(81)=P(IT)
          NEW(82)=WP
          NEW(83)=DBLE(MODE)
          NEW(84)=202609D0
      END IF
      END

      SUBROUTINE WCM_PLY(MODE,P,C0,OLD,E,DT,LC,IADV,
     1                   NEW,S,C,IERR)
      IMPLICIT NONE
      INTEGER MODE,IADV,IERR,IUP
      DOUBLE PRECISION P(*),C0(6,6),OLD(40),E(6),DT,LC
      DOUBLE PRECISION NEW(40),S(6),C(6,6),D(4),FI(4),EFF(6)
      IERR=0
      NEW=OLD
      D=0D0
      IUP=25
      IF(MODE.EQ.2) IUP=30
      IF(MODE.EQ.1) THEN
          D=OLD(7:10)*P(17:20)
      ELSE
          D=OLD(21:24)
      END IF

C--- Trial stress and frozen-damage tangent.
      CALL WCM_STIFF(MODE,P,D,C0,C)
      S=MATMUL(C,E)
      EFF=S
      IF(MODE.EQ.2) EFF=MATMUL(C0,E)

C--- Four Hashin initiation indices.
      CALL WCM_HASHIN(EFF,P,FI)
      IF(IADV.EQ.1) THEN
          IF(MODE.EQ.1) THEN
              CALL WCM_DAMAGE_A(P,FI,NEW)
              D=NEW(7:10)*P(17:20)
          ELSE
              CALL WCM_DAMAGE_B(P,OLD,E,EFF,FI,DT,LC,NEW,IERR)
              IF(IERR.NE.0) RETURN
              D=NEW(21:24)
          END IF
C--- Current-damage feedback if UPDATE=1.
          IF(NINT(P(IUP)).EQ.1) THEN
              CALL WCM_STIFF(MODE,P,D,C0,C)
              S=MATMUL(C,E)
          END IF
          NEW(25:28)=FI
          NEW(29:34)=S
          NEW(35:40)=E
      END IF
      END

      SUBROUTINE WCM_DAMAGE_A(P,FI,NEW)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION P(*),FI(4),NEW(40),D(4)
C--- A: latch failure and apply fixed stiffness factors.
      DO I=1,4
          IF(FI(I).GE.1D0) NEW(I+6)=1D0
      END DO
      D=NEW(7:10)*P(17:20)
      NEW(1:4)=FI
      NEW(5)=1D0-(1D0-D(1))*(1D0-D(2))
      NEW(6)=1D0-(1D0-D(3))*(1D0-D(4))
      NEW(11)=1D0
      END

      SUBROUTINE WCM_DAMAGE_B(P,OLD,E,EFF,FI,DT,LC,NEW,IERR)
      IMPLICIT NONE
      INTEGER I,IERR
      DOUBLE PRECISION P(*),OLD(40),E(6),EFF(6),FI(4),DT,LC
      DOUBLE PRECISION NEW(40),DEL(4),SIG(4),DFINAL,KAP,ETA,DNEW
      IERR=0
C--- B: equivalent displacement, linear softening and viscosity.
      CALL WCM_EQ(E,EFF,LC,DEL,SIG)
      DO I=1,4
          ETA=P(24)
          IF(I.GE.3) ETA=P(25)
          IF(FI(I).GE.1D0) NEW(I+8)=1D0
          IF(NEW(I+8).GE.1D0) THEN
              IF(NEW(I+12).LE.0D0) THEN
                  IF(DEL(I).LE.0D0.OR.SIG(I).LE.0D0) THEN
                      IERR=20+I
                      RETURN
                  END IF
                  NEW(I+12)=DEL(I)
                  NEW(I+16)=SIG(I)
              END IF
              DFINAL=2D0*P(I+16)/NEW(I+16)
              IF(DFINAL.LE.NEW(I+12)) THEN
                  IERR=30+I
                  RETURN
              END IF
              KAP=MAX(OLD(I+4),DEL(I))
              NEW(I+4)=KAP
              DNEW=OLD(I)
              IF(KAP.GT.NEW(I+12)) DNEW=DFINAL
     1          *(KAP-NEW(I+12))/(KAP*(DFINAL-NEW(I+12)))
              NEW(I)=MIN(MAX(OLD(I),DNEW,0D0),1D0)
          END IF
          IF(ETA.GT.0D0.AND.DT.GT.0D0) THEN
              NEW(I+20)=(ETA*OLD(I+20)+DT*NEW(I))
     1                  /(ETA+DT)
          ELSE
              NEW(I+20)=NEW(I)
          END IF
          NEW(I+20)=MIN(MAX(NEW(I+20),0D0),1D0)
      END DO
      END
