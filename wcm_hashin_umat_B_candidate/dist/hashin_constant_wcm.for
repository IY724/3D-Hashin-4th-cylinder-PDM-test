C=======================================================================
C A（直接刚度折减） UMAT: local +/-theta plies, 3D Hashin damage.
C Material slots are read and checked by WCM_CHECK below.
C Explicit material-card order; the Fortran indices are 1-based:
C E11=PROPS(1)
C E22=PROPS(2)
C E33=PROPS(3)
C PR12=PROPS(4)
C PR13=PROPS(5)
C PR23=PROPS(6)
C G12=PROPS(7)
C G13=PROPS(8)
C G23=PROPS(9)
C XT=PROPS(10)
C XC=PROPS(11)
C YT=PROPS(12)
C YC=PROPS(13)
C S12=PROPS(14)
C S13=PROPS(15)
C S23=PROPS(16)
C DFT0=PROPS(17)
C DFC0=PROPS(18)
C DMT0=PROPS(19)
C DMC0=PROPS(20)
C SMT=PROPS(21)
C SMC=PROPS(22)
C THETA=PROPS(23)
C WPLUS=PROPS(24)
C UPDATE=PROPS(25)
C Units: modulus/strength MPa; GFT/GFC/GMT/GMC N/mm;
C THETA deg; ETA_F/ETA_M step-time units.
C S13=PROPS(15) is reserved; current tau13 uses S12.
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
      PARAMETER (MODE=1)
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
C--- B softening endpoint error: request a smaller increment.
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
C--- NLGEOM tangent correction.
      IF(JSTEP(3).EQ.1) THEN
          DO J=1,3
              DDSDDE(:,J)=DDSDDE(:,J)+STRESS
          END DO
      END IF
C--- Save damage histories and numbered SDV results.
      IF(IADV.EQ.1) THEN
          CALL WCM_PACK(MODE,PROPS,NPROPS,OLD,NEW,STATEV,NSTATV)
          SSE=0.5D0*DOT_PRODUCT(STRESS,STRAN+DSTRAN)
C--- Signed stress-work diagnostic.
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
      SUBROUTINE WCM_PACK(MODE,PROPS,NP,OLD,NEW,STATE,NS)
      IMPLICIT NONE
      INTEGER MODE,NP,NS,K,IB,IUSED
      DOUBLE PRECISION PROPS(NP),OLD(84),NEW(84),STATE(NS),W
      DOUBLE PRECISION D(4),V(4),U(4),DF,DM,UF,UM
      DOUBLE PRECISION WPLUS,UPDATE,CAPF,CAPM
      DOUBLE PRECISION DFT0,DFC0,DMT0,DMC0,DLOSS(4)
      STATE=0D0
      STATE(1:4)=-1D100
      IUSED=7
      IF(MODE.EQ.1) THEN
          DFT0=PROPS(17)
          DFC0=PROPS(18)
          DMT0=PROPS(19)
          DMC0=PROPS(20)
          DLOSS=(/DFT0,DFC0,DMT0,DMC0/)
          WPLUS=PROPS(24)
          UPDATE=PROPS(25)
          STATE(9:12)=NEW(7:10)
          STATE(13:16)=NEW(47:50)
      ELSE
          WPLUS=PROPS(27)
          CAPF=PROPS(28)
          CAPM=PROPS(29)
          UPDATE=PROPS(30)
          IUSED=9
          STATE(11:34)=NEW(1:24)
          STATE(35:58)=NEW(41:64)
      END IF
      STATE(NS-3:NS)=NEW(81:84)
      DO K=1,2
          IB=40*(K-1)
          W=WPLUS
          IF(K.EQ.2) W=1D0-W
          IF(W.EQ.0D0) CYCLE
          STATE(1:4)=MAX(STATE(1:4),NEW(IB+25:IB+28))
          IF(MODE.EQ.1) THEN
              D=NEW(IB+7:IB+10)
              U=OLD(IB+7:IB+10)*DLOSS
              IF(NINT(UPDATE).EQ.1) U=D*DLOSS
          ELSE
              D=NEW(IB+1:IB+4)
              V=NEW(IB+21:IB+24)
              STATE(7)=MAX(STATE(7),
     1                    1D0-(1D0-V(1))*(1D0-V(2)))
              STATE(8)=MAX(STATE(8),
     1                    1D0-(1D0-V(3))*(1D0-V(4)))
              U=OLD(IB+21:IB+24)
              IF(NINT(UPDATE).EQ.1) U=V
          END IF
          DF=1D0-(1D0-D(1))*(1D0-D(2))
          DM=1D0-(1D0-D(3))*(1D0-D(4))
          UF=1D0-(1D0-U(1))*(1D0-U(2))
          UM=1D0-(1D0-U(3))*(1D0-U(4))
          IF(MODE.EQ.2) THEN
              UF=MIN(UF,CAPF)
              UM=MIN(UM,CAPM)
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
      SUBROUTINE WCM_UPDATE(MODE,PROPS,NP,OLD,E0,DE,DT,LC,IADV,
     1                     NEW,S,C,IERR)
      IMPLICIT NONE
      INTEGER MODE,NP,IADV,IERR,I
      DOUBLE PRECISION PROPS(NP),OLD(84),NEW(84),E0(6),DE(6)
      DOUBLE PRECISION DT,LC,S(6),C(6,6),C0(6,6),E1(6),UPDATE
      IERR=0
      NEW=OLD
      S=0D0
      C=0D0

C--- Material card and undamaged stiffness.
      CALL WCM_CHECK(MODE,PROPS,NP,OLD,DT,LC,IERR)
      IF(IERR.NE.0) RETURN
      CALL WCM_ELASTIC(PROPS,C0,IERR)
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
C     E0/E1 are passed explicitly so B can locate the first onset along
C     the real increment path instead of trusting the endpoint (plan 6.1).
      CALL WCM_EVAL(MODE,PROPS,NP,C0,OLD,E0,E1,DT,LC,IADV,1,
     1              NEW,S,C,IERR)
      IF(IERR.NE.0) RETURN

C--- B current feedback: DDSDDE from the same OLD history.
      IF(MODE.EQ.2.AND.IADV.EQ.1) THEN
          UPDATE=PROPS(30)
          IF(NINT(UPDATE).EQ.1)
     1      CALL WCM_TANGENT_B(PROPS,NP,C0,OLD,E0,E1,DT,LC,C,IERR)
      END IF
      END

      SUBROUTINE WCM_TANGENT_B(PROPS,NP,C0,OLD,E0,E,DT,LC,C,IERR)
      IMPLICIT NONE
      INTEGER NP,IERR,I,J
      DOUBLE PRECISION PROPS(NP),C0(6,6),OLD(84),E0(6),E(6),DT,LC,C(6,6)
      DOUBLE PRECISION SP(6),SM(6),EP(6),EM(6),H
      DOUBLE PRECISION TMP(84),CT(6,6)
      DO J=1,6
          H=MAX(1D-9,1D-6*ABS(E(J)))
          EP=E
          EM=E
          EP(J)=EP(J)+H
          EM(J)=EM(J)-H
          CALL WCM_EVAL(2,PROPS,NP,C0,OLD,E0,EP,DT,LC,1,0,
     1                  TMP,SP,CT,IERR)
          IF(IERR.NE.0) RETURN
          CALL WCM_EVAL(2,PROPS,NP,C0,OLD,E0,EM,DT,LC,1,0,
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

      SUBROUTINE WCM_ELASTIC(PROPS,C,IERR)
      IMPLICIT NONE
      DOUBLE PRECISION PROPS(*),C(6,6),V21,V31,V32,DEN,FAC
      DOUBLE PRECISION E11,E22,E33,PR12,PR13,PR23,G12,G13,G23
      INTEGER IERR
      E11=PROPS(1)
      E22=PROPS(2)
      E33=PROPS(3)
      PR12=PROPS(4)
      PR13=PROPS(5)
      PR23=PROPS(6)
      G12=PROPS(7)
      G13=PROPS(8)
      G23=PROPS(9)
      C=0D0
      IERR=0
      V21=E22*PR12/E11
      V31=E33*PR13/E11
      V32=E33*PR23/E22
      DEN=1D0-PR12*V21-PR13*V31-PR23*V32
     1    -2D0*V21*V32*PR13
      IF (1D0-PR12*V21.LE.1D-12.OR.DEN.LE.1D-12) THEN
          IERR=12
          RETURN
      END IF
      FAC=1D0/DEN
      C(1,1)=E11*FAC*(1D0-PR23*V32)
      C(2,2)=E22*FAC*(1D0-PR13*V31)
      C(3,3)=E33*FAC*(1D0-PR12*V21)
      C(1,2)=E11*FAC*(V21+V31*PR23)
      C(1,3)=E11*FAC*(V31+V21*V32)
      C(2,3)=E22*FAC*(V32+PR12*V31)
      C(2,1)=C(1,2)
      C(3,1)=C(1,3)
      C(3,2)=C(2,3)
      C(4,4)=G12
      C(5,5)=G13
      C(6,6)=G23
      END

      SUBROUTINE WCM_STIFF(MODE,PROPS,D,C0,C)
      IMPLICIT NONE
      INTEGER MODE,I,J
      DOUBLE PRECISION PROPS(*),D(4),C0(6,6),C(6,6)
      DOUBLE PRECISION DF,DM,RF,RM,RS,RMT,RMC,RN,RMIN
      DOUBLE PRECISION SMT,SMC,CAPF,CAPM
      PARAMETER (RMIN=1D-6)
      DF=1D0-(1D0-D(1))*(1D0-D(2))
      DM=1D0-(1D0-D(3))*(1D0-D(4))
      IF (MODE.EQ.1) THEN
C--- A fixed degradation factors.
          SMT=PROPS(21)
          SMC=PROPS(22)
          RF=MAX(1D0-DF,1D-6)
          RM=MAX(1D0-DM,1D-6)
          RMT=MAX(1D0-SMT*D(3),1D-6)
          RMC=MAX(1D0-SMC*D(4),1D-6)
          RS=RF*RMT*RMC
      ELSE
C--- B stiffness caps and viscous damage.
          SMT=PROPS(22)
          SMC=PROPS(23)
          CAPF=PROPS(28)
          CAPM=PROPS(29)
          DF=MIN(MAX(DF,0D0),CAPF)
          DM=MIN(MAX(DM,0D0),CAPM)
          RF=MAX(1D0-DF,RMIN)
          RM=MAX(1D0-DM,RMIN)
          RS=(1D0-SMT*D(3))*(1D0-SMC*D(4))
          RS=MIN(MAX(RS,RMIN),1D0)
      END IF
C--- Positive residual stiffness for the tangent.
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

      SUBROUTINE WCM_HASHIN(S,PROPS,FI)
      IMPLICIT NONE
      DOUBLE PRECISION S(6),PROPS(*),FI(4),Q
      DOUBLE PRECISION XT,XC,YT,YC,S12,S23
C--- 3D Hashin indices; tau13 uses S12 (P14).
      XT=PROPS(10)
      XC=PROPS(11)
      YT=PROPS(12)
      YC=PROPS(13)
      S12=PROPS(14)
      S23=PROPS(16)
      FI=0D0
      IF(S(1).GE.0D0) THEN
          FI(1)=(S(1)/XT)**2+(S(4)/S12)**2
     1         +(S(5)/S12)**2
      ELSE
          FI(2)=(-S(1)/XC)**2
      END IF
      Q=S(2)+S(3)
      IF(Q.GE.0D0) THEN
          FI(3)=(Q/YT)**2+(S(6)**2-S(2)*S(3))/S23**2
     1         +(S(4)/S12)**2+(S(5)/S12)**2
      ELSE
          FI(4)=((YC/(2D0*S23))**2-1D0)*Q/YC
     1         +Q**2/(2D0*S23)**2
     2         +(S(6)**2-S(2)*S(3))/S23**2
     3         +(S(4)/S12)**2+(S(5)/S12)**2
      END IF
      END

      SUBROUTINE WCM_EQ(PROPS,E,S,LC,DEL,SIG)
      IMPLICIT NONE
      DOUBLE PRECISION PROPS(*),E(6),S(6),LC,DEL(4),SIG(4),SH
      DOUBLE PRECISION E11
C--- Equivalent displacement and stress for B softening.
      E11=PROPS(1)
      DEL(1)=LC*SQRT(MAX(E(1),0D0)**2+E(4)**2+E(5)**2)
C--- FC is stress-triggered. Use the matching stress-equivalent strain:
C    the previous -E11 measure could be zero under Poisson coupling
C    even while S11<-XC, leading to unseedable error 1022.
      DEL(2)=LC*MAX(-S(1),0D0)/E11
      SH=E(4)**2+E(5)**2+E(6)**2
      DEL(3)=LC*SQRT(MAX(E(2),0D0)**2+MAX(E(3),0D0)**2+SH)
      DEL(4)=LC*SQRT(MAX(-E(2),0D0)**2+MAX(-E(3),0D0)**2+SH)
      SIG=0D0
      SH=S(4)*E(4)+S(5)*E(5)+S(6)*E(6)
      IF(DEL(1).GT.0D0) SIG(1)=LC*(MAX(S(1),0D0)
     1  *MAX(E(1),0D0)+S(4)*E(4)+S(5)*E(5))/DEL(1)
      IF(DEL(2).GT.0D0) SIG(2)=MAX(-S(1),0D0)
      IF(DEL(3).GT.0D0) SIG(3)=LC*(MAX(S(2),0D0)
     1  *MAX(E(2),0D0)+MAX(S(3),0D0)*MAX(E(3),0D0)+SH)/DEL(3)
      IF(DEL(4).GT.0D0) SIG(4)=LC*(MAX(-S(2),0D0)
     1  *MAX(-E(2),0D0)+MAX(-S(3),0D0)*MAX(-E(3),0D0)+SH)/DEL(4)
      END

      SUBROUTINE WCM_EVAL(MODE,PROPS,NP,C0,OLD,E0,E,DT,LC,IADV,IONSET,
     1                    NEW,S,C,IERR)
      IMPLICIT NONE
      INTEGER MODE,NP,IADV,IONSET,IERR,K,IB
      DOUBLE PRECISION PROPS(NP),C0(6,6),OLD(84),E0(6),E(6),DT,LC
      DOUBLE PRECISION NEW(84),S(6),C(6,6),B(6,6),EL(6),EL0(6)
      DOUBLE PRECISION SL(6),CL(6,6),LOCAL(40),TH,W,WP,THETA
      NEW=OLD
      S=0D0
      C=0D0
      IERR=0
      IF(MODE.EQ.1) THEN
          THETA=PROPS(23)
          WP=PROPS(24)
      ELSE
          THETA=PROPS(26)
          WP=PROPS(27)
      END IF
      DO K=1,2
          IB=40*(K-1)
          TH=THETA
          W=WP
          IF(K.EQ.2) THEN
              TH=-TH
              W=1D0-WP
          END IF
C--- Skip a zero-weight family.
          IF(W.EQ.0D0) CYCLE
C--- W axes to the local +/-theta fiber axes.
          CALL WCM_ROTATE(TH,B)
          EL=MATMUL(B,E)
          EL0=MATMUL(B,E0)
C--- Local stress, Hashin criterion and damage update.
          CALL WCM_PLY(MODE,PROPS,C0,OLD(IB+1:IB+40),EL0,EL,DT,LC,
     1                 IADV,IONSET,LOCAL,SL,CL,IERR)
          IF(IERR.NE.0) THEN
C--- Encode which angle family failed.
              IERR=IERR+1000*K
              RETURN
          END IF
          NEW(IB+1:IB+40)=LOCAL
C--- Rotate stress and tangent back; sum the two families.
          S=S+W*MATMUL(TRANSPOSE(B),SL)
          C=C+W*MATMUL(TRANSPOSE(B),MATMUL(CL,B))
      END DO
      IF(IADV.EQ.1) THEN
          NEW(81)=THETA
          NEW(82)=WP
          NEW(83)=DBLE(MODE)
          NEW(84)=202609D0
      END IF
      END

      SUBROUTINE WCM_PLY(MODE,PROPS,C0,OLD,E0,E,DT,LC,IADV,IONSET,
     1                   NEW,S,C,IERR)
      IMPLICIT NONE
      INTEGER MODE,IADV,IONSET,IERR
      DOUBLE PRECISION PROPS(*),C0(6,6),OLD(40),E0(6),E(6),DT,LC
      DOUBLE PRECISION NEW(40),S(6),C(6,6),D(4),FI(4),EFF(6),EFF0(6)
      DOUBLE PRECISION DFT0,DFC0,DMT0,DMC0,DLOSS(4),UPDATE
      IERR=0
      NEW=OLD
      D=0D0
      IF(MODE.EQ.1) THEN
          DFT0=PROPS(17)
          DFC0=PROPS(18)
          DMT0=PROPS(19)
          DMC0=PROPS(20)
          UPDATE=PROPS(25)
          DLOSS=(/DFT0,DFC0,DMT0,DMC0/)
          D=OLD(7:10)*DLOSS
      ELSE
          UPDATE=PROPS(30)
          D=OLD(21:24)
      END IF

C--- Trial stress and frozen-damage tangent.
      CALL WCM_STIFF(MODE,PROPS,D,C0,C)
      S=MATMUL(C,E)
      EFF=S
      IF(MODE.EQ.2) EFF=MATMUL(C0,E)

C--- Four Hashin initiation indices.
      CALL WCM_HASHIN(EFF,PROPS,FI)
      IF(IADV.EQ.1) THEN
          IF(MODE.EQ.1) THEN
              CALL WCM_DAMAGE_A(PROPS,FI,NEW)
              D=NEW(7:10)*DLOSS
          ELSE
C--- B: locate the first onset along the real increment path.
              EFF0=MATMUL(C0,E0)
              CALL WCM_DAMAGE_B(PROPS,OLD,E0,E,EFF0,EFF,FI,DT,LC,IONSET,
     1                              NEW,IERR)
              IF(IERR.NE.0) RETURN
              D=NEW(21:24)
          END IF
C--- Current-damage feedback if UPDATE=1.
          IF(NINT(UPDATE).EQ.1) THEN
              CALL WCM_STIFF(MODE,PROPS,D,C0,C)
              S=MATMUL(C,E)
          END IF
          NEW(25:28)=FI
          NEW(29:34)=S
          NEW(35:40)=E
      END IF
      END

      SUBROUTINE WCM_DAMAGE_A(PROPS,FI,NEW)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION PROPS(*),FI(4),NEW(40),D(4)
      DOUBLE PRECISION DFT0,DFC0,DMT0,DMC0,DLOSS(4)
C--- A: latch failure and apply fixed stiffness factors.
      DFT0=PROPS(17)
      DFC0=PROPS(18)
      DMT0=PROPS(19)
      DMC0=PROPS(20)
      DLOSS=(/DFT0,DFC0,DMT0,DMC0/)
      DO I=1,4
          IF(FI(I).GE.1D0) NEW(I+6)=1D0
      END DO
      D=NEW(7:10)*DLOSS
      NEW(1:4)=FI
      NEW(5)=1D0-(1D0-D(1))*(1D0-D(2))
      NEW(6)=1D0-(1D0-D(3))*(1D0-D(4))
      NEW(11)=1D0
      END

      SUBROUTINE WCM_DAMAGE_B(PROPS,OLD,E0,E,EFF0,EFF,FI,DT,LC,IONSET,
     1                        NEW,IERR)
      IMPLICIT NONE
      INTEGER I,IERR,IONSET
      DOUBLE PRECISION PROPS(*),OLD(40),E0(6),E(6),EFF0(6),EFF(6),FI(4)
      DOUBLE PRECISION DT,LC
      DOUBLE PRECISION NEW(40),DEL(4),SIG(4),DFINAL,KAP,ETA,DNEW
      DOUBLE PRECISION GS,EPSF
      DOUBLE PRECISION GFT,GFC,GMT,GMC,GCRIT(4),ETA_F,ETA_M
      LOGICAL GV
      IERR=0
C--- B: equivalent displacement, linear softening and viscosity.
C--- Endpoint equivalents drive evolution only; onset quantities come
C    from the path search below (plan 6.2), so a crossing that happens
C    inside the increment is no longer seeded with endpoint values.
      GFT=PROPS(17)
      GFC=PROPS(18)
      GMT=PROPS(19)
      GMC=PROPS(20)
      GCRIT=(/GFT,GFC,GMT,GMC/)
      ETA_F=PROPS(24)
      ETA_M=PROPS(25)
      CALL WCM_EQ(PROPS,E,EFF,LC,DEL,SIG)
      DO I=1,4
          ETA=ETA_F
          IF(I.GE.3) ETA=ETA_M
          IF(NEW(I+8).LT.1D0.AND.IONSET.EQ.1) THEN
              CALL WCM_ONSET_B(PROPS,OLD,E0,E,EFF0,EFF,I,LC,NEW,IERR)
              IF(IERR.NE.0) RETURN
C--- Endpoint fallback only when the endpoint branch is meaningful:
C    with FE rounding the tension/compression sign may be pure noise.
              GS=EFF(1)
              IF(I.GE.3) GS=EFF(2)+EFF(3)
              EPSF=1D-10*MAX(1D0,MAX(ABS(EFF(1)),ABS(EFF(2)),
     1            ABS(EFF(3)),ABS(EFF(4)),ABS(EFF(5)),ABS(EFF(6))))
              GV=.FALSE.
              IF(I.EQ.1.OR.I.EQ.3) GV=GS.GE.-EPSF
              IF(I.EQ.2.OR.I.EQ.4) GV=GS.LT.-EPSF
              IF(NEW(I+8).LT.1D0.AND.FI(I).GE.1D0.AND.GV) THEN
C--- Endpoint criterion met but no interior onset located (rounding
C    edge): seed at the endpoint with strict validity checks.
                  IF(DEL(I).LE.0D0.OR.SIG(I).LE.0D0) THEN
                      IERR=20+I
                      RETURN
                  END IF
                  NEW(I+8)=1D0
                  NEW(I+12)=DEL(I)
                  NEW(I+16)=SIG(I)
              END IF
          END IF
          IF(NEW(I+8).GE.1D0) THEN
              IF(NEW(I+12).LE.0D0) THEN
C--- Latched without onset quantities: corrupt or foreign state.
                  IERR=20+I
                  RETURN
              END IF
              DFINAL=2D0*GCRIT(I)/NEW(I+16)
              KAP=MAX(OLD(I+4),NEW(I+12),DEL(I))
              NEW(I+4)=KAP
              DNEW=OLD(I)
              IF(DFINAL.LE.NEW(I+12)) THEN
C--- User policy: an incompatible MATRIX mode fails locally.
C    This limiting branch is not fracture-energy-exact softening.
C    Only this mode completes; fiber history remains independent.
C    Viscosity below and UPDATE retain their existing meanings.
                  IF(I.GE.3) THEN
                      DNEW=1D0
                  ELSE
C--- Fiber incompatibility still requires a valid softening law.
                      WRITE(7,*) 'WCM_FIBER_SOFTEN mode/G/LC: ',
     1                    I,GCRIT(I),LC
                      WRITE(7,*) 'WCM_FIBER_SOFTEN d0/s0/df: ',
     1                    NEW(I+12),NEW(I+16),DFINAL
                      IERR=30+I
                      RETURN
                  END IF
              ELSE
C--- Compatible modes retain the original linear softening law.
                  IF(KAP.GT.NEW(I+12)) DNEW=DFINAL
     1              *(KAP-NEW(I+12))/(KAP*(DFINAL-NEW(I+12)))
              END IF
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

C=======================================================================
C B onset location along the real increment path (plan 6.2).
C Path: E(s)=E0+s*(E-E0), EFF(s)=EFF0+s*(EFF-EFF0), s in [0,1].
C For an uninitiated mode I, find the earliest s with FI_i(s)>=1 on its
C active sign branch and DEL_i(s)>0, SIG_i(s)>0 there; branch entries
C are used as-is (no F=1 fabrication). Roots of the quadratic F_i(s)-1
C are enumerated with a stability guard, deduped, slope-filtered to
C crossings into FI>=1 and bisection-polished when bracketed.
C Errors: 40+I criterion already met at E0 without history (plan 6.2.9),
C 20+I crossing exists but no seedable onset point.
      SUBROUTINE WCM_ONSET_B(PROPS,OLD,E0,E,EFF0,EFF,I,LC,NEW,IERR)
      IMPLICIT NONE
      INTEGER I,IERR,NB,NS,K,M,NC,IT
      DOUBLE PRECISION PROPS(*),OLD(40),E0(6),E(6),EFF0(6)
      DOUBLE PRECISION EFF(6),NEW(40)
      DOUBLE PRECISION LC,F0(4),X0(6),DX(6),G0,GD,SS,BP(3),A,B,GA
      DOUBLE PRECISION SMAG,EPSG
      DOUBLE PRECISION A2,A1,A0,SC,FA,DISC,SQ,RQ,RL(2),RT,DER
      DOUBLE PRECISION EP(6),SP(6),DUM(4),SG(4),CAND(2),S,U,V
      DOUBLE PRECISION PU,PV,PW,XM,H,FVAL,GC
      LOGICAL SOMEINV
      LOGICAL SKIP40
      IERR=0
C--- Effective stress is linear along the path: X(s)=X0+s*DX.
      DO K=1,6
          X0(K)=EFF0(K)
          DX(K)=EFF(K)-EFF0(K)
      END DO
C--- Noise scale: FE rounding leaves O(1e-19) strain noise, so stress
C    components below ~1e-10 of the path scale carry no sign
C    information (plan 7: pure-shear q drifts across zero per IP).
      SMAG=0D0
      DO K=1,6
          SMAG=MAX(SMAG,ABS(X0(K)),ABS(DX(K)))
      END DO
      EPSG=1D-10*MAX(SMAG,1D0)
C--- Sign-branch switch: EFF11 for fiber modes, EFF22+EFF33 for matrix
C    modes (plan 6.2.3). In the noise band the WCM_HASHIN convention
C    (Q>=0, S1>=0 -> tension branch) is enforced: the tension mode owns
C    the path and the compression mode cannot start; this keeps every
C    integration point deterministic.
      IF(I.LE.2) THEN
          G0=X0(1)
          GD=DX(1)
      ELSE
          G0=X0(2)+X0(3)
          GD=DX(2)+DX(3)
      END IF
      IF(ABS(G0).LE.EPSG.AND.ABS(GD).LE.EPSG
     1   .AND.(I.EQ.2.OR.I.EQ.4)) THEN
          WRITE(99,*) 'ONS DEGEN-RET I=',I
          RETURN
      END IF
C--- An uninitiated mode already MEANINGFULLY over the surface at E0
C    has no locatable onset history (plan 6.2.9); within the noise band
C    the search below seeds at s=0 instead of inventing history.
      CALL WCM_HASHIN(EFF0,PROPS,F0)
      IF(F0(I).GE.1D0) THEN
          SKIP40=(I.EQ.2.OR.I.EQ.4).AND.G0.GE.-EPSG
          IF(.NOT.SKIP40) THEN
              WRITE(99,*) 'ONS 40 I=',I
              IERR=40+I
              RETURN
          END IF
      END IF
      BP(1)=0D0
      BP(2)=1D0
      NB=2
      IF(ABS(G0).GT.EPSG.AND.ABS(GD).GT.EPSG) THEN
          SS=-G0/GD
          IF(SS.GT.1D-12.AND.SS.LT.1D0-1D-12) THEN
              NB=3
              BP(3)=SS
          END IF
      END IF
      DO K=1,NB-1
          DO M=K+1,NB
              IF(BP(M).LT.BP(K)) THEN
                  SS=BP(K)
                  BP(K)=BP(M)
                  BP(M)=SS
              END IF
          END DO
      END DO
      SOMEINV=.FALSE.
      DO NS=1,NB-1
          A=BP(NS)
          B=BP(NS+1)
          GA=G0+0.5D0*(A+B)*GD
C--- Skip intervals where this mode's sign branch is not active; the
C    boundary band +/-EPSG counts as active for both tension modes
C    (FI is continuous across the branch boundary).
          IF(I.EQ.1.OR.I.EQ.3) THEN
              IF(GA.LT.-EPSG) CYCLE
          ELSE
              IF(GA.GT.EPSG) CYCLE
          END IF
          CALL WCM_FICOEF(I,PROPS,X0,DX,A2,A1,A0)
          SC=MAX(ABS(A2),ABS(A1),ABS(A0))
          IF(SC.LE.0D0) CYCLE
C--- Coefficient scale-degradation guard (plan 6.2.6).
          IF(ABS(A2).LE.1D-14*SC) A2=0D0
          NC=0
          FA=A2*A*A+A1*A+A0
          IF(FA.GE.0D0) THEN
C--- Branch entry: over the surface as soon as the branch opens.
              NC=1
              CAND(1)=A
              WRITE(99,*) 'ONS BRENTRY I=',I,' A=',A,' FA=',FA
          ELSE IF(A2.NE.0D0) THEN
              DISC=A1*A1-4D0*A2*A0
              IF(DISC.GE.0D0) THEN
                  SQ=SQRT(DISC)
                  RQ=-0.5D0*(A1+SIGN(SQ,A1))
                  IF(ABS(RQ).GT.0D0) THEN
                      RL(1)=RQ/A2
                      RL(2)=A0/RQ
                      IF(ABS(RL(1)-RL(2)).LE.1D-12) THEN
                          NC=1
                          CAND(1)=0.5D0*(RL(1)+RL(2))
                      ELSE
                          NC=2
                          CAND(1)=MIN(RL(1),RL(2))
                          CAND(2)=MAX(RL(1),RL(2))
                      END IF
                  END IF
              END IF
          ELSE IF(A1.NE.0D0) THEN
              NC=1
              CAND(1)=-A0/A1
          END IF
          DO M=1,NC
              SS=CAND(M)
              IF(SS.LT.A-1D-12.OR.SS.GT.B+1D-12) CYCLE
              SS=MIN(MAX(SS,A),B)
C--- The candidate must sit meaningfully inside this mode's own sign
C    branch: the shared shear surface at the tension/compression
C    boundary belongs to the tension modes only, otherwise a tangent
C    perturbation would let MC seed on it and fail deltaf<=delta0.
              GC=G0+SS*GD
              IF(I.EQ.1.OR.I.EQ.3) THEN
                  IF(GC.LT.-EPSG) CYCLE
              ELSE
                  IF(GC.GT.-EPSG) CYCLE
              END IF
              IF(FA.LT.0D0) THEN
C--- Keep only crossings into FI>=1 (non-negative slope at the root;
C    tangent double roots count as conservative candidates).
                  IF(A2.NE.0D0) THEN
                      DER=A1+2D0*A2*SS
                      IF(DER.LT.0D0) CYCLE
                  ELSE IF(A1.LT.0D0) THEN
                      CYCLE
                  END IF
                  FVAL=A2*SS*SS+A1*SS+A0+1D0
                  IF(ABS(FVAL-1D0)/MAX(1D0,ABS(FVAL)).GT.1D-8) THEN
C--- Bisection polish when a sign bracket exists (plan 6.2.6).
                      H=(B-A)*1D-3
                      U=MAX(A,SS-H)
                      V=MIN(B,SS+H)
                      PU=A2*U*U+A1*U+A0
                      PV=A2*V*V+A1*V+A0
                      IF(PU*PV.LE.0D0) THEN
                          DO IT=1,80
                              XM=0.5D0*(U+V)
                              PW=A2*XM*XM+A1*XM+A0
                              IF(PU*PW.LE.0D0) THEN
                                  V=XM
                                  PV=PW
                              ELSE
                                  U=XM
                                  PU=PW
                              END IF
                          END DO
                          SS=0.5D0*(U+V)
                      END IF
                      FVAL=A2*SS*SS+A1*SS+A0+1D0
                      IF(ABS(FVAL-1D0)/MAX(1D0,ABS(FVAL)).GT.1D-8)
     1                    CYCLE
                  END IF
              END IF
C--- Seed with the original WCM_EQ at the onset point (plan 6.2.8).
              DO K=1,6
                  EP(K)=E0(K)+SS*(E(K)-E0(K))
                  SP(K)=X0(K)+SS*DX(K)
              END DO
              CALL WCM_EQ(PROPS,EP,SP,LC,DUM,SG)
              IF(DUM(I).GT.0D0.AND.SG(I).GT.0D0) THEN
                  NEW(I+8)=1D0
                  NEW(I+12)=DUM(I)
                  NEW(I+16)=SG(I)
                  WRITE(99,*) 'ONS SEED I=',I,' SS=',SS,' GC=',GC,
     1                ' D0=',DUM(I),' SG0=',SG(I)
                  RETURN
              END IF
              SOMEINV=.TRUE.
          END DO
      END DO
      IF(SOMEINV) IERR=20+I
      END

C=======================================================================
C Quadratic coefficients of FI_i(s)-1 for the linear effective-stress
C path X(s)=X0+s*DX; coefficients cross-checked against WCM_HASHIN.
      SUBROUTINE WCM_FICOEF(I,PROPS,X0,DX,A2,A1,A0)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION PROPS(*),X0(6),DX(6),A2,A1,A0
      DOUBLE PRECISION XT,XC,YT,YC,S12,S23
      DOUBLE PRECISION X1,D1,X2,D2,X3,D3,X4,D4,X5,D5,X6,D6,Q0,QD
      XT=PROPS(10)
      XC=PROPS(11)
      YT=PROPS(12)
      YC=PROPS(13)
      S12=PROPS(14)
      S23=PROPS(16)
      X1=X0(1)
      D1=DX(1)
      X2=X0(2)
      D2=DX(2)
      X3=X0(3)
      D3=DX(3)
      X4=X0(4)
      D4=DX(4)
      X5=X0(5)
      D5=DX(5)
      X6=X0(6)
      D6=DX(6)
      Q0=X2+X3
      QD=D2+D3
      A2=0D0
      A1=0D0
      A0=-1D0
      IF(I.EQ.1) THEN
          A2=(D1/XT)**2+(D4/S12)**2+(D5/S12)**2
          A1=2D0*(X1*D1/XT**2+X4*D4/S12**2+X5*D5/S12**2)
          A0=(X1/XT)**2+(X4/S12)**2+(X5/S12)**2-1D0
      ELSE IF(I.EQ.2) THEN
          A2=(D1/XC)**2
          A1=2D0*X1*D1/XC**2
          A0=(X1/XC)**2-1D0
      ELSE IF(I.EQ.3) THEN
          A2=(QD/YT)**2+(D6*D6-D2*D3)/S23**2
     1       +(D4/S12)**2+(D5/S12)**2
          A1=2D0*Q0*QD/YT**2+(2D0*X6*D6-X2*D3-X3*D2)/S23**2
     1       +2D0*(X4*D4+X5*D5)/S12**2
          A0=(Q0/YT)**2+(X6*X6-X2*X3)/S23**2+(X4/S12)**2
     1       +(X5/S12)**2-1D0
      ELSE
          A2=(QD/(2D0*S23))**2+(D6*D6-D2*D3)/S23**2
     1       +(D4/S12)**2+(D5/S12)**2
          A1=QD*((YC/(2D0*S23))**2-1D0)/YC+2D0*Q0*QD/(2D0*S23)**2
     1       +(2D0*X6*D6-X2*D3-X3*D2)/S23**2
     1       +2D0*(X4*D4+X5*D5)/S12**2
          A0=Q0*((YC/(2D0*S23))**2-1D0)/YC+(Q0/(2D0*S23))**2
     1       +(X6*X6-X2*X3)/S23**2+(X4/S12)**2
     1       +(X5/S12)**2-1D0
      END IF
      END
