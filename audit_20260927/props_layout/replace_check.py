from pathlib import Path
p=Path('wcm_hashin_umat/src/wcm_core.for');s=p.read_text(encoding='utf-8');a=s.index('      SUBROUTINE WCM_CHECK(');b=s.index('      SUBROUTINE WCM_ROTATE(',a)
new='''      SUBROUTINE WCM_CHECK(MODE,PROPS,NPROPS,OLD,DT,LC,IERR)
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

'''
s=s[:a]+new+s[b:];p.write_text(s,encoding='utf-8')
