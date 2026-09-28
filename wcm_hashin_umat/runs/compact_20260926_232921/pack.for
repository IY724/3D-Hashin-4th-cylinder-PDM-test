C Persistent history is separate for the two weighted families.
      SUBROUTINE WCM_UNPACK(MODE,STATE,NS,OLD)
      IMPLICIT NONE
      INTEGER MODE,NS
      DOUBLE PRECISION STATE(NS),OLD(84)
      OLD=0D0
      IF(MODE.EQ.1) THEN
          OLD(7:10)=STATE(7:10)
          OLD(47:50)=STATE(11:14)
      ELSE
          OLD(1:24)=STATE(9:32)
          OLD(41:64)=STATE(33:56)
      END IF
      OLD(81:84)=STATE(NS-3:NS)
      END

C Output: FI(FT,FC,MT,MC), combined DF,DM; B adds viscous DF,DM.
C Combine tension/compression within each family BEFORE maximum.
      SUBROUTINE WCM_PACK(MODE,P,NP,NEW,STATE,NS)
      IMPLICIT NONE
      INTEGER MODE,NP,NS,K,IB,IT
      DOUBLE PRECISION P(NP),NEW(84),STATE(NS),W,D(4),V(4)
      DOUBLE PRECISION DF,DM,VF,VM
      STATE=0D0
      IT=23
      IF(MODE.EQ.1) THEN
          STATE(7:10)=NEW(7:10)
          STATE(11:14)=NEW(47:50)
      ELSE
          IT=26
          STATE(9:32)=NEW(1:24)
          STATE(33:56)=NEW(41:64)
      END IF
      STATE(NS-3:NS)=NEW(81:84)
      DO K=1,2
          IB=40*(K-1)
          W=P(IT+1)
          IF(K.EQ.2) W=1D0-W
          IF(W.EQ.0D0) CYCLE
          STATE(1:4)=MAX(STATE(1:4),NEW(IB+25:IB+28))
          IF(MODE.EQ.1) THEN
              D=NEW(IB+7:IB+10)*P(17:20)
          ELSE
              D=NEW(IB+1:IB+4)
              V=NEW(IB+21:IB+24)
              VF=MIN(1D0-(1D0-V(1))*(1D0-V(2)),P(28))
              VM=MIN(1D0-(1D0-V(3))*(1D0-V(4)),P(29))
              STATE(7)=MAX(STATE(7),VF)
              STATE(8)=MAX(STATE(8),VM)
          END IF
          DF=1D0-(1D0-D(1))*(1D0-D(2))
          DM=1D0-(1D0-D(3))*(1D0-D(4))
          IF(MODE.EQ.2) THEN
              DF=MIN(DF,P(28))
              DM=MIN(DM,P(29))
          END IF
          STATE(5)=MAX(STATE(5),DF)
          STATE(6)=MAX(STATE(6),DM)
      END DO
      END
