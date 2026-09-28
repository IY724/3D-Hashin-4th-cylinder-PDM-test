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
      INTEGER MODE,NP,NS,K,IB,IT,IUP,IEN,IUSED
      DOUBLE PRECISION P(NP),OLD(84),NEW(84),STATE(NS),W
      DOUBLE PRECISION D(4),V(4),U(4),DF,DM,UF,UM
      STATE=0D0
      STATE(1:4)=-1D100
      IT=23
      IUP=26
      IEN=25
      IUSED=7
      IF(MODE.EQ.1) THEN
          STATE(9:12)=NEW(7:10)
          STATE(13:16)=NEW(47:50)
      ELSE
          IT=26
          IUP=31
          IEN=30
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
          IF(NINT(P(IEN)).EQ.0) U=0D0
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
