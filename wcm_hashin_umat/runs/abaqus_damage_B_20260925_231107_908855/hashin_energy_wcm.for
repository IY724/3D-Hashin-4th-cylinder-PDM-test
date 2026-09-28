C=======================================================================
C WCM 双角度 UMAT：版本 B（能量线性软化）；原单层公式作者 IY。
C 由 tools/build.py 打包；两份dist文件只选择其一提交。
C 本构公式保持原项目，新增旋转层位于 WCM_ROTATE/WCM_EVAL。
C 仅支持NDI=3、NSHR=3、NTENS=6的三维实体；需84个SDV。
C 不支持旧SDV重启动、热耦合、壳或平面应力。
C 适用小材料应变及随动坐标大转动，不是任意有限应变模型。
C 默认保持旧分步时序，不以计算收敛声称真实爆破压力已验证。
C=======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER
      INTEGER KSPT,JSTEP(4),KINC,I,J,IERR,IADV,MODE
      PARAMETER (MODE=2)
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
      DOUBLE PRECISION DDSDDE(NTENS,NTENS),DDSDDT(NTENS)
      DOUBLE PRECISION DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS)
      DOUBLE PRECISION TIME(2),PREDEF(*),DPRED(*),PROPS(NPROPS)
      DOUBLE PRECISION COORDS(3),DROT(3,3),DFGRD0(3,3)
      DOUBLE PRECISION DFGRD1(3,3),SSE,SPD,SCD,RPL,DRPLDT
      DOUBLE PRECISION DTIME,TEMP,DTEMP,PNEWDT,CELENT
      DOUBLE PRECISION NEW(84),S(6),C(6,6),SOLD(6),UOLD,WORK
      RPL=0D0
      DRPLDT=0D0
      DDSDDT=0D0
      DRPLDE=0D0
      SCD=0D0
      IF(NDI.NE.3.OR.NSHR.NE.3.OR.NTENS.NE.6.OR.
     1   NSTATV.NE.84) THEN
          WRITE(7,*) 'WCM_ERROR dimensions: ',NDI,NSHR,NTENS,
     1                NSTATV
          CALL XIT
          RETURN
      END IF
      SOLD=STRESS
      UOLD=SSE
      IADV=1
      IF(KINC.EQ.0.OR.DTIME.EQ.0D0) IADV=0
      CALL WCM_UPDATE(MODE,PROPS,NPROPS,STATEV,STRAN,DSTRAN,
     1                DTIME,CELENT,IADV,NEW,S,C,IERR)
      IF(IERR.NE.0) THEN
C 数值错误码含义见 tools/model.py；不静默改变断裂能或长度。
          WRITE(7,*) 'WCM_ERROR code/mode/element/point: ',
     1               IERR,MODE,NOEL,NPT
          WRITE(7,*) 'WCM_ERROR material/length: ',CMNAME,CELENT
          CALL XIT
          RETURN
      END IF
      STRESS=S
      DDSDDE=C
C NLGEOM中的体积项来自J*sigma的线性化；不对sigma重复旋转。
C 此处仍采用原项目随动小材料应变框架，非有限应变超弹性。
      IF(JSTEP(3).EQ.1) THEN
          DO J=1,3
              DDSDDE(:,J)=DDSDDE(:,J)+STRESS
          END DO
      END IF
      IF(IADV.EQ.1) THEN
          STATEV=NEW
          SSE=0.5D0*DOT_PRODUCT(STRESS,STRAN+DSTRAN)
C SPD为离散应力功减弹性能变化的有符号诊断，含算法/黏性影响。
C 保留旧本构与分步时序时不能将它当作精确断裂耗散。
C 不强行截掉负值来掩盖能量误差；基准测试检查步长敏感性。
          WORK=0.5D0*DOT_PRODUCT(SOLD+STRESS,DSTRAN)
          SPD=SPD+WORK-(SSE-UOLD)
      END IF
      RETURN
      END

C=======================================================================
C 原始单层公式：IY，中国石油大学（华东）。本文件增加 WCM 角度层。
C 不替换 Hashin、Vo 型刚度、等效量、固定损伤值或线性软化公式。
C
C 坐标约定（十分重要）：
C Abaqus 已按 Orientation/AddRot 将全局量变到 WCM 局部坐标 W。
C W 的 1/2/3 为子午/环/厚度向，不是纤维方向！本程序仅做 W->L。
C L+ 与 L- 的纤维方向相对 W1 分别为 +theta、-theta。
C eL=B*eW；sW=B^T*sL；CW=B^T*CL*B。后两式来自功共轭。
C 绝对不能用 B 直接转应力，也不能再做一次全局圆柱坐标变换。
C 工程 Voigt 顺序固定为 11,22,33,12,13,23；剪应变为 gamma。
C 每族40个状态，两族加权后返回完整6x6矩阵（包括拉剪耦合）。
C
C A材料：原22项 + theta,wplus,enabled,update,202609,1 = 28项。
C B材料：原25项 + theta,wplus,capf,capm,enabled,update,
C                     202609,2 = 33项。
C update=0保留原分步算法：旧损伤算应力，新损伤下增量生效。
C update=1本次损伤立即反馈；B用数值算法切线，需UNSYMM。
C 默认0便于隔离角度改动；两种算法的有限步长结果并不相同。
C
C 每族1:12(A)/1:24(B)保留原SDV语义；25:28=FI；
C 29:34=实际返回的单层应力；35:40=单层工程应变。
C 总81=theta,82=wplus,83=版本(A1/B2),84=202609。
C 旧模型的restart不可直接用于此布局。
C=======================================================================
      SUBROUTINE WCM_UPDATE(MODE,P,NP,OLD,E0,DE,DT,LC,IADV,
     1                     NEW,S,C,IERR)
      IMPLICIT NONE
      INTEGER MODE,NP,IADV,IERR,I,J,IU
      DOUBLE PRECISION P(NP),OLD(84),NEW(84),E0(6),DE(6)
      DOUBLE PRECISION DT,LC,S(6),C(6,6),C0(6,6),E1(6)
      DOUBLE PRECISION SP(6),SM(6),EP(6),EM(6),H
      DOUBLE PRECISION TMP(84),CT(6,6)
      IERR=0
      NEW=OLD
      S=0D0
      C=0D0
      CALL WCM_CHECK(MODE,P,NP,OLD,DT,LC,IERR)
      IF (IERR.NE.0) RETURN
      CALL WCM_ELASTIC(P,C0,IERR)
      IF (IERR.NE.0) RETURN
      IU=26
      IF (MODE.EQ.2) IU=31
      E1=E0+DE
      DO I=1,6
          IF (E1(I).NE.E1(I).OR.ABS(E1(I)).GT.1D10) THEN
              IERR=11
              RETURN
          END IF
      END DO
      CALL WCM_EVAL(MODE,P,NP,C0,OLD,E1,DT,LC,IADV,
     1              NEW,S,C,IERR)
      IF (IERR.NE.0) RETURN
C 当前反馈B版的切线：所有扰动都从同一OLD重算，不能污染历史。
C 旧分步算法的C就是严格的冻结损伤导数，无需数值扰动。
C Fortran逻辑表达式不保证短路；必须先退出A再访问B的P(30)。
      IF (MODE.NE.2) RETURN
      IF (NINT(P(IU)).EQ.1.AND.IADV.EQ.1
     1    .AND.NINT(P(30)).EQ.1) THEN
          DO J=1,6
              H=MAX(1D-9,1D-6*ABS(E1(J)))
              EP=E1
              EM=E1
              EP(J)=EP(J)+H
              EM(J)=EM(J)-H
              CALL WCM_EVAL(MODE,P,NP,C0,OLD,EP,DT,LC,IADV,
     1                      TMP,SP,CT,IERR)
              IF (IERR.NE.0) RETURN
              CALL WCM_EVAL(MODE,P,NP,C0,OLD,EM,DT,LC,IADV,
     1                      TMP,SM,CT,IERR)
              IF (IERR.NE.0) RETURN
              DO I=1,6
                  C(I,J)=(SP(I)-SM(I))/(2D0*H)
              END DO
          END DO
      END IF
      END

      SUBROUTINE WCM_CHECK(MODE,P,NP,OLD,DT,LC,IERR)
      IMPLICIT NONE
      INTEGER MODE,NP,IERR,I,IT,IEN,IUP,ISCH,IMOD
      DOUBLE PRECISION P(NP),OLD(84),DT,LC
      IERR=0
      IF ((MODE.EQ.1.AND.NP.NE.28).OR.
     1    (MODE.EQ.2.AND.NP.NE.33).OR.
     2    (MODE.NE.1.AND.MODE.NE.2)) THEN
          IERR=1
          RETURN
      END IF
      DO I=1,NP
          IF(P(I).NE.P(I).OR.ABS(P(I)).GT.1D100) IERR=2
      END DO
      DO I=1,84
          IF(OLD(I).NE.OLD(I).OR.ABS(OLD(I)).GT.1D100) IERR=3
      END DO
      IF (IERR.NE.0) RETURN
      IT=23
      IEN=25
      IUP=26
      ISCH=27
      IMOD=28
      IF (MODE.EQ.2) THEN
          IT=26
          IEN=30
          IUP=31
          ISCH=32
          IMOD=33
      END IF
      IF (P(ISCH).NE.202609D0.OR.P(IMOD).NE.DBLE(MODE))
     1    IERR=4
      IF (OLD(84).NE.0D0) THEN
          IF (OLD(84).NE.202609D0.OR.OLD(83).NE.DBLE(MODE))
     1        IERR=5
          IF(ABS(OLD(81)-P(IT)).GT.1D-12.OR.
     1       ABS(OLD(82)-P(IT+1)).GT.1D-12) IERR=5
      ELSE IF (ANY(OLD(1:80).NE.0D0)) THEN
          IERR=5
      END IF
      IF (MINVAL(P(1:3)).LE.0D0.OR.
     1    MINVAL(P(7:16)).LE.0D0) IERR=6
      IF (P(IT).LT.0D0.OR.P(IT).GT.90D0.OR.
     1    P(IT+1).LT.0D0.OR.P(IT+1).GT.1D0) IERR=7
      IF ((P(IEN).NE.0D0.AND.P(IEN).NE.1D0).OR.
     1    (P(IUP).NE.0D0.AND.P(IUP).NE.1D0)) IERR=8
      IF (MODE.EQ.1) THEN
          IF (MINVAL(P(17:22)).LT.0D0.OR.
     1        MAXVAL(P(17:22)).GT.1D0) IERR=9
      ELSE
          IF (MINVAL(P(17:20)).LE.0D0.OR.
     1        MINVAL(P(22:25)).LT.0D0.OR.
     2        MAXVAL(P(22:23)).GT.1D0.OR.
     3        MINVAL(P(28:29)).LE.0D0.OR.
     4        MAXVAL(P(28:29)).GE.1D0) IERR=9
C PROPS(21)仍是原版保留参数；不擅自把它加入Hashin表达式。
          IF (LC.NE.LC.OR.LC.LE.0D0.OR.LC.GT.1D100) IERR=10
      END IF
      IF (DT.NE.DT.OR.DT.LT.0D0.OR.DT.GT.1D100) IERR=10
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
      DOUBLE PRECISION DF,DM,RF,RM,RS,RMT,RMC
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
C 保留hashin.for的模式上限及组合上限；剪切不额外乘RF。
          DF=MIN(MAX(DF,0D0),P(28))
          DM=MIN(MAX(DM,0D0),P(29))
          RF=1D0-DF
          RM=1D0-DM
          RS=(1D0-P(22)*D(3))*(1D0-P(23)*D(4))
          RS=MIN(MAX(RS,0.05D0),1D0)
      END IF
      C=C0
      DO I=1,3
          DO J=1,3
              C(I,J)=C0(I,J)*RF*RM
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
C 与原版完全相同：FC为平方形式；保持alpha=1的剪切项。
      FI=0D0
      IF(S(1).GE.0D0) THEN
          FI(1)=(S(1)/P(10))**2+(S(4)/P(14))**2
     1         +(S(5)/P(15))**2
      ELSE
          FI(2)=(-S(1)/P(11))**2
      END IF
      Q=S(2)+S(3)
      IF(Q.GE.0D0) THEN
          FI(3)=(Q/P(12))**2+(S(6)**2-S(2)*S(3))/P(16)**2
     1         +(S(4)/P(14))**2+(S(5)/P(15))**2
      ELSE
          FI(4)=((P(13)/(2D0*P(16)))**2-1D0)*Q/P(13)
     1         +Q**2/(2D0*P(16))**2
     2         +(S(6)**2-S(2)*S(3))/P(16)**2
     3         +(S(4)/P(14))**2+(S(5)/P(15))**2
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
          CALL WCM_ROTATE(TH,B)
          EL=MATMUL(B,E)
          CALL WCM_PLY(MODE,P,C0,OLD(IB+1:IB+40),EL,DT,LC,
     1                 IADV,LOCAL,SL,CL,IERR)
          IF(IERR.NE.0) THEN
C 1000/2000标记出错族，余数为实际错误码。
              IERR=IERR+1000*K
              RETURN
          END IF
          NEW(IB+1:IB+40)=LOCAL
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
      INTEGER MODE,IADV,IERR,I,IEN,IUP
      DOUBLE PRECISION P(*),C0(6,6),OLD(40),E(6),DT,LC
      DOUBLE PRECISION NEW(40),S(6),C(6,6),D(4),FI(4),EFF(6)
      DOUBLE PRECISION DEL(4),SIG(4),CAP,DFINAL,KAP,ETA,DNEW
      IERR=0
      NEW=OLD
      D=0D0
      IEN=25
      IUP=26
      IF(MODE.EQ.2) THEN
          IEN=30
          IUP=31
      END IF
      IF(NINT(P(IEN)).EQ.1) THEN
          IF(MODE.EQ.1) THEN
              D=OLD(7:10)*P(17:20)
          ELSE
              D=OLD(21:24)
          END IF
      END IF
      CALL WCM_STIFF(MODE,P,D,C0,C)
      S=MATMUL(C,E)
      EFF=S
C B沿用原程序的有效应力判据(C0*e)，A沿用旧损伤名义试应力。
      IF(MODE.EQ.2) EFF=MATMUL(C0,E)
      CALL WCM_HASHIN(EFF,P,FI)
      IF(IADV.EQ.1.AND.NINT(P(IEN)).EQ.1) THEN
          IF(MODE.EQ.1) THEN
              DO I=1,4
                  IF(FI(I).GE.1D0) NEW(I+6)=1D0
              END DO
              D=NEW(7:10)*P(17:20)
              NEW(1:4)=FI
              NEW(5)=1D0-(1D0-D(1))*(1D0-D(2))
              NEW(6)=1D0-(1D0-D(3))*(1D0-D(4))
              NEW(11)=1D0
          ELSE
              CALL WCM_EQ(E,EFF,LC,DEL,SIG)
              DO I=1,4
                  CAP=P(28)
                  ETA=P(24)
                  IF(I.GE.3) THEN
                      CAP=P(29)
                      ETA=P(25)
                  END IF
                  IF(FI(I).GE.1D0) NEW(I+8)=1D0
                  IF(NEW(I+8).GE.1D0) THEN
C 保持原版在第一次触发的端点取起始量，需检查增量敏感性。
C 不改变此公式；修复旧版非法deltaf时静默不演化的问题。
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
     1                  *(KAP-NEW(I+12))/(KAP*(DFINAL-NEW(I+12)))
C 明确保证不可逆；原损伤上限仍默认纤维0.95、基体0.85。
                      NEW(I)=MIN(MAX(OLD(I),DNEW,0D0),CAP)
                  END IF
                  IF(ETA.GT.0D0.AND.DT.GT.0D0) THEN
                      NEW(I+20)=(ETA*OLD(I+20)+DT*NEW(I))
     1                          /(ETA+DT)
                  ELSE
                      NEW(I+20)=NEW(I)
                  END IF
                  NEW(I+20)=MIN(MAX(NEW(I+20),0D0),CAP)
              END DO
              D=NEW(21:24)
          END IF
          IF(NINT(P(IUP)).EQ.1) THEN
C 可选当前反馈仅改变更新时序，不替换刚度及损伤公式。
              CALL WCM_STIFF(MODE,P,D,C0,C)
              S=MATMUL(C,E)
          END IF
      END IF
      IF(IADV.EQ.1) THEN
          NEW(25:28)=FI
          NEW(29:34)=S
          NEW(35:40)=E
      END IF
      END
