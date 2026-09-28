from pathlib import Path
p=Path('wcm_hashin_umat/src/wcm_core.for');s=p.read_text(encoding='utf-8')
start=s.index('C=======================================================================')
end=s.index('      SUBROUTINE WCM_UPDATE',start)
s='''C=======================================================================
C WCM UMAT core: parameter check, local rotation, stress, Hashin, damage.
C Abaqus orientation provides W axes; +/-theta rotate W to fiber L axes.
C A: fixed degradation. B: energy softening and viscous damage.
C Units follow the input card: MPa, mm, N/mm and step time.
C=======================================================================
'''+s[end:]
s=s.replace('      IU=26\n      IF (MODE.EQ.2) IU=31','      IU=25\n      IF (MODE.EQ.2) IU=30')
s=s.replace('C Fortran逻辑表达式不保证短路；必须先退出A再访问B的P(30)。','C Both finite-difference evaluations use the same OLD history.')
s=s.replace('      INTEGER MODE,NP,IERR,I,IT,IEN,IUP,ISCH,IMOD,IB','      INTEGER MODE,NP,IERR,I,IT,IUP,IB')
s=s.replace('(MODE.EQ.1.AND.NP.NE.28)', '(MODE.EQ.1.AND.NP.NE.25)').replace('(MODE.EQ.2.AND.NP.NE.33)', '(MODE.EQ.2.AND.NP.NE.30)')
s=s.replace('      IT=23\n      IEN=25\n      IUP=26\n      ISCH=27\n      IMOD=28\n      IF (MODE.EQ.2) THEN\n          IT=26\n          IEN=30\n          IUP=31\n          ISCH=32\n          IMOD=33\n      END IF\n      IF (P(ISCH).NE.202609D0.OR.P(IMOD).NE.DBLE(MODE))\n     1    IERR=4', '      IT=23\n      IUP=25\n      IF (MODE.EQ.2) THEN\n          IT=26\n          IUP=30\n      END IF')
s=s.replace('      IF ((P(IEN).NE.0D0.AND.P(IEN).NE.1D0).OR.\n     1    (P(IUP).NE.0D0.AND.P(IUP).NE.1D0)) IERR=8','      IF(P(IUP).NE.0D0.AND.P(IUP).NE.1D0) IERR=8')
s=s.replace('      IEN=25\n      IUP=26\n      IF(MODE.EQ.2) THEN\n          IEN=30\n          IUP=31\n      END IF\n      IF(NINT(P(IEN)).EQ.1) THEN\n          IF(MODE.EQ.1) THEN\n              D=OLD(7:10)*P(17:20)\n          ELSE\n              D=OLD(21:24)\n          END IF\n      END IF', '      IUP=25\n      IF(MODE.EQ.2) IUP=30\n      IF(MODE.EQ.1) THEN\n          D=OLD(7:10)*P(17:20)\n      ELSE\n          D=OLD(21:24)\n      END IF')
s=s.replace('      IF(IADV.EQ.1.AND.NINT(P(IEN)).EQ.1) THEN','      IF(IADV.EQ.1) THEN')
# Group the solver calls by purpose without touching formulas.
s=s.replace('      CALL WCM_CHECK(MODE,P,NP,OLD,DT,LC,IERR)','C--- Check material card and retained history.\n      CALL WCM_CHECK(MODE,P,NP,OLD,DT,LC,IERR)',1)
s=s.replace('      CALL WCM_ELASTIC(P,C0,IERR)','C--- Undamaged local stiffness.\n      CALL WCM_ELASTIC(P,C0,IERR)',1)
s=s.replace('      CALL WCM_EVAL(MODE,P,NP,C0,OLD,E1,DT,LC,IADV,','C--- Rotate both angle families, update stress and damage.\n      CALL WCM_EVAL(MODE,P,NP,C0,OLD,E1,DT,LC,IADV,',1)
s=s.replace('      IF (MODE.NE.2) RETURN\n      IF (NINT(P(IU)).EQ.1', 'C--- B current-feedback option: numerical DDSDDE from fixed OLD.\n      IF (MODE.NE.2) RETURN\n      IF (NINT(P(IU)).EQ.1')
s=s.replace('          CALL WCM_ROTATE(TH,B)','C--- W axes to the local +/-theta fiber axes.\n          CALL WCM_ROTATE(TH,B)')
s=s.replace('          CALL WCM_PLY(MODE,P,C0,OLD(IB+1:IB+40),EL,DT,LC,','C--- Local stress, Hashin criterion and damage update.\n          CALL WCM_PLY(MODE,P,C0,OLD(IB+1:IB+40),EL,DT,LC,')
s=s.replace('          S=S+W*MATMUL(TRANSPOSE(B),SL)','C--- Rotate stress and tangent back; sum the two families.\n          S=S+W*MATMUL(TRANSPOSE(B),SL)')
s=s.replace('      CALL WCM_STIFF(MODE,P,D,C0,C)\n      S=MATMUL(C,E)','C--- Stress and frozen-damage Jacobian for this trial strain.\n      CALL WCM_STIFF(MODE,P,D,C0,C)\n      S=MATMUL(C,E)',1)
s=s.replace('      CALL WCM_HASHIN(EFF,P,FI)','C--- Four Hashin initiation indices.\n      CALL WCM_HASHIN(EFF,P,FI)')
s=s.replace('              CALL WCM_EQ(E,EFF,LC,DEL,SIG)','C--- B: equivalent displacement, fracture energy and viscosity.\n              CALL WCM_EQ(E,EFF,LC,DEL,SIG)')
p.write_text(s,encoding='utf-8')
