from pathlib import Path
import sys,ctypes as ct,json,datetime,numpy as np
r=Path.cwd();sys.path.insert(0,str(r/'wcm_hashin_umat/tools'))
from run_tests import Core,BASE,ptr
from model import make_props,elastic
run=r/'wcm_hashin_umat/runs'/('completion_probe_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'));core=Core(run)
metrics={};c0=elastic(BASE['B'][:9]);L=np.linalg.cholesky(c0);Li=np.linalg.inv(L)
minimum=1.
for variant in ('A','B'):
 for cap in (1.,.95):
  base=BASE[variant].copy()
  if variant=='A':base[16:22]=[1.]*6
  p=make_props(variant,base,0,capf=cap,capm=cap)
  for d in [np.ones(4),np.array([1.,0,1.,0]),np.zeros(4)]:
   c=np.zeros((6,6),order='F');m=ct.c_int(1 if variant=='A' else 2)
   core.library.wcm_stiff_(ct.byref(m),ptr(p),ptr(d),ptr(c0.copy(order='F')),ptr(c))
   eigen=np.linalg.eigvalsh(Li@c@Li.T);minimum=min(minimum,float(eigen.min()))
   assert eigen.min()>=.999999e-6 and np.isfinite(c).all()
metrics['minimum_relative_elastic_eigenvalue']=minimum
# Independent analytic law for monotonic uniaxial FT, with finite viscosity.
base=BASE['B'].copy();base[16]=4.;base[-2:]=[.02,.02]
p=make_props('B',base,0);state=np.zeros(84);e0=np.zeros(6);dt=.01;lc=.2
vexpected=0.;observed_lag=False;completed=False
for value in np.linspace(0,.03,101):
 direction=np.linalg.solve(c0,[1.,0,0,0,0,0]);direction/=direction[0];e=direction*value
 old=state.copy();state,stress,c,_=core.update('B',p,state,e0,e-e0,dt=dt,length=lc)
 if state[8]:
  dfinal=2*base[16]/state[16];k=state[4];d0=state[12]
  dref=max(old[0],min(1.,dfinal*(k-d0)/(k*(dfinal-d0))))
  np.testing.assert_allclose(state[0],dref,atol=1e-14)
 vexpected=(.02*vexpected+dt*state[0])/(.02+dt)
 np.testing.assert_allclose(state[20],vexpected,atol=1e-14)
 if state[0]==1. and not completed:
  observed_lag=state[20]<1. and old[20]<state[20]
  metrics['first_completion']={'raw_damage':float(state[0]),'updated_viscous_damage':float(state[20]),'used_viscous_damage':float(old[20]),'strain':float(value)};completed=True
 e0=e
assert completed and observed_lag
# Holding full damage: approach 1 without healing; residual remains positive.
for _ in range(200):state,stress,c,_=core.update('B',p,state,e0,np.zeros(6),dt=dt,length=lc)
assert state[0]==1. and np.linalg.eigvalsh(c).min()>0
metrics['passed']=True;metrics['scope']='Compiled full-damage residual stiffness and independent energy/viscosity recurrence; material-point verification only'
(r/'audit_20260927/damage_one/completion_probe.json').write_text(json.dumps(metrics,indent=2),encoding='utf-8');print(json.dumps(metrics))
