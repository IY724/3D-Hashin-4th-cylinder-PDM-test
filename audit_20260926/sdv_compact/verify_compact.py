from pathlib import Path
import sys,ctypes as ct,subprocess,shutil,json,datetime
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'wcm_hashin_umat/tools'))
from run_tests import Core,BASE,ptr
from model import make_props,state_count,output_indices
run=ROOT/'wcm_hashin_umat/runs'/('compact_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
core=Core(run)
s=(ROOT/'wcm_hashin_umat/src/umat_entry.for.in').read_text(encoding='utf-8')
(run/'pack.for').write_text(s[s.index('C Persistent history'):],encoding='utf-8')
subprocess.run([shutil.which('gfortran'),'-shared','-fcheck=all','-ffixed-line-length-72',str(run/'pack.for'),'-o',str(run/'pack.dll')],check=True)
lib=ct.CDLL(str(run/'pack.dll')); checks=0; maxerr=0.
for variant in ('A','B'):
 for update in (0,1):
  for weight in (0.,.35,1.):
   for angle in (0.,37.,90.):
    p=make_props(variant,BASE[variant],angle,weight=weight,update=update)
    mode=ct.c_int(1 if variant=='A' else 2);ns=ct.c_int(state_count(variant));np_=ct.c_int(len(p))
    old=np.zeros(84); packed=np.zeros(ns.value); e0=np.zeros(6)
    # Includes initiation, evolution, unloading, reversal, mixed loading.
    for amplitude in [0,.2,.5,1,1.5,2,1,.2,-.2,-.5,-1,-1.5,0,.5]:
     e=np.array([.025,.018,-.003,.012,.004,.003])*amplitude
     restored=np.zeros(84)
     lib.wcm_unpack_(ct.byref(mode),ptr(packed),ct.byref(ns),ptr(restored))
     ref,sref,cref,_=core.update(variant,p,old,e0,e-e0,length=.01)
     new,stress,c,_=core.update(variant,p,restored,e0,e-e0,length=.01)
     np.testing.assert_allclose(stress,sref,rtol=1e-12,atol=1e-10)
     np.testing.assert_allclose(c,cref,rtol=1e-12,atol=1e-8)
     lib.wcm_pack_(ct.byref(mode),ptr(p),ct.byref(np_),ptr(new),ptr(packed),ct.byref(ns))
     active=[i for i,w in ((0,weight),(40,1-weight)) if w>0]
     np.testing.assert_allclose(packed[:4],np.max([new[i+24:i+28] for i in active],axis=0))
     def combined(d,cap=None):
      vals=1-(1-d[::2])*(1-d[1::2])
      return vals if cap is None else np.minimum(vals,cap)
     ds=[combined(new[i+6:i+10]*p[16:20]) if variant=='A' else combined(new[i:i+4],p[27:29]) for i in active]
     np.testing.assert_allclose(packed[4:6],np.max(ds,axis=0))
     if variant=='B':
      np.testing.assert_allclose(packed[6:8],np.max([combined(new[i+20:i+24],p[27:29]) for i in active],axis=0))
     maxerr=max(maxerr,float(np.max(abs(stress-sref))));checks+=1;old=ref;e0=e
result={'passed':True,'increment_comparisons':checks,'max_stress_difference':maxerr,'scope':'Compiled pack/unpack versus retained 84-state core; both timing modes, angles, zero-weight families, unload/reversal','directory':str(run)}
(ROOT/'audit_20260926/sdv_compact/regression.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
