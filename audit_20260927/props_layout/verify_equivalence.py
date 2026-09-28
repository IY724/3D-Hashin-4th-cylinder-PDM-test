from pathlib import Path
import ctypes as ct, subprocess, shutil, json, datetime, sys
import numpy as np
root=Path.cwd(); sys.path.insert(0,str(root/'wcm_hashin_umat/tools'))
from run_tests import Core, BASE, ptr
from model import make_props
run=root/'wcm_hashin_umat/runs'/('props_equivalence_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
core=Core(run)
old_source=root/'audit_20260927/props_layout/backup/wcm_core.for'
subprocess.run([shutil.which('gfortran'),'-shared','-O1','-frecursive','-fcheck=all','-ffixed-line-length-72',str(old_source),'-o',str(run/'prior.dll')],check=True)
prior=ct.CDLL(str(run/'prior.dll'));rng=np.random.default_rng(20260927);comparisons=0; max_stress=0.; max_tangent=0.
for variant in ('A','B'):
 for update in (0,1):
  for angle,weight in ((0.,1.),(15.,.5),(37.,.37),(90.,1.)):
   props=make_props(variant,BASE[variant],angle,weight=weight,update=update)
   old_props=np.array(list(BASE[variant])+[angle,weight]+([] if variant=='A' else [1.,1.])+[1.,update,202609.,1. if variant=='A' else 2.],dtype=np.float64)
   old_state=np.zeros(84); new_state=np.zeros(84); e0=np.zeros(6)
   for value in np.r_[np.linspace(0,1.4,20),np.linspace(1.4,-.4,10)]:
    strain=np.array([.025,.014,-.001,.009,.003,.002])*value
    strain+=rng.normal(0,1e-6,6)
    de=np.array(strain-e0,dtype=np.float64)
    mode=ct.c_int(1 if variant=='A' else 2);n=ct.c_int(len(old_props));advance=ct.c_int(1);err=ct.c_int(0)
    dt=ct.c_double(.01);lc=ct.c_double(.01);old_out=np.zeros(84);s_old=np.zeros(6);c_old=np.zeros((6,6),order='F')
    prior.wcm_update_(ct.byref(mode),ptr(old_props),ct.byref(n),ptr(old_state),ptr(e0),ptr(de),ct.byref(dt),ct.byref(lc),ct.byref(advance),ptr(old_out),ptr(s_old),ptr(c_old),ct.byref(err))
    new_out,s_new,c_new,new_err=core.update(variant,props,new_state,e0,de,dt=.01,length=.01,allow_error=True)
    assert err.value==new_err,(variant,angle,err.value,new_err)
    if err.value==0:
     np.testing.assert_allclose(old_out,new_out,atol=1e-12,rtol=1e-12)
     np.testing.assert_allclose(s_old,s_new,atol=1e-12,rtol=1e-12)
     np.testing.assert_allclose(c_old,c_new,atol=1e-9,rtol=1e-12)
     max_stress=max(max_stress,float(np.max(abs(s_old-s_new))))
     max_tangent=max(max_tangent,float(np.max(abs(c_old-c_new))))
     old_state=old_out;new_state=new_out
    comparisons+=1;e0=strain
result={'passed':True,'comparisons':comparisons,'max_stress_difference':max_stress,'max_tangent_difference':max_tangent,'scope':'Compiled prior vs current WCM core using equivalent active material cards; A/B, both update modes, four angles and cyclic load paths','directory':str(run)}
(root/'audit_20260927/props_layout/equivalence.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
