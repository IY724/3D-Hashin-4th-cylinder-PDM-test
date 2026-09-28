from pathlib import Path
import json,shutil,sys,subprocess,hashlib
r=Path.cwd();audit=r/'audit_20260927/damage_one';b=r/'test/B_energy_softening'
for p in [r/'wcm_hashin_umat/config/project.json',b/'config/project.json']:
 m=json.loads(p.read_text(encoding='utf-8'));m['fiber_cap']=1.;m['matrix_cap']=1.;m['说明']='B非黏性/黏性损伤演化到1；cap仅限制法向刚度折减，WCM_STIFF保留1e-6残余刚度因子。';p.write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
for n in ['model.py','convert_wcm_inp.py']:shutil.copy2(r/'wcm_hashin_umat/tools'/n,b/'scripts'/n)
for n in ['wcm_core.for','umat_entry.for.in']:shutil.copy2(r/'wcm_hashin_umat/src'/n,b/'source'/n)
p=b/'source/origin_manifest.json';m=json.loads(p.read_text(encoding='utf-8'))
for e in m['files']:e['sha256']=hashlib.sha256((b/e['snapshot']).read_bytes()).hexdigest()
p.write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
stage=audit/'B_regenerated'
for rel in ['source','config','scripts']:shutil.copytree(b/rel,stage/rel,dirs_exist_ok=True)
subprocess.run([sys.executable,str(stage/'scripts/prepare_case.py')]+(['--check'] if (stage/'Job2_PDM_B.inp').exists() else []),check=True,stdout=subprocess.DEVNULL)
for n in ['Job2_PDM_B.inp','hashin_energy_wcm.for']:shutil.copy2(stage/n,b/n)
shutil.copytree(stage/'results',b/'results',dirs_exist_ok=True)
sys.path.insert(0,str(r/'wcm_hashin_umat/tools'))
from model import state_count,sdv_labels,output_indices
from convert_wcm_inp import blocks,data_end,select_sdv_outputs
p=r/'test/A_stiffness_reduction/Job2_PDM_A.inp';s=p.read_text(encoding='utf-8')
for x in reversed(blocks(s)):
 if x.key=='depvar':s=s[:x.start]+f'*Depvar\n{state_count("A")},\n'+s[data_end(x):]
s=select_sdv_outputs(s,'A');p.write_text(s,encoding='utf-8')
shutil.copy2(r/'wcm_hashin_umat/dist/hashin_constant_wcm.for',p.parent/'hashin_constant_wcm.for')
records={}
for v,d,n in [('A','A_stiffness_reduction','hashin_constant_wcm.for'),('B','B_energy_softening','hashin_energy_wcm.for')]:
 folder=r/'test'/d;current=(folder/f'Job2_PDM_{v}.inp').read_text();old=(audit/'backup/test'/d/f'Job2_PDM_{v}.inp').read_text()
 skip={'depvar','element output','user material'}
 def protected(t):return [(x.key,'\n'.join(l for l in x.text.splitlines() if not l.startswith('**')).strip()) for x in blocks(t) if x.key not in skip]
 assert protected(old)==protected(current)
 olds=[x.numbers() for x in blocks(old) if x.key=='user material'];news=[x.numbers() for x in blocks(current) if x.key=='user material']
 for op,np_ in zip(olds,news):
  if v=='B':op[27:29]=[1.,1.]
  assert op==np_
 (folder/'SDV_MAP.csv').write_text('index,name,meaning,output\n'+'\n'.join(f'{i},SDV{i},{s},{int(i in output_indices(v))}' for i,s in enumerate(sdv_labels(v),1))+'\n',encoding='utf-8-sig')
 records[v]={'depvar':state_count(v),'output_count':len(output_indices(v)),'protected_input_equal':True,'material_count':len(news),'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [folder/f'Job2_PDM_{v}.inp',folder/n]},'changes':'SDV allocation/output; B stiffness caps .95/.85 -> 1/1'}
 config={'variant':v,'fiber_threshold':1.,'matrix_threshold':1.,'tolerance':1e-6,'vessel_failure_rule':'either','step_name':None,'pressure_at_step_start_mpa':0.,'pressure_at_step_end_mpa':250.,'step_duration':1.,'pressure_assumption':'Linear ramp; update for any changed amplitude/load/step.'}
 (folder/'failure_thresholds.json').write_text(json.dumps(config,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
(audit/'delivery.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
subprocess.run([sys.executable,str(b/'scripts/prepare_case.py'),'--check'],check=True,stdout=subprocess.DEVNULL)
print(json.dumps(records))
