from pathlib import Path
import sys,shutil,json,hashlib
root=Path.cwd();sys.path.insert(0,str(root/'wcm_hashin_umat/tools'))
from build import build
from model import state_count,output_indices,sdv_labels
from convert_wcm_inp import blocks,data_end,select_sdv_outputs
build();records={}
for v,d,name in [('A','A_stiffness_reduction','hashin_constant_wcm.for'),('B','B_energy_softening','hashin_energy_wcm.for')]:
 folder=root/'test'/d;p=folder/f'Job2_PDM_{v}.inp';s=p.read_text(encoding='utf-8');edits=[]
 for b in blocks(s):
  if b.key=='depvar':edits.append((b.start,data_end(b),f'*Depvar\n{state_count(v)},\n'))
  if b.key=='element output':edits.append((b.start,data_end(b),'*Element Output, directions=YES\nS, SDV\n'))
  if b.key=='node output':edits.append((b.start,data_end(b),'*Node Output\nU\n'))
  if b.key=='energy output' or (b.key=='output' and 'history' in b.text.splitlines()[0].lower()):edits.append((b.start,data_end(b),''))
 for a,b,t in sorted(edits,reverse=True):s=s[:a]+t+s[b:]
 s=select_sdv_outputs(s,v);p.write_text(s,encoding='utf-8');shutil.copy2(root/'wcm_hashin_umat/dist'/name,folder/name)
 # Confirm everything except SDV allocation and output requests unchanged.
 old=(root/'audit_20260926/sdv_compact/backup/test'/d/p.name).read_text(encoding='utf-8')
 skip={'depvar','output','element output','node output','energy output'}
 def protected(text):return [(b.key,'\n'.join(l for l in b.text.splitlines() if not l.startswith('**')).strip()) for b in blocks(text) if b.key not in skip]
 assert protected(old)==protected(s)
 assert all(b.text.splitlines()[1]==f'{state_count(v)},' for b in blocks(s) if b.key=='depvar')
 (folder/'SDV_MAP.csv').write_text('index,name,meaning,output\n'+'\n'.join(f'{i},SDV{i},{label},{int(i in output_indices(v))}' for i,label in enumerate(sdv_labels(v),1))+'\n',encoding='utf-8-sig')
 records[v]={'depvar':state_count(v),'output_count':len(output_indices(v)),'protected_input_equal':True,'material_count':sum(b.key=='user material' for b in blocks(s)),'files':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (p,folder/name)}}
(root/'audit_20260926/sdv_compact/delivery.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
# Update package generator to the requested S/U/selected SDV output contract.
p=root/'test/scripts/templates/prepare_case.py';s=p.read_text(encoding='utf-8').replace('LE, PE, PEEQ, PEMAG, S, SDV','S, SDV').replace('U, RF','U')
s=s.replace("             (ends[0].start, ends[0].start,\n              '*Energy Output\\nALLAE, ALLIE, ALLSE, ALLPD, ALLSD, ALLWK\\n')]", "             *[(b.start, data_end(b), '') for b in items\n               if b.key == 'output' and 'history' in b.text.splitlines()[0].lower()]]")
s=s.replace('UVARM -> SDV、RF/能量输出','UVARM -> S/U/指定SDV输出');p.write_text(s,encoding='utf-8')
b=root/'test/B_energy_softening';backup=root/'audit_20260926/sdv_compact/backup/B_support';backup.mkdir(parents=True,exist_ok=True)
for rel in ['scripts','source/umat_entry.for.in','source/wcm_core.for','source/origin_manifest.json','results']:
 source=b/rel;target=backup/rel
 if source.is_dir():shutil.copytree(source,target,dirs_exist_ok=True)
 else:target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
for n in ['model.py','convert_wcm_inp.py']:shutil.copy2(root/'wcm_hashin_umat/tools'/n,b/'scripts'/n)
shutil.copy2(p,b/'scripts/prepare_case.py')
for n in ['umat_entry.for.in','wcm_core.for']:shutil.copy2(root/'wcm_hashin_umat/src'/n,b/'source'/n)
p=b/'source/origin_manifest.json';m=json.loads(p.read_text(encoding='utf-8'))
for entry in m['files']:
 entry['sha256']=hashlib.sha256((b/entry['snapshot']).read_bytes()).hexdigest()
p.write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(records))
