from pathlib import Path
import shutil,subprocess,sys,json
root=Path.cwd();b=root/'test/B_energy_softening';stage=root/'audit_20260926/sdv_compact/B_regenerated'
for rel in ('source','config','scripts'):shutil.copytree(b/rel,stage/rel,dirs_exist_ok=True)
subprocess.run([sys.executable,str(stage/'scripts/prepare_case.py')],check=True)
sys.path.insert(0,str(root/'wcm_hashin_umat/tools'))
from convert_wcm_inp import blocks
skip={'depvar','output','element output','node output','energy output'}
def protected(t):return [(x.key,'\n'.join(l for l in x.text.splitlines() if not l.startswith('**')).strip()) for x in blocks(t) if x.key not in skip]
assert protected((stage/'Job2_PDM_B.inp').read_text())==protected((b/'Job2_PDM_B.inp').read_text())
for name in ('Job2_PDM_B.inp','hashin_energy_wcm.for'):shutil.copy2(stage/name,b/name)
shutil.copytree(stage/'results',b/'results',dirs_exist_ok=True)
subprocess.run([sys.executable,str(b/'scripts/prepare_case.py'),'--check'],check=True)
