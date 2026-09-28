from pathlib import Path
import subprocess,shutil,json
r=Path.cwd();records=json.loads((r/'wcm_hashin_umat/validation/abaqus_completion.json').read_text());out=[]
for rec in records:
 v=rec['variant'];d='A_stiffness_reduction' if v=='A' else 'B_energy_softening';folder=r/'test'/d
 cfg=json.loads((folder/'failure_thresholds.json').read_text());cfg['pressure_at_step_end_mpa']=1.
 p=Path(rec['directory'])/'threshold_test_config.json';p.write_text(json.dumps(cfg),encoding='utf-8')
 result=Path(rec['directory'])/'threshold_test_result.json'
 cmd=[shutil.which('abaqus'),'python',str(folder/'extract_failure.py'),str(Path(rec['directory'])/'verify.odb'),'--config',str(p),'--output',str(result)]
 x=subprocess.run(cmd,capture_output=True);(Path(rec['directory'])/'threshold_test.log').write_bytes(x.stdout+x.stderr)
 if x.returncode or not result.exists():raise RuntimeError((x.stdout+x.stderr).decode(errors='replace'))
 m=json.loads(result.read_text(encoding='utf-8'));assert m['events']['fiber_completion'] is not None and m['events']['matrix_completion'] is not None
 out.append({'variant':v,'passed':True,'scope':'Extractor on completion test ODB; 0..1 pressure mapping is synthetic, not a vessel pressure','result':str(result)})
(r/'audit_20260927/damage_one/postprocess_tests.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out))
