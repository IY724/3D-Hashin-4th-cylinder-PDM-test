"""Read-only source review; all new evidence stays in this audit directory."""
from pathlib import Path
import sys, json, hashlib, re, unittest, datetime
from collections import Counter
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT/'wcm_hashin_umat/tools'))
from convert_wcm_inp import blocks
import run_tests as rt

def data(block):
    return [x.strip() for x in block.text.splitlines()[1:]
            if x.strip() and not x.lstrip().startswith('**')]

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def inspect_input():
    source = ROOT/'wcm_hashin_umat/generated/Job_wcm_B.inp'
    text = source.read_bytes().decode('latin-1')
    items = blocks(text)
    report = json.loads(source.with_suffix('.report.json').read_text('utf-8'))
    angles = {r['material'].lower(): r['angle'] for r in report['materials']}
    sets, elements, assignments, rotations, nodes = {}, {}, {}, {}, {}
    duplicates = []
    active = False
    for b in items:
        if b.key == 'part':
            active = b.options['name'].lower() == 'tank-1'
        if not active:
            continue
        if b.key == 'node':
            for row in data(b):
                v = row.split(',')
                nodes[int(v[0])] = np.array(list(map(float, v[1:4])))
        if b.key == 'element':
            for row in data(b):
                v = [int(x) for x in row.split(',') if x.strip()]
                elements[v[0]] = (b.options['type'], v[1:])
        if b.key == 'elset':
            name = b.options['elset'].lower()
            entries = []
            for row in data(b):
                tokens = [x.strip() for x in row.split(',') if x.strip()]
                if 'generate' in b.options:
                    vals = list(map(int, tokens))
                    entries += list(range(vals[0], vals[1]+1, vals[2]))
                else:
                    entries += [int(x) if x.isdigit() else x.lower() for x in tokens]
            sets.setdefault(name, []).extend(entries)
        if b.key == 'distribution' and b.options['name'].lower() == 'wcm_tank1_addrot':
            for row in data(b):
                vals = row.split(',')
                if vals[0].strip():
                    rotations[int(vals[0])] = float(vals[1])
        if b.key == 'end part':
            active = False
    def resolve(name):
        out = set()
        for v in sets[name]:
            if isinstance(v, int): out.add(v)
            else: out.update(resolve(v))
        return out
    active = False
    for b in items:
        if b.key == 'part': active = b.options['name'].lower() == 'tank-1'
        if active and b.key == 'solid section':
            for eid in resolve(b.options['elset'].lower()):
                if eid in assignments: duplicates.append(eid)
                assignments[eid] = (angles[b.options['material'].lower()], b.options['orientation'])
        if b.key == 'end part': active = False
    layers = []
    for name in sorted(sets):
        if not re.fullmatch('layer[0-9]+', name): continue
        ids = resolve(name)
        row = {'layer': name, 'elements': len(ids)}
        centroids = {e:np.mean([nodes[n] for n in elements[e][1]],axis=0) for e in ids}
        ys = [v[1] for v in centroids.values()]
        midpoint, span = (min(ys)+max(ys))/2., max(ys)-min(ys)
        middle = {e for e,v in centroids.items() if abs(v[1]-midpoint)<0.05*span}
        row['geometric_midbody_angles'] = sorted({assignments[e][0] for e in middle})
        center = resolve(name+'_center') if name+'_center' in sets else {e for e in ids if abs(rotations[e]-90.) < 1e-8}
        ends = resolve(name+'_endcaps') if name+'_endcaps' in sets else ids-center
        row['region_basis'] = 'named_sets' if name+'_center' in sets else 'AddRot=90 proxy; not independent geometry validation'
        for label, selection in [('all', ids), ('center', center), ('endcaps', ends)]:
            vals = sorted({assignments[e][0] for e in selection})
            rot = [rotations[e] for e in selection]
            row[label] = {'elements':len(selection), 'angles':vals,
                          'addrot_range':[min(rot),max(rot)] if rot else []}
        layers.append(row)
    evidence = {'input_sha256':digest(source),
                'input_matches_conversion_report':digest(source)==report['output_sha256'],
                'elements': len(elements), 'element_types':dict(Counter(x[0] for x in elements.values())),
                'assigned_elements':len(assignments), 'duplicates':len(duplicates),
                'missing_sections':len(set(elements)-set(assignments)),
                'missing_addrot':len(set(elements)-set(rotations)), 'layers':layers}
    return evidence

def verify_saved_runs():
    result = []
    for path in (ROOT/'wcm_hashin_umat/validation').glob('abaqus_*.json'):
        for entry in json.loads(path.read_text('utf-8')):
            run = Path(entry['directory'])
            name = 'hashin_constant_wcm.for' if entry['variant']=='A' else 'hashin_energy_wcm.for'
            dist = ROOT/'wcm_hashin_umat/dist'/name
            template = (ROOT/'wcm_hashin_umat/src/umat_entry.for.in').read_text('utf-8')
            mode = '1' if entry['variant']=='A' else '2'
            version = 'A（直接刚度折减）' if mode=='1' else 'B（能量线性软化）'
            expected = template.replace('@MODE@',mode).replace('@VERSION@',version)+'\n'+(ROOT/'wcm_hashin_umat/src/wcm_core.for').read_text('utf-8')
            content = (run/'launcher.txt').read_text(errors='replace')
            sta = (run/'verify.sta').read_text(errors='replace') if (run/'verify.sta').exists() else ''
            dat = (run/'verify.dat').read_text(errors='replace')
            result.append({'suite':entry['suite'],'variant':entry['variant'],
                'current_dist_matches_tested_source':digest(dist)==entry['umat_sha256']==digest(run/name),
                'dist_matches_current_source':dist.read_bytes()==expected.encode('utf-8'),
                'launcher_completed':'Abaqus JOB verify COMPLETED' in content,
                'analysis_completed':'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' in sta,
                'datacheck_completed':'ANALYSIS DATACHECK COMPLETE' in dat,
                'warnings':dict(Counter(re.findall(r'\*\*\*WARNING:([^\n]+)',dat))),
                'datacheck_input_matches_current':digest(run/'verify.inp')==digest(ROOT/f'wcm_hashin_umat/generated/Job_wcm_{entry["variant"]}.inp') if entry['suite']=='datacheck' else None})
    return result

def offaxis_history(core):
    results=[]
    for variant in ('A','B'):
        for angle in (15.,20.,25.,30.,35.,40.,60.):
            p=rt.make_props(variant,rt.BASE[variant],angle)
            state=np.zeros(84)
            legacy=[np.zeros(24),np.zeros(24)]
            olde=np.zeros(6)
            worst=0.
            for amplitude in np.r_[np.linspace(0,1,101),np.linspace(1,0,31),np.linspace(0,1.1,31)]:
                e=amplitude*np.array([.025,.002,-.0001,.018,.002,.001])
                state,s,c,_=core.update(variant,p,state,olde,e-olde,length=.01)
                ref=np.zeros(6)
                for k,sign in enumerate((1,-1)):
                    local=rt.rotate_strain(e,sign*angle)
                    local0=rt.rotate_strain(olde,sign*angle)
                    legacy[k],ls,lc=core.original(variant,legacy[k],local0,local-local0,length=.01)
                    ref+=.5*rt.rotate_stress_back(ls,sign*angle)
                    count=12 if variant=='A' else 24
                    np.testing.assert_allclose(state[k*40:k*40+count],legacy[k][:count],rtol=1e-9,atol=1e-9)
                worst=max(worst,float(np.max(abs(s-ref))/max(1.,np.max(abs(ref)))))
                np.testing.assert_allclose(s,ref,rtol=1e-9,atol=1e-9)
                olde=e
            results.append({'variant':variant,'angle':angle,'max_relative_stress_error':worst,
                            'family_histories_match_legacy':True,'final_plus':state[:10].tolist(),
                            'final_minus':state[40:50].tolist()})
    return results

if __name__ == '__main__':
    evidence = {'mapping':inspect_input(), 'saved_runs':verify_saved_runs()}
    rt.MaterialTests.core = rt.Core(HERE/('core_recheck_'+datetime.datetime.now().strftime('%H%M%S_%f')))
    with (HERE/'core_recheck.txt').open('w',encoding='utf-8') as log:
        tests = unittest.defaultTestLoader.loadTestsFromTestCase(rt.MaterialTests)
        result = unittest.TextTestRunner(stream=log,verbosity=2).run(tests)
    evidence['fresh_core_tests'] = {'run':result.testsRun,'passed':result.wasSuccessful(),
                                    'metrics':rt.MaterialTests.metrics}
    evidence['fresh_offaxis_history'] = offaxis_history(rt.MaterialTests.core)
    (HERE/'audit_evidence.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'core_tests':evidence['fresh_core_tests'],'offaxis_cases':len(evidence['fresh_offaxis_history'])},ensure_ascii=False))
