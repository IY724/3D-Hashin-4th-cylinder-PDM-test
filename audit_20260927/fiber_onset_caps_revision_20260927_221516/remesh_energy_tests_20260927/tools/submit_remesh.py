# -*- coding: utf-8 -*-
"""新加密网格：替换B材料卡及输出后直接提交，不另做datacheck。"""
from pathlib import Path
import datetime
import hashlib
import json
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
CAND = PROJECT / 'wcm_hashin_umat_B_candidate'
sys.path.insert(0, str(CAND / 'tools'))
from convert_wcm_inp import blocks, data_end, transform


def patch(text, edits):
    edits.sort()
    if any(a[1] > b[0] for a, b in zip(edits, edits[1:])):
        raise ValueError('转换区间重叠')
    for start, end, replacement in reversed(edits):
        text = text[:start] + replacement + text[end:]
    return text


def main():
    source = PROJECT / 'geometry_zhuning_joint_rebuild' / 'Job-remesh.inp'
    raw = source.read_bytes()
    text = raw.decode('latin-1').replace('\r\n', '\n')
    config = json.loads((CAND / 'config/project.json').read_text(encoding='utf-8'))
    for material in config['materials'].values():
        material['B_damage'][:4] = [133.0, 40.0, 0.6, 2.1]
        material['B_damage'][5:7] = [0.9, 0.5]
    config['update'] = 0
    config['fiber_cap'] = config['matrix_cap'] = 1.0

    # 合并新INP中两处UVARM场输出，保留其余模型定义。
    edits = []
    first_field = True
    for block in blocks(text):
        if block.key == 'static':
            edits.append((block.start, data_end(block), '*Static\n0.005, 1., 1e-20, 0.02\n'))
        elif block.key == 'element output':
            replacement = '*Element Output, directions=YES\nLE, PE, PEEQ, PEMAG, S, SDV\n' if first_field else ''
            first_field = False
            edits.append((block.start, data_end(block), replacement))
        elif block.key == 'node output':
            edits.append((block.start, data_end(block), '*Node Output\nU, RF\n'))
        elif block.key == 'output' and 'history' in block.options:
            edits.append((block.start, data_end(block),
                          '*Output, history, frequency=1\n*Energy Output\nALLAE, ALLIE, ALLSE, ALLPD, ALLSD, ALLVD, ALLWK\n'))
    text = patch(text, edits)
    converted, records = transform(text, 'B', config)
    targets = {record['material'].lower() for record in records}

    # 新模型减缩积分复材截面接入UMAT时启用Enhanced沙漏控制。
    edits = []
    count = 0
    for block in blocks(converted):
        if block.key == 'section controls':
            raise ValueError('源模型已有截面控制，不能重复添加')
        if block.key == 'solid section' and block.options.get('material', '').lower() in targets:
            heading = block.text.splitlines()[0]
            if 'controls' in block.options:
                raise ValueError('复材截面已有controls')
            edits.append((block.start, block.start + len(heading), heading + ', controls=HG_UMAT'))
            count += 1
    first_part = next(b for b in blocks(converted) if b.key == 'part')
    edits.append((first_part.start, first_part.start,
                  '*Section Controls, name=HG_UMAT, hourglass=ENHANCED\n'))
    converted = patch(converted, edits)

    name = 'B_remesh_M21'
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    run = ROOT / 'runs' / (name + '_' + stamp)
    run.mkdir(parents=True, exist_ok=False)
    (run / 'scratch').mkdir()
    inp = converted.replace('\n', '\r\n').encode('latin-1')
    umat = (CAND / 'dist/hashin_energy_wcm.for').read_bytes()
    (run / (name + '.inp')).write_bytes(inp)
    (run / 'hashin_energy_wcm.for').write_bytes(umat)
    command = [r'E:\ABAQUS2025\Commands\abaqus.BAT', 'job=' + name,
               'input=' + name + '.inp', 'user=hashin_energy_wcm.for',
               'cpus=12', 'scratch=' + str(run / 'scratch'), 'interactive']
    manifest = {
        'source': str(source), 'source_sha256': hashlib.sha256(raw).hexdigest(),
        'inp_sha256': hashlib.sha256(inp).hexdigest(),
        'umat_sha256': hashlib.sha256(umat).hexdigest(),
        'fracture_energies_N_per_mm': [133.0, 40.0, 0.6, 2.1],
        'energy_source': '张华伟等2021，碳纤维复合材料层合板低速冲击影响因素，表1，T700GC/M21',
        'doi': '10.3969/j.issn.1007-2012.2021.12.028',
        'source_url': 'https://qikan.cmes.org/sxgcxb/EN/PDF/10.3969/j.issn.1007-2012.2021.12.028',
        'limitation': '文献参数敏感性试算，非本材料实验标定；原文非能量软化验证，不保证收敛或设计爆压。',
        'smt_smc': [0.9, 0.5], 'viscosity': [0.0001, 0.005],
        'static': [0.005, 1.0, 1e-20, 0.02], 'pressure_mpa': 250.0,
        'update': 0, 'caps': [1.0, 1.0], 'material_bins': len(records),
        'controlled_sections': count,
        'max_elastic_error': max(r['relative_elastic_error'] for r in records),
        'command': command, 'cwd': str(run), 'status': 'SUBMITTING'
    }
    def save():
        output = json.dumps(manifest, ensure_ascii=False, indent=2) + '\n'
        (run / 'run_manifest.json').write_text(output, encoding='utf-8')
        (ROOT / 'latest_run.json').write_text(output, encoding='utf-8')
    save()
    print(json.dumps({'run_directory': str(run), 'job': name,
                      'materials': len(records), 'sections': count}, ensure_ascii=False), flush=True)
    with (run / 'launcher.txt').open('wb') as log:
        process = subprocess.Popen(command, cwd=run, stdout=log, stderr=subprocess.STDOUT)
        manifest['status'] = 'RUNNING'
        manifest['launcher_pid'] = process.pid
        save()
        returncode = process.wait()
    launcher = (run / 'launcher.txt').read_text(errors='replace')
    completed = returncode == 0 and ('Abaqus JOB ' + name + ' COMPLETED') in launcher
    manifest['returncode'] = returncode
    manifest['status'] = 'COMPLETED' if completed else 'STOPPED_CHECK_LOG'
    save()
    print(json.dumps({'status': manifest['status'], 'run_directory': str(run)}, ensure_ascii=False), flush=True)
    return 0 if completed else 1


if __name__ == '__main__':
    sys.exit(main())
