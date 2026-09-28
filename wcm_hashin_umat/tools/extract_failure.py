# -*- coding: utf-8 -*-
"""Run with Abaqus Python: extract_failure.py result.odb [--config file]."""
import argparse,json,math
from pathlib import Path

def scan(odb, cfg):
    variant=cfg['variant'];required=8 if variant=='A' else 10
    for k in ('fiber_threshold','matrix_threshold'):
        if not 0<float(cfg[k])<=1:raise ValueError(k+' must be in (0,1]')
    tol=float(cfg['tolerance'])
    if not 0<=tol<min(float(cfg['fiber_threshold']),float(cfg['matrix_threshold'])):raise ValueError('Invalid tolerance')
    if cfg['vessel_failure_rule'] not in ('either','fiber','matrix'):raise ValueError('Invalid failure rule')
    names=list(odb.steps.keys());name=cfg.get('step_name')
    if name is None:
        if len(names)!=1:raise ValueError('Multiple steps: set step_name and actual step pressure ramp')
        name=names[0]
    step=odb.steps[name];duration=float(cfg['step_duration'])
    if duration<=0:raise ValueError('step_duration must be positive')
    if hasattr(step,'timePeriod') and abs(step.timePeriod-duration)>1e-8:raise ValueError('Configured duration differs from ODB')
    p0=float(cfg['pressure_at_step_start_mpa']);p1=float(cfg['pressure_at_step_end_mpa'])
    def pressure(t):return p0+(p1-p0)*t/duration
    criteria={'fiber_initiation':([1,2],1.,0.),'matrix_initiation':([3,4],1.,0.),'fiber_completion':([5],float(cfg['fiber_threshold']),tol),'matrix_completion':([6],float(cfg['matrix_threshold']),tol)}
    events={k:None for k in criteria};prev={k:None for k in criteria};checked=0
    for frame in step.frames:
        keys={str(k) for k in frame.fieldOutputs.keys() if str(k).startswith('SDV')}
        if not keys:continue
        if keys!={'SDV'+str(i) for i in range(1,required+1)}:raise ValueError('ODB SDV layout mismatch; use the current paired UMAT/INP')
        checked+=1;t=float(frame.frameValue)
        for label,(indices,threshold,eps) in criteria.items():
            if events[label] is not None:continue
            peak=None
            for index in indices:
                for value in frame.fieldOutputs['SDV'+str(index)].values:
                    if str(value.position)!='INTEGRATION_POINT':raise ValueError('Require unaveraged integration-point SDV')
                    raw=value.dataDouble if str(getattr(value,'precision',''))=='DOUBLE_PRECISION' else value.data
                    try:number=float(raw)
                    except TypeError:number=float(raw[0])
                    if not math.isfinite(number):raise ValueError('Non-finite SDV')
                    if peak is None or number>peak[0]:peak=(number,index,value)
            if peak is None:continue
            if peak[0]>=threshold-eps:
                value=peak[2]
                events[label]={'step_time':t,'pressure_mpa':pressure(t),'previous_below_threshold_pressure_mpa':None if prev[label] is None else pressure(prev[label]),'sdv':'SDV'+str(peak[1]),'value':peak[0],'threshold':threshold,'tolerance':eps,'instance':value.instance.name,'element':value.elementLabel,'integration_point':value.integrationPoint}
            else:prev[label]=t
    if not checked:raise ValueError('No SDV frames found')
    mode=cfg['vessel_failure_rule'];selected=['fiber_completion','matrix_completion'] if mode=='either' else [mode+'_completion']
    candidates=[dict(events[k],criterion=k) for k in selected if events[k] is not None]
    failure=min(candidates,key=lambda x:x['step_time']) if candidates else None
    return {'variant':variant,'step':name,'configuration':cfg,'events':events,'vessel_failure':failure,'status':'THRESHOLD_REACHED' if failure else 'NOT_REACHED_IN_RECORDED_FRAMES','pressure_assumption':'User-configured linear ramp; not inferred from ODB','scope':'Recorded integration-point threshold events; not experimental burst validation'}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('odb');parser.add_argument('--config',default=str(Path(__file__).with_name('failure_thresholds.json')));parser.add_argument('--output');args=parser.parse_args()
    cfg=json.loads(Path(args.config).read_text(encoding='utf-8'))
    from odbAccess import openOdb
    odb=openOdb(args.odb,readOnly=True)
    try:result=scan(odb,cfg)
    finally:odb.close()
    output=Path(args.output) if args.output else Path(args.odb).with_name(Path(args.odb).stem+'_failure.json')
    with output.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(str(output))
if __name__=='__main__':main()
