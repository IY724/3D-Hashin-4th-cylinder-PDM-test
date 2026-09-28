from pathlib import Path
import shutil,json
r=Path.cwd();audit=r/'audit_20260927/damage_one'
def write(p,s):
 q=audit/'backup_docs'/p.relative_to(r);q.parent.mkdir(parents=True,exist_ok=True)
 if p.exists() and not q.exists():shutil.copy2(p,q)
 p.write_text(s,encoding='utf-8')
common='''## 怎么判断失效和压力

- **损伤起始**：SDV1–4 中对应模式第一次达到 1；纤维看 SDV1/2，树脂基体看 SDV3/4。
- **纤维失效**：SDV5 达到 `fiber_threshold`，默认 **1**。
- **树脂基体失效**：SDV6 达到 `matrix_threshold`，默认 **1**。
- 所有输出均在权重非零的 ±θ 族中取最大值；任一角度族达到即可，不等待另一族。组合损伤先在每族内算 `1−(1−d拉)(1−d压)`，再取两族最大值。螺旋层保留独立内部历史，环向 ±90° 表示同一纤维轴向。
- 按此次“纤维或树脂任一达到即判废”的要求，`vessel_failure_rule="either"`；可改为 `"fiber"` 或 `"matrix"`。这是本项目可配置判废规则，实际爆压需试验对照。

后处理参数在本目录 **failure_thresholds.json**，两个阈值分别设置，默认均为 1，数值比较容差为 1e-6。阈值可降低，但这表示提前判废，不能把阈值 0.8 描述为损伤已经演化完毕。

运行：`abaqus python extract_failure.py 新结果.odb`。脚本只读 ODB，输出 `新结果_failure.json`，包含纤维/基体起始、达到阈值的压力、上一未触发帧压力、单元和积分点位置及按选定规则得到的判废压力。仅检查积分点原始值，不使用节点平均云图；无达到阈值的帧则报告未达到，不把最后收敛压力当作爆压。输出文件已存在时拒绝覆盖。

当前单步从零线性加载到 250 MPa、步时长1：**p=250×帧时间 MPa**。例如 t=0.60 对应150 MPa。输出给出离散帧的压力区间；改变幅值、步或载荷后必须同步修改后处理配置，多步需指定 step_name。脚本不从 ODB 自动推断实际压力载荷。

## 残余刚度与兼容性

`WCM_STIFF` 的 `RMIN=1e-6` 对法向组合折减因子与剪切折减因子设置下限。因此损伤/标志可以显示1，材料弹性刚度不会直接归零；这是刚度因子的下限，不是把损伤输出强行限制为0.999999。非零残余刚度不保证全局雅可比良态或非线性一定收敛，软化、黏性和增量仍会影响求解。

新 INP/UMAT 必须配对；旧 18/60/84 项以及任何旧编号 ODB/restart 不兼容，不用旧结果解释新 SDV。本次未提交整瓶完整求解。
'''
a='''# A 版 SDV：直接刚度折减

2026-09-27 当前接口：Job2_PDM_A.inp + hashin_constant_wcm.for。输出 S、U、**SDV1–8**，*Depvar=20。

| 输出 | 含义 |
|---|---|
| SDV1 | 纤维拉伸 Hashin 判据 |
| SDV2 | 纤维压缩 Hashin 判据 |
| SDV3 | 基体拉伸 Hashin 判据 |
| SDV4 | 基体压缩 Hashin 判据 |
| SDV5 | 纤维失效标志：0未触发，1已触发 |
| SDV6 | 基体失效标志：0未触发，1已触发 |
| SDV7 | 本帧应力实际采用的纤维组合损伤系数 |
| SDV8 | 本帧应力实际采用的基体组合损伤系数 |

A 不引入黏性或连续软化。SDV5/6 是不可逆0/1标志，不是刚度损失比例。触发后仍按现有 FT/FC/MT/MC 固定损伤量0.07/0.14/0.2/0.4折减。当前 update=0，本帧新触发标志从下一增量改变刚度；因此 SDV5=1 时 SDV7 可以仍为0。下一增量仅纤维拉伸触发时 SDV7=0.07，表示该族纵向刚度因子约为0.93，不能称纵向刚度已经丧失。

SDV1–4 来自旧损伤下名义试应力，失效后可能下降；SDV5/6 保留已触发记录。失效压力按标志阈值判断，无需从末帧 FI 反推。

'''+common+'''
## 内部历史与验证

9–12/13–16：正/负族 FT、FC、MT、MC 起始标志；17–20：角度、权重、模式、schema。完整槽位见 SDV_MAP.csv。内部项不请求为场输出。

原A固定折减响应在残余下限未介入时保持一致。当前版通过编译回归、状态恢复对照及真实 Abaqus 失效小模型。此前旧版 A 整瓶已完成至250 MPa，见 RUN_RECORD.md，该记录不代表当前版已重跑整瓶。
'''
b='''# B 版 SDV：能量软化与黏性滞后

2026-09-27 当前接口：Job2_PDM_B.inp + hashin_energy_wcm.for。输出 S、U、**SDV1–10**，*Depvar=62。

| 输出 | 含义 |
|---|---|
| SDV1 | 纤维拉伸 Hashin 判据 |
| SDV2 | 纤维压缩 Hashin 判据 |
| SDV3 | 基体拉伸 Hashin 判据 |
| SDV4 | 基体压缩 Hashin 判据 |
| SDV5 | 非黏性纤维组合损伤：0→1，达到1表示软化完成 |
| SDV6 | 非黏性基体组合损伤：0→1，达到1表示软化完成 |
| SDV7 | 本增量更新后的黏性纤维组合损伤 |
| SDV8 | 本增量更新后的黏性基体组合损伤 |
| SDV9 | 本帧应力实际采用的纤维组合损伤系数 |
| SDV10 | 本帧应力实际采用的基体组合损伤系数 |

**判废看 SDV5/6；SDV7–10 用于解释黏性滞后与本帧刚度。**非黏性损伤由原能量线性软化公式演化，现允许达到1，不再在0.95/0.85截断。FI=1仅表示损伤起始，通常不是软化完成。

每模式黏性更新为 `dv新=(η·dv旧+Δt·d新)/(η+Δt)`。当前 update=0，本帧应力用 dv旧，因此 SDV9/10 对应旧黏性损伤，SDV7/8 对应新黏性损伤；update=1 时本次反馈，二者在未受刚度上限截顶时一致。η=0 时新黏性损伤等于新非黏性损伤，但 update=0 仍有一个增量的反馈滞后。

真实编译材料点例：纤维非黏性损伤首次达到1时，新黏性损伤约0.911，本帧实际采用约0.866。不要等 SDV7 精确等于1才认为 SDV5 的软化完成；黏性量可能渐近趋于1。

PROPS(28/29) 现在仅限制法向组合刚度折减，当前材料卡均改为1；不再限制非黏性/黏性损伤历史。保留0.95/0.85可保留更多法向刚度，但不是恢复旧版本的完整响应，因为历史演化已改变，剪切仍按各模式黏性损伤计算。当前程序对组合刚度另设1e-6下限。材料参数上限与 failure_thresholds.json 的后处理阈值相互独立。

'''+common+'''
## 内部历史与验证

11–34/35–58：正/负族历史；每族依次为D、KAPPA、INIT、DELTA0、SIGMA0、DV（各4项，FT/FC/MT/MC顺序）；59–62：角度、权重、模式、schema。完整槽位见 SDV_MAP.csv。

当前版通过编译回归、能量/黏性公式独立核对、完整损伤残余刚度检查及真实 Abaqus 失效小模型。小模型采用专门测试断裂能以到达软化终点，不代表整瓶材料标定。B 尚未完成整瓶求解，当前变化可能改变原有收敛性和预测压力。独立核对：`python scripts/prepare_case.py --check`。
'''
for d,t in [('A_stiffness_reduction',a),('B_energy_softening',b)]:
 folder=r/'test'/d;write(folder/'SDV说明.md',t)
 v=d[0];name='hashin_constant_wcm.for' if v=='A' else 'hashin_energy_wcm.for'
 write(folder/'README.md',f'''# {v} 版当前运行入口

配对：Job2_PDM_{v}.inp + {name}。输出 S、U、SDV1–{'8' if v=='A' else '10'}；内部状态{'20' if v=='A' else '62'}项。

[SDV含义与失效压力判断](SDV说明.md) · [完整槽位](SDV_MAP.csv) · [后处理阈值](failure_thresholds.json)

SDV1–4为起始判据；SDV5/6为纤维/基体{'失效标志' if v=='A' else '非黏性损伤'}，默认达到1判该模式失效。{'SDV7/8给出本帧采用的损伤。' if v=='A' else 'SDV7/8为新黏性损伤，SDV9/10为本帧采用的损伤。'}任一角度族达到即可；默认纤维或基体任一达到就按项目规则判废，配置可改。

`abaqus python extract_failure.py 新结果.odb` 生成失效压力记录。仅适用于当前配对的新结果。

当前版已通过编译、回归和 Abaqus 失效小模型；整瓶 Data Check 记录见 ../../audit_20260927/damage_one/datacheck.json。本次未重跑整瓶。旧输出编号不再适用。
''')
write(r/'wcm_hashin_umat/SDV说明.md','# 当前 SDV 接口\n\n'+a+'\n---\n\n'+b)
write(r/'test/B_energy_softening/AGENT_HANDOFF.md','''# B 当前交接

Job2_PDM_B.inp + hashin_energy_wcm.for；33个常数、62个内部状态、SDV1–10。以 SDV说明.md 为准。非黏性损伤完整演化到1；保留新黏性损伤及本帧实际采用的损伤；刚度因子下限1e-6，当前P28/P29=1。

scripts/prepare_case.py --check 核对独立包，--name 生成新配对。后处理使用包根 extract_failure.py 与 failure_thresholds.json；新候选目录需复制这两个文件或明确指定配置路径。计算产物放独立 runs 目录。旧B求解响应不可假定沿用；当前完成材料点、小模型及Data Check验证，未进行整瓶完整求解。证据见 ../../audit_20260927/damage_one/。
''')
write(r/'test/README.md','''# 当前两版测试入口

|版本|输出|内部状态|主要判读|
|---|---|---|---|
|[A 直接折减](A_stiffness_reduction/README.md)|SDV1–8|20|SDV5/6=1为对应模式触发；7/8为本帧采用的损伤|
|[B 能量软化](B_energy_softening/README.md)|SDV1–10|62|SDV5/6=1为软化完成；7/8为新黏性损伤；9/10为本帧采用的损伤|

四模式起始判据均在SDV1–4。各版输出S/U和指定SDV；损伤完成指标与非零残余刚度分开。阈值与判废规则在各包 failure_thresholds.json，可通过 extract_failure.py 输出压力记录。

2026-09-27本轮证据为 ../audit_20260927/damage_one/，旧sdv_compact与sdv_definition为历史。A旧整瓶计算已完成至250 MPa，当前两版未重跑整瓶。旧84、18/60项结果和新接口不可混用。
''')
p=r/'README.md';s=p.read_text(encoding='utf-8');original=s[s.index('## 原项目介绍'):].split('## 2026-09-27 SDV 判读约定')[0]
write(p,'# 3D Hashin Type IV 气瓶渐进损伤\n\n## 当前版本（2026-09-27）\n\nA直接折减输出8项；B能量软化输出10项。SDV5/6分别表示纤维/基体失效标志或非黏性损伤，默认1为阈值；B同时输出新黏性损伤与本帧实际损伤。刚度保留1e-6残余因子。\n\n[两版运行入口](test/README.md) · [SDV完整说明](wcm_hashin_umat/SDV说明.md)\n\nA旧版整瓶已完成至250 MPa；当前版完成回归、失效小模型和Data Check验证，本次未重跑整瓶。\n\n'+original)
p=r/'wcm_hashin_umat/CHANGELOG.md';old=p.read_text(encoding='utf-8');write(p,'''# 2026-09-27 损伤完成值与黏性输出

当前A输出8项/20状态，B输出10项/62状态。旧条目为历史，不作为当前接口。

- B非黏性/黏性损伤历史允许到1；P28/P29只限制法向刚度折减，默认改为1/1。与旧0.95/0.85截断历史的响应不同。
- 刚度法向组合与剪切因子下限1e-6，防止损伤完成后归零；不保证全局收敛。
- A保留固定折减参数，以0/1标志单独输出失效；两版增加本帧实际采用的损伤，B另保留新黏性损伤。
- 新增可配置的积分点阈值后处理；默认纤维/基体阈值1，规则either。旧restart不兼容。
- 证据见 ../audit_20260927/damage_one/；未重跑整瓶。

---

'''+old)
write(audit/'更新记录.md','''# 2026-09-27 损伤完成与黏性输出交付

本轮是算法修改，不是重命名。B能量软化损伤继续到1，黏性历史独立保留；旧0.95/0.85上限改为法向刚度折减上限，当前正式卡取1/1。两版刚度因子下限1e-6。A固定折减参数保持。

源模板、核心、工具、dist、两版INP/UMAT、B独立生成工具、槽位CSV、README、SDV说明、后处理配置与脚本已同步。旧文件在backup和backup_docs。旧84、18/60状态布局不可接续。

验证：13项编译核心回归（B与旧版仅在旧截顶前保持一致）；504次新核心与紧凑历史恢复对照应力差0；独立能量/黏性递推和残余刚度检查；A/B各8个单元的真实Abaqus完成失效试验；后处理读取真实试验ODB通过。该试验B使用专用断裂能，不能充当整瓶物理标定。整瓶Data Check记录在datacheck.json。本次未求解整瓶。
''')
# Point active verifier to current delivery instead of modifying historical evidence.
p=r/'test/scripts/verify_delivery.py';s=p.read_text(encoding='utf-8').replace('audit_20260926/sdv_compact','audit_20260927/damage_one');write(p,s)
