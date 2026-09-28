from pathlib import Path
import shutil,json,hashlib
r=Path.cwd();backup=r/'audit_20260926/sdv_compact/backup/docs'
for rel in ['README.md','test/README.md','test/A_stiffness_reduction/README.md','test/B_energy_softening/README.md','test/B_energy_softening/AGENT_HANDOFF.md','test/A_stiffness_reduction/RUN_RECORD.md','wcm_hashin_umat/CHANGELOG.md']:
 p=r/rel;q=backup/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
common='''# SDV 精简接口与失效压力判读

更新：2026-09-26。A 为直接刚度折减，B 为能量软化。0.95 的最终损伤阈值仅用于 B。

## 可见输出

| 输出 | A：直接折减 | B：能量软化 |
|---|---|---|
| SDV1 | 纤维拉伸 Hashin 判据 FI_FT | 同左 |
| SDV2 | 纤维压缩 Hashin 判据 FI_FC | 同左 |
| SDV3 | 基体拉伸 Hashin 判据 FI_MT | 同左 |
| SDV4 | 基体压缩 Hashin 判据 FI_MC | 同左 |
| SDV5 | 纤维组合损伤 df | 非黏性纤维组合损伤 df；0.95 为本项目纤维完全失效阈值 |
| SDV6 | 基体组合损伤 dm | 非黏性基体组合损伤 dm；当前上限 0.85 |
| SDV7 | 不输出 | 黏性纤维组合损伤 dfv |
| SDV8 | 不输出 | 黏性基体组合损伤 dmv |

仅输出 S、U 与上表 SDV；不定义 *Depvar 自定义标签，所以 ODB 直接显示 SDV1、SDV2 等。
A 的 *Depvar=18，B 的 *Depvar=60。内部历史不请求为场输出，不能把 *Depvar 改成可见输出个数。

## 两个模式的区别

A 的判据采用旧损伤刚度下的名义试应力；达到 FI>=1 时，对应模式起始标志锁定，按 PROPS(17:20) 的固定损伤量直接折减。当前 FT/FC/MT/MC 分别为 0.07/0.14/0.2/0.4；这些是损伤量，不是剩余刚度比例。纤维单独拉伸触发后 df=0.07，C11 保留约 93%；不能将其解释成损伤已达到 93%。A 没有独立的渐进完全损伤阶段，不使用 0.95。

B 的判据采用未损伤刚度 C0*应变的有效应力。FI>=1 是损伤起始；起始帧的损伤仍可能为 0，之后按能量软化增长。SDV5/6 表示非黏性损伤；SDV7/8 反映黏性正则化后的损伤。当前 update=0 保留原分步时序：本帧更新的损伤用于下一次增量刚度，本帧应力采用上一步损伤。因此 SDV5 达到 0.95 不代表本帧返回应力已使用同等折减。黏性量可能滞后，SDV7 精确等于 0.95 也不是必然，不能据此否认 SDV5 已达到规定阈值。

每个量只在权重非零的 +θ/-θ 两族中取最大值。FI 分模式取最大值；损伤先在每一族内计算 df=1-(1-dFT)(1-dFC)、dm=1-(1-dMT)(1-dMC)，B 再按材料 capf/capm 截顶，然后取两族最大值。黏性组合使用相同公式。禁止先跨族合并拉压损伤再组合，否则会制造不存在于任何一族的损伤。最大值表示至少一族达到，不代表两族全部失效。为精简输出，不提供族身份；需要区分族时须另开诊断输出。

## 如何确定开始失效压力

在复材单元的积分点上逐帧读取最大值，保留单元号、积分点号和帧时间；不要用节点平均后的云图颜色判断阈值。

- **首个纤维失效起始**：第一次 max(SDV1, SDV2)>=1；只关注拉伸断纤时用 SDV1>=1。
- **首个基体失效起始**：第一次 max(SDV3, SDV4)>=1。
- **首个任意模式失效**：第一次 max(SDV1, SDV2, SDV3, SDV4)>=1。
- A 起始后应力折减可能使 FI 再次低于 1，故须扫描历史，不能只看最后一帧。SDV5>0、SDV6>0 分别证明纤维/基体已发生过损伤（当前材料卡损伤系数均非零）。B 的 FI 起始与损伤开始增长也可能分属相邻帧。

当前两套 INP 为单步从零线性加载至 250 MPa，步长 T=1，无另挂载荷幅值，因此 **p=250×帧的 step time MPa**。例：首次纤维 FI>=1 出现在 t=0.60，则首次记录的纤维起始压力为 150 MPa。若上一帧 t=0.59 尚未触发，则事件实际位于 (147.5,150] MPa；报告此区间，不把离散帧当作精确临界点。若改了幅值、多步或非零初压，应按实际幅值和继承载荷重算，不能继续直接乘 250。

## 如何确定最终失效压力

### B：按本项目 SDV5=0.95 阈值

1. 逐帧计算复材积分点 SDV5 的最大值。
2. 找到首次 SDV5>=0.95 的帧；考虑 ODB 单精度，可用 0.95-1e-6，并在结果中注明容差。
3. 以该帧压力作为“**首个局部纤维达到完全失效阈值的压力 P_f,0.95**”，同时报告前一帧压力构成的区间、位置，以及该点 SDV1/2 和 SDV7。若尚未触发过目标模式，不能只凭异常初值判废。
4. 若研究预先规定“任一积分点纤维达到 0.95 即判整瓶失效”，这个 P_f,0.95 就是该约定下的数值失效压力；报告中明确这项局部判废规则。SDV5 是拉压组合，若研究只允许纤维拉伸破坏判废，应结合该点 SDV1/2 的历程确认模式，不能将压缩损伤混算为拉伸爆破。
5. 若采用更强的“整瓶实际爆破”定义，还需检查损伤是否形成关键承载路径的破坏、位移响应及试验对照。仅一个积分点达到上限不足以自动证明实际爆破。

### A：直接折减没有第二个完全损伤阈值

A 用 FI>=1 定位首次模式触发，并用 SDV5/6 确认锁定后的折减。若采用“首个纤维触发即判废”的规则，可将首次 max(SDV1,SDV2)>=1 的压力作为该规则下的失效压力，并明确它是首纤维失效判据。若要得到后续整瓶极限承载压力，A 的某一个 SDV 不会自动提供答案，需依据预先定义的承载路径破坏或全局极限响应判定。当前 A 求解成功至 250 MPa 只说明数值求解完成，不能据此写爆压等于或大于 250 MPa。

两版均不得把最后收敛帧、第一处基体开裂、迭代不收敛直接等同于最终爆破。若计算终止前 B 从未达到 0.95，应报告“在已获得的收敛帧中未达到此阈值”，不要外推最终压力。

## 内部状态与兼容性

- A：1–6 可见；7–10 为正族 FT/FC/MT/MC 起始标志，11–14 为负族对应标志；15–18 为角度、正族权重、模式、schema。
- B：1–8 可见；9–32 为正族，33–56 为负族；每族依次为 D(4)、KAPPA(4)、INIT(4)、DELTA0(4)、SIGMA0(4)、DV(4)；57–60 为角度、正族权重、模式、schema。
- 逐槽索引见各版本 SDV_MAP.csv。这些历史用于不可逆性、能量软化、黏性与材料身份检查，不可为了减少 ODB 项目而删除。
- 核心计算仍使用临时 84 槽数组，已移除持久存储的重复单层应力/应变等诊断量；现有公式、参数和更新时序不变。
- 新 INP 与新 UMAT 必须成对使用。旧 84 项 ODB 不会自动改名，旧 restart 不能接续新布局。A 已完成的旧 ODB 保留，备份子程序位于 audit_20260926/sdv_compact/backup/。

## 本次验证

- 原核心 13 项回归通过。
- 精简历史与 84 槽核心的 504 次编译执行增量对照通过，最大应力差 0；覆盖双版本、两种更新时序、角度、零权重族、卸载和反向加载。
- 最终交付子程序各完成 Abaqus 损伤小模型验证（每版 14 个单元）；真实 ODB 仅有 A 的 SDV1–6、B 的 SDV1–8。
- 各版 131 张材料卡已配对；受保护的材料参数、网格、载荷、边界、取向和求解设置核对一致。B 独立生成器 --check 通过。
- 整瓶新布局 Data Check 状态见 audit_20260926/sdv_compact/datacheck.json；本次未重跑整瓶求解。A 的旧布局整瓶已成功计算至 250 MPa，B 尚无整瓶完整求解验证。
'''
(r/'wcm_hashin_umat/SDV说明.md').write_text(common,encoding='utf-8')
for v,d,name in [('A','A_stiffness_reduction','hashin_constant_wcm.for'),('B','B_energy_softening','hashin_energy_wcm.for')]:
 f=r/'test'/d
 (f/'SDV说明.md').write_text(common,encoding='utf-8')
 mode='直接刚度折减' if v=='A' else '能量软化'
 (f/'README.md').write_text(f'''# {v}：{mode}

当前运行配对：`Job2_PDM_{v}.inp` + `{name}`。

输出为 S、U、{'SDV1–SDV6' if v=='A' else 'SDV1–SDV8'}；*Depvar={'18' if v=='A' else '60'}，内部历史不请求为场输出。

**[SDV 含义、开始失效与最终失效压力判读](SDV说明.md)**；逐槽映射见 [SDV_MAP.csv](SDV_MAP.csv)。

最终精简版通过编译回归与真实 Abaqus 损伤小模型；整瓶 Data Check 证据见项目 audit_20260926/sdv_compact/datacheck.json。本次未重跑整瓶。

''' + ('旧布局 A 整瓶已完成 105 增量至 250 MPa，原 ODB 与日志保留；详见 [RUN_RECORD.md](RUN_RECORD.md)。该 ODB 对应旧 84 项布局，不能按新 SDV 序号解释。原 UMAT 已备份。\n' if v=='A' else 'B 尚无整瓶完整求解验证。独立生成/核对工具已同步，可运行 `python scripts/prepare_case.py --check`；生成新候选用 `--name 新作业名`，产物在 runs/。\n'),encoding='utf-8')
 if v=='B':
  (f/'AGENT_HANDOFF.md').write_text('''# B 能量软化：当前交接

以 README.md 和 SDV说明.md 为当前接口依据。主配对 Job2_PDM_B.inp + hashin_energy_wcm.for；33 个材料常数、60 个内部状态，仅输出 SDV1–8 与 S/U。B 的 SDV5 达到 0.95 为约定的局部纤维完全失效阈值。

scripts/prepare_case.py --check 已通过；--name 可生成独立候选。源几何/角度/材料配置保留，source 下 UMAT 模板与计算核心已同步；更新前快照保存在项目 audit_20260926/sdv_compact/backup/B_support。

每次整瓶运行在新 runs 子目录中进行。当前编译、小模型及状态映射验证不能代替整瓶求解或实际爆压验证。保留现有参数、取向、网格、载荷和本构时序；详细证据见项目 audit_20260926/sdv_compact/。
''',encoding='utf-8')
p=r/'test/A_stiffness_reduction/RUN_RECORD.md';s=p.read_text(encoding='utf-8');s=s.replace('## 运行目录','> 历史记录：本记录和已有 ODB 对应旧 84 项布局。本轮根目录 UMAT/INP 已升级为精简版，不能再按下文历史哈希认定当前文件。旧 UMAT 备份见 `../../audit_20260926/sdv_compact/backup/test/A_stiffness_reduction/`；原计算结果未修改。\n\n## 运行目录',1);p.write_text(s,encoding='utf-8')
p=r/'README.md';s=p.read_text(encoding='utf-8');a=s.index('## 当前测试入口');b=s.index('## 原项目介绍');s=s[:a]+'''## 当前进度与运行入口

- A 直接刚度折减：旧 84 项布局整瓶已成功完成 105 增量至 250 MPa。
- B 能量软化：尚未完成整瓶求解。
- 两版已精简 SDV：A 输出 SDV1–6，B 输出 SDV1–8；分别使用 18/60 个内部状态。新旧布局不可混用。
- 精简版已通过原核心回归、504 次状态映射对照和两版真实 Abaqus 损伤小模型；本次未重跑整瓶。

[A 运行包](test/A_stiffness_reduction/README.md) · [B 运行包](test/B_energy_softening/README.md) · [SDV 与失效压力判读](wcm_hashin_umat/SDV说明.md) · [测试入口](test/README.md)

''' +s[b:];p.write_text(s,encoding='utf-8')
(r/'test/README.md').write_text('''# 两版本运行包

| 版本 | 输入 / 子程序 | 可见 SDV | 内部状态 |
|---|---|---|---|
| [A 直接刚度折减](A_stiffness_reduction/README.md) | Job2_PDM_A.inp / hashin_constant_wcm.for | 1–6 | 18 |
| [B 能量软化](B_energy_softening/README.md) | Job2_PDM_B.inp / hashin_energy_wcm.for | 1–8 | 60 |

A 原布局整瓶已经计算成功至 250 MPa，日志与 ODB 保留在 A 根目录；B 尚无完整整瓶求解。当前精简版已完成两版真实 Abaqus 损伤小模型；本次未重跑整瓶。整瓶 Data Check、回归及交付哈希见 ../audit_20260926/sdv_compact/。

只输出 S、U 和指定 SDV。SDV1–4 为四模式判据，SDV5/6 为纤维/基体损伤，B 的 SDV7/8 为黏性损伤。0.95 阈值仅用于 B 的 SDV5，A 按触发判据直接折减。

完整使用说明与压力判读见各套 SDV说明.md。旧结果不随源码修改而改名；禁止用新布局解释原 84 项结果或进行旧 restart。

B 的独立脚本、源模板、状态字典和结果清单已经同步，scripts/prepare_case.py --check 可核对；A 保留简洁的 INP/UMAT 配对，旧生成工具仍在历史 Git 中。本轮替换前文件已备份到 ../audit_20260926/sdv_compact/backup/。

历史初始化记录 results/packages.json、results/delivery_checks.json 对应旧包，不作为当前 SDV 接口依据；当前清单为 ../audit_20260926/sdv_compact/delivery.json。
''',encoding='utf-8')
p=r/'wcm_hashin_umat/CHANGELOG.md';s=p.read_text(encoding='utf-8');p.write_text('# 2026-09-26 SDV 精简\n\nA 持久状态 84→18、可见 6；B 84→60、可见 8。保留四模式 FI、组合纤维/基体损伤，B 另保留黏性组合损伤。移除自定义 Depvar 标签，明确仅 B 使用 0.95 阈值。504 次状态恢复回归及两版真实 Abaqus 损伤模型通过。详见 SDV说明.md。\n\n'+s,encoding='utf-8')
# Update final delivery hashes after independent B regeneration.
p=r/'audit_20260926/sdv_compact/delivery.json';m=json.loads(p.read_text())
for v,d in [('A','A_stiffness_reduction'),('B','B_energy_softening')]:
 for n in m[v]['files']:m[v]['files'][n]=hashlib.sha256((r/'test'/d/n).read_bytes()).hexdigest()
p.write_text(json.dumps(m,indent=2),encoding='utf-8')
