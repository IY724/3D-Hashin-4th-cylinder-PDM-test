# B 版待测包（断裂能线性软化 + 黏性正则化）

配对：`Rebuild_PDM_B.inp` + `hashin_energy_wcm.for`。本目录只保留这两个待测文件与本说明，求解产物提交时另建 `runs/<作业名>/`，不与本目录混放。

## 计算基准

- 源输入：`../../geometry_zhuning_joint_rebuild/WCM_Job.inp`（朱宁 rebuild 模型）。旧 `geometry_zhuning_10deg` 的 Job-2/Job-multi/WCM_multi 基准连同其结果已全部废弃。
- 复材 part `Tank-1`：51,558 个单元（C3D8R 50,916、C3D8 636、C3D6 6），135 个 WCM 等效 Bin 全部换成 `*User Material`（constants=30, unsymm）。
- 沙漏控制：C3D8R 配 UMAT，本套对 391 个复材截面挂 `*Section Controls, name=HG_UMAT, hourglass=ENHANCED`；A 版按本轮口径不加，两版单元控制不同，跨版对照时须计入该差异。Enhanced 只保证无需手填沙漏刚度，不等于精度保证，仍须查变形模式与 ALLAE。
- 内衬 `PA6`、阀座 `al-6061` 与备用单层 `t-700` 卡片逐字保留，仍用内置本构。
- 分析步与 A 完全相同：`*Static` 0.01、1、1e-20、0.1；`*Dsload, Surf-2, P, 250.`（p=250×t MPa）；约束与循环对称未改动。两套 INP 除 Bin 材料卡与上述沙漏控制外逐字一致（生成脚本以保护摘要证明）。
- 场输出：`*Node Output` = U；`*Element Output, directions=YES` = LE, PE, PEEQ, PEMAG, S, SDV1–SDV10。UVARM 场请求与 `*User Output Variables` 已移除。

## 接口

| 项 | 值 |
|---|---|
| PROPS 常数 | 30 项（前 9 项为原始单层工程常数；PROPS(24/25)=1e-4/5e-3 为纤维/基体黏性时间；PROPS(28/29)=1/1 只限制法向刚度折减上限） |
| `*Depvar` | 62（内部状态，±θ 两族各自保存 D、KAPPA、INIT、DELTA0、SIGMA0、DV） |
| 场输出 SDV | SDV1–SDV10 |

| 输出 | 含义 |
|---|---|
| SDV1–4 | 纤维拉伸/纤维压缩/基体拉伸/基体压缩 Hashin 判据（两族取最大） |
| SDV5 / SDV6 | 非黏性纤维 / 基体组合损伤，0→1，达到 1 表示软化完成 |
| SDV7 / SDV8 | 本增量更新后的黏性纤维 / 基体组合损伤 |
| SDV9 / SDV10 | 本帧应力实际采用的纤维 / 基体组合损伤系数 |

判废看 SDV5/6；SDV7–10 用于解释黏性滞后与本帧刚度。当前 `update=0`，本帧应力用旧黏性损伤，故 SDV9/10 与 SDV7/8 相差一个增量。LE 为对数应变；UMAT 不返回非弹性应变，PE/PEEQ/PEMAG 预期恒为 0。完整槽位含义见 `../../wcm_hashin_umat/SDV说明.md`。

## 提交与后处理

```bat
abaqus job=Rebuild_PDM_B user hashin_energy_wcm.for cpus=8 interactive
abaqus python ../../wcm_hashin_umat/tools/extract_failure.py runs/<作业名>/Rebuild_PDM_B.odb --config ../../wcm_hashin_umat/config/failure_thresholds_B.json
```

先做 datacheck 再正式提交；阈值与判废规则在 `failure_thresholds_B.json` 中改，默认纤维/基体阈值均为 1、规则 either。黏性影响随 η/Δt 变化，最大增量 0.1 时基体 η=5e-3 的权重明显低于初始增量 0.01，故步长敏感的起始量结论须做增量敏感性对照。

## 生成与校验

`python ../../wcm_hashin_umat/tools/prepare_test_cases.py`（加 `--check` 只复算比对、不写文件）。离线证据写入 `../../wcm_hashin_umat/validation/test_cases.json`。UMAT 已通过编译、材料点回归及 Abaqus 小模型；整瓶 Data Check 通过。当前正式版整瓶曾因 `DFINAL≤DEL0` 在约 0.170 步时间停止；直接置损伤 1 的实验版又在 0.4428 停止，且大范围触发即时损伤。实验版及证据见 `../../audit_20260927/instant_fallback/`，当前文件已恢复正式版，尚无完整整瓶求解。

## 材料参数读法

当前子程序在 `WCM_CHECK` 中逐项读取并命名 PROPS；完整槽位、单位和 Bin0 实际数值见 [PROPS说明](../../wcm_hashin_umat/PROPS说明.md)。已去掉材料卡中的损伤开关、schema 和 model_id，当前材料卡分别为 A25/B30 项。`UPDATE` 及 B 的刚度折减上限仍可配置。
