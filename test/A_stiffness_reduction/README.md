# A 版待测包（直接刚度折减）

配对：`Rebuild_PDM_A.inp` + `hashin_constant_wcm.for`。本目录只保留这两个待测文件与本说明，求解产物提交时另建 `runs/<作业名>/`，不与本目录混放。

## 计算基准

- 源输入：`../../geometry_zhuning_joint_rebuild/WCM_Job.inp`（朱宁 rebuild 模型，单 part 合并内衬+阀座+WCM 复材）。旧 `geometry_zhuning_10deg` 的 Job-2/Job-multi/WCM_multi 基准连同其结果已全部废弃，不再作为输入或对照。
- 复材 part `Tank-1`：51,558 个单元（C3D8R 50,916、C3D8 636、C3D6 6），135 个 WCM 等效 Bin 全部换成 `*User Material`（constants=25, unsymm），391 个复材截面。
- 内衬 `PA6`、阀座 `al-6061` 与备用单层 `t-700` 卡片逐字保留，仍用内置本构。
- 分析步：`*Step, nlgeom=YES, inc=1000000` + `*Static` 初值 0.01、步时长 1、最小 1e-20、最大 0.1；载荷 `*Dsload, Surf-2, P, 250.`（线性升压，p=250×t MPa）；36 份循环对称与位移约束未改动。
- 场输出：`*Node Output` = U；`*Element Output, directions=YES` = LE, PE, PEEQ, PEMAG, S, SDV1–SDV8。WCM 的 UVARM 场请求与 135 处 `*User Output Variables` 已移除。
- 沙漏控制：A 版当前不加 `*Section Controls`。大模型 Data Check 因 50,916 个 C3D8R 单元的沙漏刚度为零而终止（`ErrElemZeroHourGlassStiffness`）；当前 INP 不能直接用于整瓶求解。

## 接口

| 项 | 值 |
|---|---|
| PROPS 常数 | 25 项（前 9 项为原始单层工程常数，非 Bin 等效值） |
| `*Depvar` | 20（内部状态，含 ±θ 两族独立历史） |
| 场输出 SDV | SDV1–SDV8 |

| 输出 | 含义 |
|---|---|
| SDV1–4 | 纤维拉伸/纤维压缩/基体拉伸/基体压缩 Hashin 判据（两族取最大） |
| SDV5 / SDV6 | 纤维 / 基体失效标志，0 未触发、1 已触发 |
| SDV7 / SDV8 | 本帧应力实际采用的纤维 / 基体组合损伤系数 |

LE 为对数应变；UMAT 不返回非弹性应变，PE/PEEQ/PEMAG 预期恒为 0，仅保留以便与 WCM 弹性基准同口径对照。完整槽位含义见 `../../wcm_hashin_umat/SDV说明.md`。

## 提交与后处理

```bat
abaqus job=Rebuild_PDM_A user hashin_constant_wcm.for cpus=8 interactive
abaqus python ../../wcm_hashin_umat/tools/extract_failure.py runs/<作业名>/Rebuild_PDM_A.odb --config ../../wcm_hashin_umat/config/failure_thresholds_A.json
```

先做 datacheck 再正式提交；阈值与判废规则在 `failure_thresholds_A.json` 中改，默认纤维/基体阈值均为 1、规则 either（任一达到即判废）。

## 生成与校验

`python ../../wcm_hashin_umat/tools/prepare_test_cases.py`（加 `--check` 只复算比对、不写文件）。离线证据写入 `../../wcm_hashin_umat/validation/test_cases.json`：Bin 数、单元数、PROPS/Depvar 数量、SDV 清单、逐 Bin 角度与弹性往返误差、保护摘要。UMAT 已通过编译、材料点回归及 Abaqus 小模型；整瓶 Data Check 报上述沙漏刚度错误，尚未提交整瓶求解。

## 材料参数读法

当前子程序在 `WCM_CHECK` 中逐项读取并命名 PROPS；完整槽位、单位和 Bin0 实际数值见 [PROPS说明](../../wcm_hashin_umat/PROPS说明.md)。已去掉材料卡中的损伤开关、schema 和 model_id，当前材料卡分别为 A25/B30 项。`UPDATE` 及 B 的刚度折减上限仍可配置。
