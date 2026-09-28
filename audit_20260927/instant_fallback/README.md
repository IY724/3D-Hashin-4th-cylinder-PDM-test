# B 版即时完成分支的整瓶试跑（2026-09-27）

按用户要求，在 `DFINAL≤DEL0` 时令对应模式非黏性损伤立即为 1，保留原黏性滞后及残余刚度。未更改网格、INP、PROPS、SDV 数量。试验代码见 `experimental_wcm_core.for`、`experimental_umat_entry.for.in` 和 `experimental_hashin_energy_wcm.for`；提交用的精确副本、输入及完整 Abaqus 文件保存在 `../../test/B_energy_softening/runs/B_20260927_153131_fallback/`。日常使用的 `wcm_hashin_umat/src/` 与 `test/B_energy_softening/` 已恢复为试验前版本。

验证：13 项材料点测试通过，A/B Abaqus 小模型通过。整瓶作业编译通过并越过旧作业的 0.170 时间点，但最终只收敛到 0.4428125（线性载荷下 110.7031 MPa）；随后出现大量 `1022` 和 `TOO MANY ATTEMPTS`，整瓶求解未完成。

试算调用日志共有 709,867 条 `WCM_FALLBACK`，涉及 28,393 个不同单元；包括被拒绝增量及重复 Newton 调用，不能当作最终损伤单元数。最后收敛 ODB 帧中，56,016 个含 SDV 的积分点有 14,378 个 `SDV6=1`。只读部分 ODB 的阈值扫描给出基体 SDV6 首达 1 的帧压力 35.4687 MPa，这是即时完成分支造成的候选阈值，**不能作为可信爆压**。

其后的 20,857 条 `1022` 是纤维压缩起始等效量不合法的试算日志。源码中压缩 Hashin 判据由 `S1<0` 触发，而压缩等效位移取 `LC*max(-E1,0)`；当多轴试应力为压缩但该应变非负时，等效位移和等效应力为零。此问题不会仅因缩小 LC 而消失。

详见 `run_summary.json`、`sdv_counts.json`、`partial_failure.json`。重构 B 版软化量定义或更改断裂能含义前，需要核对材料标定；单纯放开上述错误继续计算会让压力阈值失去解释性。
