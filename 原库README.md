# 3D Hashin Type IV 气瓶渐进损伤

## 当前版本（2026-09-27）

计算基准已切换到朱宁 rebuild 模型 `geometry_zhuning_joint_rebuild`（WCM 单 part、135 个等效 Bin、51,558 个复材单元）；旧 `geometry_zhuning_10deg` 的 Job-2 等基准及其结果全部废弃，旧编号的 ODB/restart 不得与新接口混用。

A 直接折减输出 8 项；B 能量软化输出 10 项。SDV5/6 分别表示纤维/基体失效标志或非黏性损伤，默认 1 为阈值；B 同时输出新黏性损伤与本帧实际损伤。刚度保留 1e-6 残余因子。两套场输出为 LE、PE、PEEQ、PEMAG、S、U 加编号 SDV，*Static 取 0.01, 1, 1e-20, 0.1。只有 B 的复材截面加 `hourglass=ENHANCED`（A 按本轮口径不加）。

[两版运行入口](test/README.md) · [SDV完整说明](wcm_hashin_umat/SDV说明.md) · [PROPS逐项说明](wcm_hashin_umat/PROPS说明.md) · [待测包生成与校验](wcm_hashin_umat/tools/prepare_test_cases.py)

两版 UMAT 均已编译并通过 Abaqus 小模型；新整瓶 INP 的 B 版 Data Check 通过，A 版因 C3D8R 沙漏刚度为零而失败。B 正式版整瓶在约 0.170 步时间停止，`DFINAL≤DEL0` 直接置损伤 1 的试验版在 0.4428 停止，详见 [试跑证据](audit_20260927/instant_fallback/README.md)。当前代码已恢复正式版，两版均无完整整瓶求解结果；此前旧基准的 A 整瓶结果不作为当前版结论。

## 原项目介绍
Progressive failure analysis of Type IV composite hydrogen storage cylinders based on 3D Hashin criteria and energy-method linear softening for stiffness degradation — from China University of Petroleum (East China)

more information about the umat subroutine（helped by codex） and the inp introduction in readme.pdf

Because I am also a graduate student working on this topic, it’s hard to avoid having some questions.

if any question about this project,send related questions to the email Z25150070@s.upc.edu.cn
