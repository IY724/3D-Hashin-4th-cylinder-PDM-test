# 当前两版测试入口

计算基准统一为朱宁 rebuild 模型 `../geometry_zhuning_joint_rebuild/WCM_Job.inp`（135 个 WCM 等效 Bin、51,558 个复材单元、静压 0→250 MPa）。旧 `geometry_zhuning_10deg` 的 Job-2 等基准及其结果已全部废弃，旧 84/18/60 项编号的 ODB 与 restart 不可与新接口混用。

|版本|待测文件|输出|内部状态|主要判读|
|---|---|---|---|---|
|[A 直接折减](A_stiffness_reduction/README.md)|Rebuild_PDM_A.inp + hashin_constant_wcm.for|SDV1–8|20|SDV5/6=1 为对应模式触发；7/8 为本帧采用的损伤|
|[B 能量软化](B_energy_softening/README.md)|Rebuild_PDM_B.inp + hashin_energy_wcm.for|SDV1–10|62|SDV5/6=1 为软化完成；7/8 为新黏性损伤；9/10 为本帧采用的损伤|

四模式起始判据均在 SDV1–4。两套场输出同为 LE、PE、PEEQ、PEMAG、S、U 加各自编号 SDV；UMAT 不返回非弹性应变，PE/PEEQ/PEMAG 预期为 0。两套 INP 由同一源 deck 生成，除 Bin 材料卡外逐字一致；差异只在 B 版复材截面加 `hourglass=ENHANCED`（A 版按本轮口径不加），以及 PROPS 项数（A 25 / B 30）与 `*Depvar`（A 20 / B 62）。

每套目录只保留 UMAT、配套 INP 和本说明。生成与离线校验：`python ../wcm_hashin_umat/tools/prepare_test_cases.py`（`--check` 只复算比对），证据在 `../wcm_hashin_umat/validation/test_cases.json`；后处理脚本与阈值配置在 `../wcm_hashin_umat/tools/extract_failure.py` 与 `../wcm_hashin_umat/config/failure_thresholds_{A,B}.json`。

2026-09-27 取证见 `../audit_20260927/damage_one/`、`../audit_20260927/props_layout/` 与 `../audit_20260927/instant_fallback/`。两版 UMAT 已编译并通过小模型；整瓶 B 版 Data Check 通过，A 版报 C3D8R 零沙漏刚度。B 正式版与即时完成试验版都已提交过整瓶，但均未完成；当前 test 包恢复正式版。

两版材料参数逐项含义见 [PROPS说明](../wcm_hashin_umat/PROPS说明.md)。
