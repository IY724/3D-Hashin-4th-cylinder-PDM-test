# 3D Hashin IV 型气瓶渐进损伤 · 测试库（test）

本仓库是 [IY724/3D-Hashin-4th-cylinder-PDM](https://github.com/IY724/3D-Hashin-4th-cylinder-PDM) 的**独立公开测试库**。原库继续作为工作库，这里按三类内容归档，便于下载复现与审查：

1. **旧的代码案例说明** —— 原始单层 3D Hashin 子程序及其整瓶输入与说明；
2. **新的渐进损伤代码及其修改说明** —— WCM 等效层 + 双角度 Hashin 的 A/B 两套 UMAT、逐项变更记录；
3. **已做的案例（重点为基于断裂能的渐进损伤）** —— rebuild 基准、整瓶与加密网格的能量法试跑、网格相关性研究。

目录即原工程路径，筛选后结构保持一致，文内相对链接全部可用。

---

## 一、旧的代码案例说明（原始单层 Hashin 阶段）

项目起点：作者自写的**单层 3D Hashin UMAT** 与配套整瓶输入卡，不含 WCM 等效层与双角度处理。

| 路径 | 说明 |
|---|---|
| [readme.pdf](./readme.pdf) | 原始交付说明：UMAT 接口与 INP 使用方式（英文，作者自述） |
| [hashin_self/](./hashin_self/) | 旧版单层 Hashin 子程序 4 份：`hashin.for`、`hashin2.for`、`hashin_constant.for`（固定刚度折减）、`onlyhahsin.for` |
| `Job-hashinnewcdm.inp` | 旧案例一：CDM 能量法整瓶输入（根目录，18.8 MB） |
| `Job-hashinnewconstant.inp` | 旧案例二：固定折减系数整瓶输入（根目录，18.8 MB） |
| [geometry_zhuning_10deg/](./geometry_zhuning_10deg/) | 旧几何基准（朱宁 10° 扇区，`Job-2` / `Job-multi` / `WCM_multi`，250 MPa 三案例跑通记录）。**已废弃**，保留作对照；分层为 `scripts/` 生成脚本、`archive/` 历史件、`results/` 验证产物 |
| [wcm_hashin_umat/baseline/](./wcm_hashin_umat/baseline/) | 旧程序在新工程中的**只读基准副本**，由 `source_manifest.json` 做哈希保护，用于追溯改动 |

旧基准的结果（例如 A 版整瓶跑到 250 MPa）不作为当前结论，仅作历史对照。

---

## 二、新的渐进损伤代码及其修改说明

当前唯一在用的实现：WCM 等效层 + 双角度（±θ 两族独立损伤历史）3D Hashin 渐进损伤 UMAT，两条技术路线并行：

| 版本 | 子程序 | PROPS | `*Depvar` | 可见 SDV | 沙漏控制 |
|---|---|---|---|---|---|
| **A · 直接刚度折减** | `hashin_constant_wcm.for` | 25 | 20 | SDV1–8 | 不加（本轮口径） |
| **B · 断裂能线性软化** | `hashin_energy_wcm.for` | 30 | 62 | SDV1–10 | 复材截面 `hourglass=ENHANCED` |

四模式起始判据在 SDV1–4；SDV5/6 为纤维/基体失效（A：标志；B：软化完成），B 的 SDV7/8 为新黏性损伤、SDV9/10 为本帧实际采用损伤。刚度保留 1e-6 残余因子。

| 路径 | 说明 |
|---|---|
| [wcm_hashin_umat/](./wcm_hashin_umat/) | 正式工程：`src/wcm_core.for` + `src/umat_entry.for.in`（模板）→ `dist/*.for`（构建产物，可直接交给 Abaqus 编译） |
| [wcm_hashin_umat/CHANGELOG.md](./wcm_hashin_umat/CHANGELOG.md) | **修改说明主文档**：相对原始单层 UMAT 的逐项有意改动（判据分母改为 S12、参数槽位命名与精简、SDV 精简、基准切换到 rebuild、软化区间试验分支及回退） |
| [wcm_hashin_umat/PROPS说明.md](./wcm_hashin_umat/PROPS说明.md) | 材料参数逐项含义、单位、取值来源 |
| [wcm_hashin_umat/SDV说明.md](./wcm_hashin_umat/SDV说明.md) | 状态变量编号、显示名与判读规则 |
| [wcm_hashin_umat/tools/](./wcm_hashin_umat/tools/) | `prepare_test_cases.py`（由源 deck 生成 A/B 两套待测 INP，带 `--check` 幂等复算）、`extract_failure.py`（积分点阈值后处理）等 |
| [wcm_hashin_umat/config/](./wcm_hashin_umat/config/) | `failure_thresholds_A.json` / `failure_thresholds_B.json` 失效判定阈值 |
| [wcm_hashin_umat/validation/](./wcm_hashin_umat/validation/) · [tests/](./wcm_hashin_umat/tests/) | 离线校验证据（含 `test_cases.json`）与回归测试 |
| [test/](./test/) | **当前运行入口**：`A_stiffness_reduction/`、`B_energy_softening/` 各只保留「UMAT + 配套 INP + README」三件，见 [test/README.md](./test/README.md) |

### B 版候选工程（最新一轮改动都在这里）

| 路径 | 说明 |
|---|---|
| [wcm_hashin_umat_B_candidate/](./wcm_hashin_umat_B_candidate/) | 隔离候选工程：`src/`、`dist/`（当前最新 `hashin_energy_wcm.for`）、`tools/`（build / run_job / 网格研究等 15 个脚本）、`tests/`、`validation/`（18 份证据 JSON）、`results/constitutive_contract.md`（本构契约）、`package/`（三件套诊断候选包，标签：**不可用于可信爆压**） |
| [audit_20260927/B版能量渐进损伤_实施交接规划.md](./audit_20260927/B%E7%89%88%E8%83%BD%E9%87%8F%E6%B8%90%E8%BF%9B%E6%8D%9F%E4%BC%A4_%E5%AE%9E%E6%96%BD%E4%BA%A4%E6%8E%A5%E8%A7%84%E5%88%92.md) | 冻结→隔离→复现核对→起始点定位→材料点/单元验证→网格试验→受控整瓶试跑的路线与授权边界 |
| [audit_20260927/B版能量渐进损伤_首轮实施复盘.md](./audit_20260927/B%E7%89%88%E8%83%BD%E9%87%8F%E6%B8%90%E8%BF%9B%E6%8D%9F%E4%BC%A4_%E9%A6%96%E8%BD%AE%E5%AE%9E%E6%96%BD%E5%A4%8D%E7%9B%98.md) | 首轮结论分层：**实现错误（已修）/ 本构定义问题（记录不擅改）/ 长度范围不相容**，每项注明由哪个试验复核 |
| [audit_20260927/matrix_failure_revision_20260927_215103/修改说明.md](./audit_20260927/matrix_failure_revision_20260927_215103/%E4%BF%AE%E6%94%B9%E8%AF%B4%E6%98%8E.md) | 基体模式无有效软化区间时的局部失败处理 |
| [audit_20260927/fiber_onset_caps_revision_20260927_221516/修改说明.md](./audit_20260927/fiber_onset_caps_revision_20260927_221516/%E4%BF%AE%E6%94%B9%E8%AF%B4%E6%98%8E.md) | 纤维压缩起始量改用应力口径 `delta_fc=LC·max(−S11,0)/P1`，消除 1022 不可播种 |
| [audit_20260927/damage_one/](./audit_20260927/damage_one/) · [props_layout/](./audit_20260927/props_layout/) · [sdv_definition/](./audit_20260927/sdv_definition/) · [instant_fallback/](./audit_20260927/instant_fallback/) | 逐轮取证记录（更新记录、Data Check 摘要、参数槽位排版、`DFINAL≤DEL0` 即时完成试验版） |
| [audit_20260926/](./audit_20260926/) | 早期代码审查：[审查结论.md](./audit_20260926/%E5%AE%A1%E6%9F%A5%E7%BB%93%E8%AE%BA.md)、[WCM封头角度与UMAT接入验证计划.md](./audit_20260926/WCM%E5%B0%81%E5%A4%B4%E8%A7%92%E5%BA%A6%E4%B8%8EUMAT%E6%8E%A5%E5%85%A5%E9%AA%8C%E8%AF%81%E8%AE%A1%E5%88%92.md)、`core_recheck/` 判据复算、`sdv_compact/` SDV 精简记录 |

---

## 三、已做的案例（重点：基于断裂能的渐进损伤）

计算基准为朱宁 rebuild 模型：WCM 单 part、135 个等效 Bin、51,558 个复材单元（C3D8R 50,916 / C3D8 636 / C3D6 6）、391 个复材截面，静压 0→250 MPa，`*Static` 取 0.01, 1, 1e-20, 0.1。

| 案例 | 路径 | 内容与状态 |
|---|---|---|
| 基准几何与未损伤静压 | [geometry_zhuning_joint_rebuild/](./geometry_zhuning_joint_rebuild/) | `WCM_Job.inp`（内置弹性，跑到 t=1.0 成功，作未损伤基准）、`Job-remesh.inp`（加密网格）、`*_WindAngles.ang` 缠绕角数据、`scripts/` 生成脚本、`results/` 验证产物、[建模记录.md](./geometry_zhuning_joint_rebuild/%E5%BB%BA%E6%A8%A1%E8%AE%B0%E5%BD%95.md) |
| 待测入口 A / B | [test/A_stiffness_reduction/](./test/A_stiffness_reduction/) · [test/B_energy_softening/](./test/B_energy_softening/) | `Rebuild_PDM_A.inp` / `Rebuild_PDM_B.inp` + 配套 `.for`，两套除 Bin 本构卡与 B 的沙漏控制外逐字一致 |
| **B 版能量软化整瓶试跑** | [test/B_energy_softening/runs/](./test/B_energy_softening/runs/) | `B_20260927_150554/`（正式版，约 0.170 停止）、`B_20260927_153131_fallback/`（即时完成试验版，0.4428 停止）；保留 `.sta` / `.dat` 收敛历史与提交清单，**均未完成整瓶** |
| **B 版加密网格断裂能对照** | [remesh_energy_tests_20260927/](./remesh_energy_tests_20260927/) | `runs/B_remesh_M21_*`（M21 断裂能组三次提交 214056 / 215350 / 221738）、`runs/B_remesh_Gft300_*`、`runs/B_remesh_span_stab_*`（跨度口径软化 + 稳定化，最新一次）；`latest_run.json` 记录四组 G 值、黏性、刚度上限、SHA-256 与终止状态；`tools/submit_remesh.py` 为提交脚本 |
| 网格相关性与特征长度研究 | [wcm_hashin_umat_B_candidate/runs/](./wcm_hashin_umat_B_candidate/runs/) | `mesh_block_*`（G=40/150 与不同特征长度组合的单块试验）、`mesh_survey_20260927/`（临界 ρ 扫描 `rho_critical.csv`）、`material_onset/matrix/energy/reproduce`、`core_*`、`g4_fixture_*`、`abaqus_elastic/rotation/damage/completion_*` 等材料点与单元级验证运行 |
| 编译与小模型接口验证 | `wcm_hashin_umat/runs/` | 两版 UMAT 的 Abaqus 小模型与整瓶 Data Check 记录（`tank_B_*` 等），以及旧 84 项接口时期的历史运行 |

**当前结论（务必连同状态一起读）**：两版 UMAT 均已编译并通过小模型；整瓶 B 版 Data Check 通过，A 版因 C3D8R 零沙漏刚度失败；**尚无任何一版取得完整整瓶求解结果**，本库不提供"标定后的爆破压力"。能量法遗留的核心障碍是断裂能与单元特征长度的相容性（在首轮复盘中列为长度范围不相容），而非单纯的求解器参数问题。

---

## 四、本库不包含的内容

- **插件工程 `abaqus_plugin_type4/`**（IV 型瓶体参数化建模工具，按要求不上传）。
- Abaqus 求解二进制与大型缓存：`.odb`、`.stt`、`.mdl`、`.sim`、`.res`、`.cax`、`.prt`、`.com`、`.msg`、`.dll`、`.env`、`.cae`、`.jnl`、`.sat`、`.rpy`、`abq.app_cache/`。
- 各轮快照中的重复输入卡：同一 deck 的多份拷贝已按 SHA-256 去重，只保留正式入口那一份。
- 逐单元场数据导出（>5 MB 的 `actual.json` / `strain_field.json`），可用 `tools/` 内脚本重算。
- 纯历史副本目录：`backup/`、`backup_docs/`、`*_regenerated/`、`datacheck_*/`。

明细清单：

| 文件 | 含义 |
|---|---|
| [_excluded/uploaded_files.txt](./_excluded/uploaded_files.txt) | 内容筛选阶段的上传文件清单（按工程原路径逐条列出，不含后续补入的 `.gitattributes`、`原库README.md` 与本目录自身） |
| [_excluded/duplicate_decks.txt](./_excluded/duplicate_decks.txt) | 被去重的输入卡及其保留位置 |
| [_excluded/oversized_derived.txt](./_excluded/oversized_derived.txt) | 被剔除的大体积派生数据 |

---

## 五、复现要点

1. **子程序与材料卡必须配对**：A 用 `hashin_constant_wcm.for` + PROPS 25 项，B 用 `hashin_energy_wcm.for` + PROPS 30 项；只换 `.for` 不换材料卡必然槽位错位。
2. **待测 INP 由同一源 deck 生成**：`python wcm_hashin_umat/tools/prepare_test_cases.py`（加 `--check` 只复算比对、不写文件）。
3. **提交示例（B 版）**：`abaqus job=<名> input=<名>.inp user=hashin_energy_wcm.for cpus=12 interactive`。
4. **后处理**：`python wcm_hashin_umat/tools/extract_failure.py`，阈值读 `wcm_hashin_umat/config/failure_thresholds_{A,B}.json`。
5. **场输出**为 LE、PE、PEEQ、PEMAG、S、U 加编号 SDV；UMAT 不返回非弹性应变，PE/PEEQ/PEMAG 预期恒为 0，不要读作"无塑性"。
6. 旧编号（84/18/60 项接口）的 ODB 与 restart 不可与新接口混用。

---

## 六、文献与出处

参数与判据来源见 [references/文献清单.md](./references/%E6%96%87%E7%8C%AE%E6%B8%85%E5%8D%95.md)，PDF 全文一并放在 `references/`（Vo 2013；Zhang 2024；朱宁 2023 浙江大学硕士；Wu 2024；Lin 2021；张哲 2025 中国石油大学（华东）硕士；以及渐进损伤附录）。这些文献版权归原出版方与作者，本库仅作参数溯源用途。

## 联系与声明

原作者：中国石油大学（华东），Z25150070@s.upc.edu.cn。本库为研究过程中的实验性代码与案例归档，判据与参数仍在修订中，请勿直接用于工程定案。
