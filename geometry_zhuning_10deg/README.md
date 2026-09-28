# geometry_zhuning_10deg — 三个跑通的成功案例（2026-09-26 整理）

本目录收敛为 **3 个 250 MPa 内压全部跑通（t = 1.0）的案例**；失败尝试与求解中间
文件已在本日整理中移除（见文末记录）。

## 一、总览

| # | 结构 | CAE 文件 | 模型名 | 作业 | 子程序（WCM 生成 UVARM） | 缠绕角数据 | 结果 |
|---|------|----------|--------|------|--------------------------|------------|------|
| 1 | **单 part**（一体件基体 Vessel_Base + 缠绕层 Tank-1，qoder 建模） | `zhuning_base_10deg.cae` | `Zhuning_Base_10deg-Copy` | `Job-2` | `WCM_Zhuning_Base_10deg-Copy_Uvarm.for` | `WCM_Zhuning_Base_10deg-Copy_WindAngles.ang` | ✔ 11 增量到 t=1.0 |
| 2 | **多 part**（PA6_Liner + Al6061_Bosses 分件，qoder 新建） | `zhuning_base_10deg_multipart.cae` | `Zhuning_Base_10deg_Multipart` | `Job-multi` | `WCM_Zhuning_Base_10deg_Multipart_Uvarm.for` | `WCM_Zhuning_Base_10deg_Multipart_WindAngles.ang` | ✔ 15 增量到 t=1.0 |
| 3 | **多 part**（原老模型 zhuning-cylinder-1/2/3，由 G 盘 STEP 重建） | `zhuning_base_10deg_multipart.cae` | `Model-2` | `WCM_multi` | `WCM_Model-2_Uvarm.for` | `WCM_Model-2_WindAngles.ang` | ✔ 15 增量到 t=1.0 |

- 两个 CAE 均保留：`zhuning_base_10deg.cae`（案例1）与 `zhuning_base_10deg_multipart.cae`（案例2、3 两个 Model）。
- 配套 journal 保留：`zhuning_base_10deg.jnl`、`zhuning_base_10deg_multipart.jnl`（含 Job-multi / WCM_multi 提交记录）。
- 共同工况：10° 扇区 + 循环对称（n=36），内压 **250 MPa**，`*Static 0.01, 1., 1e-20`，`nlgeom=YES`。
- 案例 1 提交时曾报 “More than one coordinate transformations … on 56 nodes”，将 BC-2/BC-3 的
  局部坐标系引用统一后提交通过（提交过程记录于 `zhuning_base_10deg.jnl`）。

### 作业文件构成（每个作业保留 5 件）
`*.inp`（输入）、`*.sta`（进度/完成证据）、`*.msg`（求解消息）、`*.log`（提交日志）、`*.odb`（结果）。

### ⚠ 子程序路径警示
三个 `*_Uvarm.for` 内部以 `fullPathName` 硬编码了对应 `*_WindAngles.ang` 的**绝对路径**
（当前均指向本目录）。目录整体搬迁后必须同步修改该路径并重新编译，否则作业启动即失败。

## 二、材料属性（三个案例一致）

| 材料 | 用途 | 密度 (t/mm³) | 弹性 | 塑性（两点硬化） |
|------|------|--------------|------|------------------|
| `PA6`（老模型同义） | 塑料内衬 | 1.13e-9 | E = 1880 MPa, ν = 0.4 | 47 @ 0 + 55 @ 0.17 |
| `Al6061_T6`（老模型为 `al-6061`） | 金属瓶口/阀座 | 2.7e-9 | E = 70000 MPa, ν = 0.32 | 306 @ 0 + 340 @ 0.12 |
| `t700`（老模型为 `t-700`） | T700 单层（供 WCM 等效化，Tank-1 实际用 Bin 卡） | 1.75e-9 | E1 = 141000, E2 = E3 = 11400, ν12=ν13=ν23 = 0.28, G12=G13=G23 = 7100 MPa | — |

失效卡（`t700` 及全部 WCM Bin 卡）：

- `*Fail Stress`（顺序 Xt, Xc, Yt, Yc, S）：`2080., 1250., 60., 290., 110., 0., 0.`
  —— 压缩值按约定填正数，Abaqus 自动取负（提交时提示 “COMPRESSIVE FAILURE VALUES MUST BE NEGATIVE. NEGATIVE OF INPUT TAKEN.” 属正常）。
- `*Fail Strain`：
  - 案例 2、3：`0.015, 0.009, 0.005, 0.025, 0.15`
  - 案例 1：`0.01475, 0.00887, 0.00526, 0.02544, 0.01549`（按 ε = σ/E 精确换算的版本）

WCM 等效层材料（Tank-1）：每模型 **131 个** `WCM_Tank1_Mat1_BinNN` 卡，各自带
`*Elastic, type=ENGINEERING CONSTANTS`（±θ 均匀化等效常数，随缠绕角 β 变化）
+ 上述失效卡 + `*User Output Variables 39`。示例（β = 15° 的 Bin0）：
`116940., 11903.3, 11528.5, 0.829782, 0.112683, 0.262969, 14513.3, 7100., 7100.`

## 三、铺层（Tank-1 缠绕层，三个案例一致）

由 `scripts/extract_tank1_layup.py` 对三个 `*.inp` 的 Tank-1 筒身（|y| < 5 mm）节点半径
实测：

- 内半径 **140.0 mm** → 外半径 **161.6 mm**，总厚 **21.6 mm**，共 **24 个径向实体单元层**；
- 三个厚度带（内→外）：8 × 0.9 mm + 8 × 1.2 mm + 8 × 0.6 mm = 7.2 + 9.6 + 4.8 = 21.6 mm；
- 换算 0.3 mm 名义单层：0.9 = 3 层、1.2 = 4 层、0.6 = 2 层打包 ⇒ 物理层数 8×(3+4+2) = **72 层**，
  与朱宁论文表 3.2（72 层 × 0.3 mm，总厚 21.6 mm）一致；
- 缠绕角主值（`*_WindAngles.ang`，弧度制）：90°、15°、20°、25°、35°、40°；
  其中 90°（π/2）环向单元 23478 个，约占 46%；
- 对照：原案例 `Job-hashinnewcdm.inp` 筒身为 74 层 22.2 mm（历史记录），当前重建按 72 层 21.6 mm。

纵向 131 个 Bin 分区与角度集合对应材料卡的逐 Bin 等效常数；每单元角度由 `.ang` 在
`UEXTERNALDB` 中读入并映射（`NumElemsWithUvar = 51078`）。

## 四、边界条件与载荷（三个案例一致）

- `*Transform, nset=_T-Tank-1-TankCenter, type=C`：圆柱坐标变换，轴为全局 Y 轴；
- BC-2：`YSYMM`（对称面）；BC-3：`U3 = 0`（端面轴向约束）——均经上述圆柱坐标系施加；
- `*Cyclic Symmetry Model, n = 36`（10° 扇区循环对称）；
- `*Dsload`：加载面（案例 1 名 `load`，案例 2/3 名 `Surf-2`）**P = 250.**
- 分析步：`*Step, nlgeom=YES, inc=100000`（案例 3 为 1000000）；
  `*Static 0.01, 1., 1e-20, 0.1`（案例 1 最大增量 1.0）；
- Tie：`Constraint-*`（案例 2/3 为 liner↔boss、Surf-1↔Tank-1.TieSurf；案例 1 为 Surf-1↔Tank-1.TieSurf）。

## 五、复现步骤

1. CAE 中打开 `zhuning_base_10deg_multipart.cae`（列表含 `Zhuning_Base_10deg_Multipart` 与 `Model-2`）
   或 `zhuning_base_10deg.cae`（含 `Zhuning_Base_10deg-Copy`）；
2. 提交对应作业时选择 user subroutine：
   - `Job-2` → `WCM_Zhuning_Base_10deg-Copy_Uvarm.for`
   - `Job-multi` → `WCM_Zhuning_Base_10deg_Multipart_Uvarm.for`
   - `WCM_multi` → `WCM_Model-2_Uvarm.for`
3. 后处理：UVARM 共 39 项（角度、主方向应力、各强度理论比值、首层失效判据入口），
   本目录模型为 WCM 等效层 + `*Fail Stress/Strain`，**不含 Hashin 渐进损伤 UMAT**。

## 六、目录文件清单（整理后）

| 类 | 文件 |
|----|------|
| 模型 | `zhuning_base_10deg.cae` / `.jnl` / `.sat`；`zhuning_base_10deg_multipart.cae` / `.jnl` |
| 案例 1 | `Job-2.{inp,sta,msg,log,odb}` + `WCM_Zhuning_Base_10deg-Copy_{Uvarm.for,WindAngles.ang}` |
| 案例 2 | `Job-multi.{inp,sta,msg,log,odb}` + `WCM_Zhuning_Base_10deg_Multipart_{Uvarm.for,WindAngles.ang}` |
| 案例 3 | `WCM_multi.{inp,sta,msg,log,odb}` + `WCM_Model-2_{Uvarm.for,WindAngles.ang}` |
| 参数/材料 | `case_materials_and_layup.json`、`case_material_cards.inp`、`parameters.json`、`profiles.json`、`wcm_meridian.csv`、`material_update_validation.json`、`latest_run.json` |
| 子目录 | `scripts/`（建模/渲染/铺层提取脚本，含 `extract_tank1_layup.py`）；`results/`（渲染图）；`archive/`（历史 runs 与取证，仅供参考） |

## 七、本次整理移出的内容（记录备查）

- 失败作业（未到达 t=1.0）：`WCM_Job`（t=0.765 退出）、`Job-newpu`（t=0.147）、
  `Job-oldmodel`（t=0.217）、`WCM_old`（未进入迭代）、`old/WCM_J`；
- 失败作业所用子程序 `WCM_Zhuning_Base_10deg_Uvarm.for` / `WCM_Zhuning_Base_10deg_WindAngles.ang`；
- 成功作业的求解中间文件（`*.dat/*.com/*.env/*.ipm/*.prt`）；
- 会话与临时文件（`abaqus.rpy*`、`abaqus.guiState*`、`abq.app_cache`、`temp-zhuning-cylinder-new.sat`）；
- `old/` 目录全部内容（2026-03-27 老模型 inp、失败导入会话 `old.cae/old.jnl`、`WCM_Model-1_*` 等）；
- `scripts/` 内一次性诊断脚本与日志（`_probe_*`、`_verify_*`、`_scan_*`、`probe_failed_run*`、`check_g_odb_frames.py`）；
- 归档至 `archive/`（保留备查，不再位于根目录）：`abaqus_plugins/test.lib`（Al6061_T6 材料库测试文件）
  与 `runs/20260926_143805_multipart/`（案例 2 构建验证记录 `multipart_validation.json`）。
