# B版能量渐进损伤候选包（诊断候选）

日期：2026-09-27（阶段4后更新）。**标签：诊断候选，不可用于可信爆压。**

## 内容（仅三件）

| 文件 | 说明 | SHA-256 |
|---|---|---|
| `hashin_energy_wcm.for` | B版UMAT（起始定位修正+FE噪声鲁棒分支+切线不播种） | 48148693270ee5175187020b397136b12766fc6cfaf52c8a38a08ff7e3f35121 |
| `Rebuild_PDM_B.inp` | 朱宁rebuild整瓶输入（51558单元/135 Bin/391 Enhanced截面，62 SDV） | 435e31cf314e767a53770044be5f62dcbc563b42165303105c3cc0fa9d9befc3 |
| `README.md` | 本文件 | — |

## 相对正式版（test/B_energy_softening）的差异

1. `WCM_DAMAGE_B` 起始量定位：由"首次FI≥1的**增量端点**"改为"沿增量应变路径
   EL(s)=B·(E0+s·(E1−E0)) 定位**最早有效起始点**"（FI≥1 且 DEL/SIG 非正者跳过，
   分支进入按边界值处理，不伪造F=1）；新增错误码 41–44；kappa 下限加入 delta0。
2. 噪声鲁棒分支（阶段4发现）：FE舍入使纯剪 Q=S2+S3 带各IP不同的 ±1e-14 噪声，
   旧式符号判分支会被噪声随机化（同一单元各积分点起始量发散）。现按路径噪声尺度
   EPSG=1e-10·max|EFF| 判退化：Q/S1 为数值零时张力支拥有整个路径、压缩支不起始；
   候选点须落在本模式有意义的分支内（GC 过滤）。
3. 切线不再播种（阶段4发现）：WCM_TANGENT_B 的差分调用传 IONSET=0，不执行
   新起始/端点回退/40+I 检查——起始决策只发生在主调用；已起始模式在切线中正常演化。
   否则切线扰动会把 Q 推出噪声带，使 MC 在共享剪切面上虚播种并报 1034。

材料卡、PROPS 30项、62 SDV 布局、A版行为、生产 UPDATE=0 配置均不变。

## 验证状态（详见 ../validation/status.json 与 ../validation/abaqus_candidate.json）

- G0 隔离 PASS；G1 复现 PASS；G2 起始修正 PASS；G3 = LENGTH_LIMITED。
- **G4 PASS（2026-09-27）**：Intel Fortran 2021.7.0 实际编译+Standard求解。
  既有4套件（elastic 4.5e-8 / rotation 9.1e-12 / damage / completion）+
  三种单元能量 fixture（C3D8/C3D8R/C3D6，实测 CELENT=0.1/0.1/0.0794 mm，
  与材料点解析参考 τ/状态差 ≤1e-4（实测 ~8e-8），能量账本 0.4498–0.4521 N/mm
  vs 0.45 目标）全部通过。
- **LENGTH_LIMITED 仍成立**：纯剪基体拉伸路径 Lc_crit≈0.587 mm、MC单轴≈0.158 mm；
  生产气瓶 LC≈2.23–2.34 mm 超域，1033 是固有不相容，非本修复可消除。
- **UNRESOLVED_1022 保留**。
- 阶段5（网格）、阶段6（整瓶）未授权执行。

## 使用约束

1. 不得用于爆压预测或报告"已修复"。
2. 有效长度域内（如 LC≤0.5 mm 的小模型）可用于继续本构验证。
3. 正式 `test/A`、`test/B` 未被本候选修改（SHA-256 见 ../validation/source_manifest.json）。
4. 重启/续算不得混用旧版与新版的 SDV 状态（起始量语义不同）。
