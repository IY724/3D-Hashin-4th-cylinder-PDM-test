# B 版 GPT 收敛控制补丁整瓶计算记录

- 创建时间：2026-09-28 19:48:46 +08:00
- 来源：Qoder 克隆目录中 2026-09-28 19:40:51 的第二次作业文件快照。
- INP SHA256：D94AE54824587191DD7D53BEA7838B5425A0749550C558535EFB771D298C3F47
- UMAT SHA256：35CCE3E8A930853B71A4FE06042B20536BA8A4668AA7D33E50E1F3221034577B
- 模型：朱宁重建后的加密网格整瓶；B 能量软化，Gft/Gfc/Gmt/Gmc=133/40/0.6/1.6 N/mm。
- GPT 补丁：UPDATE=1、黏性时间 0.001/0.001、最大黏性损伤跳变 0.02、静力初始增量 0.002/最大增量 0.005。
- 载荷：0→250 MPa；Abaqus/Standard 2025，12 核。
- 最终状态：用户于 2026-09-28 20:09 主动停止；未完成整瓶分析。

## 后续判读

- 以本目录 .sta/.msg/.dat/launcher 文件与实时进程核对状态。
- 未完成作业的最后收敛帧不作爆破压力。

## 提交

- 提交时间：2026-09-28 19:49:01 +08:00
- 命令：abaqus job=B_remesh_convctrl input=B_remesh_convctrl.inp user=hashin_energy_wcm.for cpus=12 scratch=<本目录>/scratch interactive
- 启动进程 PID：39592
- 当前状态：已完成子程序编译与链接、输入处理；求解进行中，尚未完成。


## 实时核对（2026-09-28 19:50:35 +08:00）

- 编译/链接：成功；输入文件处理：成功。
- STA 最后收敛增量：1     5   1     0     1     1  0.0165     0.0165     0.005000
- `standard.exe` 仍在运行；尚未出现完成/失败终态。

## Qoder 在克隆目录的前两次提交

- 19:30 作业：launcher 明确记录 Abaqus analysis interrupted by the user；属于外部中断。
- 19:40 作业：STA 至第 8 增量、步时间 0.0315；进程已退出，但无分析完成/失败尾行，且遗留锁文件。记录为异常停止，现有日志不能归因为本构失效。
- 末次复核（2026-09-28 19:51:01 +08:00）：第 7 增量收敛，步时间 0.0265，standard.exe PID 29300 仍运行。


## 用户停止与切步记录（2026-09-28 20:12 +08:00，按用户 20:xx 的说明修正）

- `.sta` 首次未接受尝试：第 22 增量，步时间 0.0965；`.msg` 明确为 `USER SUBROUTINE REQUESTS A TIME INCREMENT CHANGE`。
- 最后已收敛：第 92 增量，步时间约 0.110，增量约 6.993e-5；截至该点 UMAT 主动切步 71 次，通常下一次尝试只需 1 次平衡迭代。
- 本次补丁子程序相对上一轮只新增“任一角度族任一损伤模式的黏性损伤增量 > 0.02 时设置 PNEWDT<1”的全局切步控制；未记录 `WCM_ERROR` 或 `TOO MANY ATTEMPTS`。
- 用户确认是本人主动停止。20:09 `launcher.txt` 的 Abaqus/Standard exited with an error 是停止后的程序退出记录；`.msg` 在第 93 增量首次尝试开始处结束。
- Windows Application Error 事件 ID 1000 同时记录了 `SMAEqsDirSolverUnsymmetric.exe`/`ucrtbase.dll` 和 `0xc0000409`。在已知用户主动停止的条件下，不将该退出事件解释为求解器自然崩溃或本构错误。
- ODB 为未完成作业的部分结果，不据此判定爆破压力。


