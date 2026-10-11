# 多频相位测距的门限：数学专文

*本页只讲数学：符号、界、拐点、格点模型、分账、频点设计的多目标形式。叙事版见 [README.md](README.md)，推导与实验全文见 [deep.zh.md](deep.zh.md)。*

**数字来源标记。**
〔N〕＝[`numbers.json`](numbers.json)（由 `computation/article_numbers.py` 从归档算出）；
〔V2:eXX〕＝`computation/data/v2/eXX_*.json`；
〔VER〕＝本文新增的核验脚本 [`computation/verify_math_idents.py`](../computation/verify_math_idents.py) 的输出：它只读归档与库函数，不写任何文件；运行 `python3 computation/verify_math_idents.py`，全部 PASS 约 45 s（4 进程）。
本文用到的新库函数：`mfpr/bounds.py` 的 `wwb`（2.5 节）、`mfpr/multivalley.py`（4.6 节）、`mfpr/design.py`（格点/非格点谷、近区、O1 上界，第 6 节）；回归测试在 `computation/tests/test_models.py`（`pytest computation/tests -q`）。
与已有文字不一致之处在文中用 **〔不一致〕** 标出；已在 deep 中改正的标 **〔已修〕**。

---

## 0 符号与观测

**观测。** 真实距离 $d$，$N$ 个频点 $f_1<\dots<f_N$，每个频点只给出包裹相位

$$
\varphi_i=\mathrm{wrap}\!\left(\frac{2\pi f_i d}{c}+e_i\right),\qquad e_i\stackrel{\text{iid}}{\sim}\mathcal N(0,\sigma^2),\qquad \sigma^2\simeq\frac{1}{2\,\mathrm{SNR}} .
$$

**频率矩与派生量。**

| 记号 | 定义 | 含义 |
| --- | --- | --- |
| $\bar f$ | $\frac1N\sum f_i$ | 平均频率 |
| $F$ | $\sum f_i^2$ | 载波（细尺）信息 |
| $S$ | $\sum (f_i-\bar f)^2$ | 包络（粗尺）信息 |
| $R$ | $F/S = 1+N\bar f^2/S$ | 两把尺子的精度比 |
| $\lambda$ | $c/\bar f$ | 平均波长（一次滑周的步长） |
| $B$ | $f_N-f_1$ | 占用带宽 |
| $\Delta_{\max}$ | $\max_i (f_{i+1}-f_i)$ | 最大相邻间隔 |
| $\sigma_{\mathrm{jump}}$ | $\pi\sqrt S/\bar f$ | 门限的相位尺度 |
| $u$ | $\sigma/\sigma_{\mathrm{jump}}$ | 归一化噪声 |
| $x$ | $\sigma_{\mathrm{jump}}/\sigma = 1/u$ | 脚方程的变量 |

**无噪声相位代价。** 距离偏差 $\delta$ 处

$$
J(\delta)=\sum_{i=1}^N \mathrm{wrap}\!\left(\frac{2\pi f_i\delta}{c}\right)^2 .
$$

它的谷就是似然的峰。两个参照：

- **一波长代价** $J(\lambda)=\sum_i \bigl(2\pi(f_i-\bar f)/\bar f\bigr)^2 = 4\pi^2S/\bar f^2 = 4\sigma_{\mathrm{jump}}^2$。只要 $\max|f_i-\bar f|<\bar f/2$，这对任何排布都精确成立〔VER〕。
- **两尺度代价** $J_2(\delta)=N\,\mathrm{wrap}(2\pi\delta/\lambda)^2+(2\pi/c)^2S\delta^2$，即载波梳子加包络抛物线。

**主例。** 等间距 2400:10:2470 MHz，$N=8$，$S=4200\ \mathrm{MHz}^2$，$R=40.5288$ dB，$\sigma_{\mathrm{jump}}=0.0836133$ rad，$\lambda=0.123118$ m〔N〕。贝叶斯比较用长度 $T=14.9896$ m 的均匀先验〔N〕，恰好 $T=c/(2\cdot 10\ \mathrm{MHz})$。

**Q 函数。** $Q(t)=1-\Phi(t)$，其中 $\Phi$ 是标准正态分布函数。

---

## 1 两个局部界：载频直线与带宽直线

### 1.1 两条直线

$d$ 的 Fisher 信息是 $\sum_i g_i^2/\sigma^2$，其中 $g_i=2\pi f_i/c$。若整周已知（相位不包裹），克拉美–罗界为

$$
\mathrm{HMSE}=\frac{c^2\sigma^2}{4\pi^2F}\qquad(\text{载频直线}).
$$

若把所有频点共有的相位 $\theta$ 当作冗余参数（$\varphi_i=g_id+\theta+e_i$），投影掉 $\theta$ 之后剩下的信息是 $\sum(g_i-\bar g)^2/\sigma^2=(2\pi/c)^2S/\sigma^2$，于是

$$
\mathrm{MMSE}=\frac{c^2\sigma^2}{4\pi^2S}\qquad(\text{带宽直线}).
$$

*物理意义：* HMSE 是“细尺 + 正确整周”的精度；MMSE 是“只用相位差（包络）”的精度，也就是非相干估计的精度。两者都是**局部**界，都不随整周结构变化，所以都不跳。

### 1.2 R

$$
\frac{\mathrm{MMSE}}{\mathrm{HMSE}}=R=\frac FS=1+\frac{N\bar f^2}{S},\qquad
\text{等间距：}\ R=1+\frac{12(N-1)}{N+1}\Bigl(\frac{\bar f}{B}\Bigr)^2 .
$$

*物理意义：* $R$ 是两条平行直线在双对数图上的高度差，主例中为 40.53 dB〔N〕。它约等于相对带宽 $B/\bar f$ 平方的倒数。

在固定频段内，$F=N\bar f^2+S\approx N\bar f^2$ 几乎不随排布变化：本文表 6.6 中所有设计的 HMSE 与等间距相差都不到 0.05 dB〔VER〕。**排布不改变载频直线的高度，只改变何时离开它。**

### 1.3 两条直线之间的恒等式

用格点模型的参数（第 4 节；$s$ 是整数判决的标准差，$\lambda_{\mathrm{eff}}=\lambda(1-1/R)$）：

$$
\mathrm{HMSE}+\lambda_{\mathrm{eff}}^2 s^2=\mathrm{MMSE}\quad(\text{精确}),\qquad
\sqrt{\mathrm{MMSE}(\sigma_{\mathrm{jump}})}=\lambda/2 .
$$

*物理意义：* 粗尺的误差可以精确拆成两部分：细尺误差，加上“把整数当连续量估计”时的误差。带宽直线就是**不取整**的格点估计；第 3 节的一切都来自“取整”这一步〔VER〕。

---

## 2 全局界：半窗 Ziv–Zakai 与逐谷分解

### 2.1 定义

先验取 $d\sim U[0,T]$。对等间距网格，$J$ 有精确周期 $P=c/\gcd(f_i)=29.98$ m，并且关于 $P/2$ 对称：$J(P-\delta)=J(\delta)$。所以 $T=P/2$ 已经覆盖了所有互不相同的偏差，这就是“半窗”。

把真值 $d$ 与 $d+h$ 当作两个等概假设。若相位不包裹（高斯相位向量），最小错误概率为

$$
P_e(h)=Q\!\left(\frac{\sqrt{J(h)}}{2\sigma}\right).
$$

谷填充的 Ziv–Zakai 下界为

$$
\mathrm{ZZB}=\frac1T\int_0^T h\;\mathcal V\!\left[(T-h)P_e(h)\right]dh,\qquad \mathcal V g(h)=\max_{v\ge h}g(v).
$$

实现见 `mfpr/bounds.py`。计算时每个谷底都单独精细化；否则在低噪声端会错 30 dB 以上（deep 附录 C）。

实际观测是包裹相位，上式用的却是高斯 $P_e$。下面的引理说明这样做合法。

> **引理 2.1（包裹只会让二元检验更难）。** 对包裹观测 $\varphi=\mathrm{wrap}(gd+e)$，区分 $d$ 与 $d+h$（等先验）的最小错误概率满足
> $$P_e^{\mathrm{wrap}}(h)\ \ge\ Q\!\left(\frac{\lVert\mathrm{wrap}(gh)\rVert}{2\sigma}\right)=Q\!\left(\frac{\sqrt{J(h)}}{2\sigma}\right).$$
>
> *证明。* 取整数向量 $m$ 使 $r\equiv gh-2\pi m=\mathrm{wrap}(gh)$。由于加 $2\pi m$ 不改变包裹结果，$H_1$ 下 $\varphi=\mathrm{wrap}(gd+r+e)$。所以包裹问题是高斯问题 $y\sim\mathcal N(gd,\sigma^2I)$ 对 $\mathcal N(gd+r,\sigma^2I)$ 经过同一个确定映射 $\mathrm{wrap}(\cdot)$ 之后的问题。任何基于 $\varphi$ 的检验都是基于 $y$ 的检验（先包裹再判），所以它的错误概率不低于高斯最优值 $Q(\lVert r\rVert/2\sigma)$。∎

**推论。** 在均匀先验下 $P_e^{\mathrm{wrap}}$ 与 $d$ 无关（平移 $d$ 是环面上的旋转）。ZZB 的被积式对 $P_e$ 单调不减（$\mathcal V$ 和乘以 $T-h\ge0$ 都保序），所以把真实的 $P_e^{\mathrm{wrap}}$ 换成 $Q(\sqrt J/2\sigma)$ 只会让界变低。**上式对包裹模型是合法的下界，代价是可能偏松。**

数值〔VER〕（每点 $2\times10^5$ 或 $4\times10^4$ 次似然比检验）：

| $\sigma$ | $h$ | $\max_i\lvert r_i\rvert$ | $P_e^{\mathrm{wrap}}$ | $Q(\sqrt J/2\sigma)$ | 比 |
| --- | --- | --- | --- | --- | --- |
| $\sigma_{\mathrm{jump}}=0.0836$ | 0.123 m（$1\lambda$ 谷） | 0.09 rad | 0.1591 ± 0.0006 | 0.1587 | 1.003 |
| $10\sigma_{\mathrm{jump}}=0.836$ | 5.602 m | 2.56 rad | 0.0071 ± 0.0003 | 0.0047 | 1.5 |

*物理意义：* 只有当某个频点的残差 $r_i$ 靠近 $\pm\pi$、噪声又大到能把相位绕过 $\pm\pi$ 时，包裹才让两个假设更难分（两条路径会落到同一个观测），真实 $P_e$ 才明显高于高斯值，界才偏松。第一道门限上起作用的格点谷残差都很小，所以那里这一步不松。

### 2.2 高信噪比极限

主瓣上 $J(h)\approx(2\pi/c)^2Fh^2$。令 $a=\pi\sqrt F/(c\sigma)$，利用 $\int_0^\infty tQ(t)\,dt=\tfrac14$：

$$
\mathrm{ZZB}\to\int_0^\infty hQ(ah)\,dh=\frac1{4a^2}=\mathrm{HMSE}.
$$

*物理意义：* 主瓣内这个全局界退化成局部界，数值上差 0.00 dB（deep §6）〔VER〕。

### 2.3 逐谷分解

“谷填充”把每一段 $h$ 的积分高度交给它右侧第一个记录谷。按这个“所有权”切分积分，就得到每个谷对下界的贡献。

- **等间距，$\sigma=0.0836$ rad**：$1\lambda$ 占 68.4%，$2\lambda$ 占 28.7%，$3\lambda$ 占 2.8%，$4\lambda$ 占 0.1%。累计 50% 和 90% 分别在 0.106 m 和 0.224 m 处达到〔N〕〔VER〕。
- **几何间隔（公比 1.02）**：6.211 m 处的远谷占 96%〔N〕。

*物理意义：* 全局界“知道”误差是从哪个谷来的。门限上，等间距的误差几乎全部来自前两个整波长谷。

### 2.4 格点闭式（新）

在格点主导的设计里，谷底位于 $k\lambda_{\mathrm{eff}}$，深度为 $J_k=k^2J(\lambda)(1-1/R)$，于是 $\sqrt{J_k}/(2\sigma)=k/(2s)$。在区间 $((k-1)\lambda_{\mathrm{eff}},k\lambda_{\mathrm{eff}}]$ 上，填充高度近似为 $(T-k\lambda_{\mathrm{eff}})P_{e,k}\approx T P_{e,k}$，所以

$$
\boxed{\ \mathrm{ZZB}\approx\mathrm{HMSE}+\frac{\lambda_{\mathrm{eff}}^2}{2}\sum_{k\ge1}(2k-1)\,Q\!\left(\frac{k}{2s}\right)\ }
$$

在 $u\in[0.15,3]$ 上，它与数值 ZZB 的差不超过 0.09 dB〔VER〕。把它与相干极大似然的格点式并排放：

$$
\mathrm{MSE}_{\mathrm{ML}}=\mathrm{HMSE}+\lambda_{\mathrm{eff}}^2\sum_{k\ge1}(2k-1)\cdot 2Q\!\left(\frac{k-\frac12}{s}\right).
$$

两式结构相同，差别只有两处：

1. 每一项的系数，ZZB 是 $\tfrac12$，ML 是 $2$；
2. Q 的参数，ZZB 是“两假设的半距离” $k/2s$，ML 是“越过判决门槛” $(k-\tfrac12)/s$。

**低信噪比渐近。** 由 Euler–Maclaurin 公式，$\sum_{k\ge1}(2k-1)Q(k/2s)=2s^2-2s/\sqrt{2\pi}+\tfrac16+O(1/s)$〔VER〕，于是

$$
\frac{\mathrm{ZZB}}{\mathrm{MMSE}}\approx 1-\frac{2}{\sqrt{2\pi}\,u}\qquad(u\gg1,\ \text{第二道门限之前}).
$$

*物理意义：* 下界本身也有“两条直线”的结构，只是它的过渡被抹宽了，按 $1/u$ 的速度从下方回到带宽直线。$u=5$ 时 $0.853$（渐近式给 $0.84$），$u=10$ 时 $0.923$（渐近式给 $0.920$）〔VER〕。

### 2.5 单测试点 Weiss–Weinstein 界（新）

**定义。** 先验仍为 $d\sim U[0,T]$。对测试点 $h>0$ 和指数 $\nu\in(0,1)$，令

$$
\eta(\nu,h)=\ln\frac{T-h}{T}+\ln C_\nu(h),\qquad C_\nu(h)=\int p(\varphi\mid d+h)^{\nu}\,p(\varphi\mid d)^{1-\nu}\,d\varphi .
$$

$C_\nu$ 只依赖 $r=\mathrm{wrap}(gh)$，与 $d$ 无关。单测试点 WWB 为

$$
\mathrm{WWB}=\max_{h,\nu}\ \frac{h^2\,e^{2\eta(\nu,h)}}{e^{\eta(2\nu,h)}+e^{\eta(2-2\nu,-h)}-2e^{\eta(\nu,2h)}} .
$$

高斯相位向量下 $\ln C_\nu=-\nu(1-\nu)J(h)/2\sigma^2$。包裹模型下 $C_\nu$ 是逐频点包裹正态 Chernoff 系数之积（周期中点公式，谱收敛）。WWB 对任何模型都成立，只要 $\eta$ 用该模型的真值，所以 `bounds.wwb(..., wrapped=True)` 用高斯 $\eta$ 选出最好的 48 个测试点，再用精确的包裹 $\eta$ 重算并取最大。在 $u=1,3,5$ 上，包裹值与高斯值之差不超过 $2\times10^{-10}$ dB，最优测试点是格点谷（$u=1$ 时恰为 $\lambda_{\mathrm{eff}}$）；在 $u=0.3,1,3$ 上，$\nu\in[0.05,0.95]$ 的最优值都是 $\nu=\tfrac12$〔VER〕。

**格点闭式。** 取 $h=k\lambda_{\mathrm{eff}}$、$\nu=\tfrac12$。在近区内 $J(k\lambda_{\mathrm{eff}})=k^2J(\lambda_{\mathrm{eff}})$，$J(2k\lambda_{\mathrm{eff}})=4k^2J(\lambda_{\mathrm{eff}})$，而 $J(\lambda_{\mathrm{eff}})/\sigma^2=1/s^2$。记 $t=k^2/4s^2$，$a=(T-h)/T$，$b=(T-2h)/T$：

$$
\mathrm{WWB}_k=\frac{h^2a^2e^{-t}}{2a-2b\,e^{-2t}}\ \xrightarrow{\ T\to\infty\ }\ \boxed{\ \lambda_{\mathrm{eff}}^2s^2\,\frac{t}{\sinh t}=(\mathrm{MMSE}-\mathrm{HMSE})\,\frac{t}{\sinh t}\ }
$$

在 $u=1$ 上它与数值 WWB 一致到 $10^{-13}$ dB〔VER〕。三个区段：

| 区段 | 形式 | 读法 |
| --- | --- | --- |
| $u\approx1$（$k=1$，$t\approx1$） | $(\mathrm{MMSE}-\mathrm{HMSE})/\sinh1$ | 比带宽直线低 0.70 dB；有限 $T$ 时 −0.743 dB〔VER〕 |
| $u\gg1$（$t\to0$） | $\to\mathrm{MMSE}-\mathrm{HMSE}$ | 回到带宽直线；有限 $T$ 的 $a,b$ 使它低 0.14–0.28 dB（$u=2$–5） |
| $u\ll1$（脚） | $\approx2\lambda_{\mathrm{eff}}^2s^2t\,e^{-t}\propto e^{-1/4s^2}$ | 指数是 ZZB 滑周项 $Q(1/2s)\propto e^{-1/8s^2}$ 的 2 倍，在脚附近没用 |

**与 ZZB、格点 PM/ML 对比**（dB，相对 MMSE；ZZB、WWB 为数值界，PM、ML 为格点模型）〔VER〕：

| $u$ | 0.2 | 0.3 | 0.5 | 0.7 | 1.0 | 1.5 | 2.0 | 3.0 | 5.0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ZZB | −39.88 | −20.17 | −7.37 | −4.68 | −3.29 | −2.26 | −1.74 | −1.23 | −0.87 |
| WWB | −40.53 | −34.82 | −8.37 | −2.72 | **−0.74** | −0.20 | −0.14 | −0.19 | −0.28 |
| PM | −38.75 | −15.30 | −2.61 | −0.35 | −0.00 | −0.00 | −0.00 | −0.00 | −0.00 |
| ML | −38.36 | −14.17 | −1.38 | +0.97 | +1.14 | +0.60 | +0.35 | +0.16 | +0.06 |
| PM − max(ZZB, WWB) | 1.12 | 4.87 | 4.76 | 2.37 | **0.74** | 0.20 | 0.14 | 0.19 | 0.28 |

WWB 在 $u>0.542$ 时超过 ZZB〔VER〕；每个 $u$ 上 WWB ≤ PM〔VER〕。

*物理意义：* ZZB 只做二元检验，滑周的“代价”只积到 $\int_0^\lambda h\,dh=\lambda^2/2$（5.5 节），在门限上因此偏松 3.3 dB。WWB 把 $h=\lambda_{\mathrm{eff}}$ 这个滑周距离直接放进分子 $h^2$，在门限及其上方几乎贴住带宽直线。代价在分子的 $e^{2\eta}$：滑周稀少时它是 Bhattacharyya 系数的**平方** $e^{-J/4\sigma^2}$，而 ZZB 的 $Q(\sqrt J/2\sigma)$ 只衰减 $e^{-J/8\sigma^2}$。所以在脚上单测试点 WWB 比 ZZB 低几个数量级。两种界各管一段；$u\approx0.3$–0.5 处 PM 与两者中较大者仍差约 4.8 dB，这个缺口两者都没填上（第 7 节）。

---

## 3 拐点：脚、中心、过冲

格点模型（第 4 节）给出 $\hat k=\mathrm{round}(z)$，$z\sim\mathcal N(0,s^2)$，$s=u/(2\sqrt{1-1/R})$，并且

$$
\frac{\mathrm{MSE}}{\mathrm{MMSE}}=\frac1R+\Bigl(1-\frac1R\Bigr)\frac{\mathbb E[\hat k^2]}{s^2},\qquad
\mathbb E[\hat k^2]=\sum_{k\ge1}(2k-1)\,2Q\!\left(\frac{k-\frac12}{s}\right).
$$

### 3.1 中心（$u=1$）

$$
\mathrm{SNR}_{\mathrm{tot}}=\frac{N}{2\sigma_{\mathrm{jump}}^2}=\frac{R-1}{2\pi^2},\qquad
10\log_{10}\mathrm{SNR}_{\mathrm{tot}}=R_{\mathrm{dB}}-12.95+10\log_{10}(1-1/R).
$$

主例中每频点 18.5442 dB，总 27.5751 dB〔N〕；$10\log_{10}(2\pi^2)=12.9533$〔VER〕。

*物理意义：* 中心是粗尺的均方根误差恰好等于 $\lambda/2$ 的地方（$s=\tfrac12$）。此时滑周概率为 $2Q(1)=0.317$〔VER〕。

### 3.2 过冲

在 $s\approx\tfrac12$ 时，

$$
\frac{\mathrm{MSE}}{\mathrm{MMSE}}\approx 4\,\mathbb E[\hat k^2]=8Q(1)+24Q(3)+40Q(5)+\dots=1.302\ \Rightarrow\ +1.145\ \mathrm{dB}.
$$

逐项计算与精确值 1.14481 dB〔N〕的差不到 $10^{-3}$ dB〔VER〕。

*物理意义：* 带宽直线对应 $\mathbb E[z^2]=s^2=\tfrac14$；ML 付出的是 $\mathbb E[\mathrm{round}(z)^2]=0.325$。$s\approx\tfrac12$ 时取整把二阶矩放大约 30%。每次滑周的代价 $\lambda^2$ 恰好是 4 倍 MMSE。过冲与 $R$ 无关，只差 $O(1/R)$。

### 3.3 脚（离开载频直线）

**定义：** 总均方误差等于 2 倍 HMSE，即“滑周贡献 = 细尺方差”。滑周稀少时 $\mathbb E[\hat k^2]\approx 2Q(1/2s)\approx 2Q(x)$。再用 $\lambda^2/\mathrm{HMSE}=4x^2R$，得到

$$
\boxed{\ 8x^2Q(x)=\frac1R\ },\qquad \sigma_{\mathrm{leave}}=\sigma_{\mathrm{jump}}/x .
$$

主例：$x=4.90828$，每频点 32.3628 dB〔N〕。在这一点上，完整格点模型给出 $\mathrm{MSE}/\mathrm{HMSE}=2.0009$，近似误差可以忽略〔VER〕。

**渐近。** 用 $Q(x)\approx\varphi(x)/x$（$\varphi$ 为标准正态密度）得 $x^2\approx 2\ln\bigl(8Rx/\sqrt{2\pi}\bigr)$，所以 $x\sim\sqrt{2\ln R}$，只对数地依赖 $R$。脚处的总信噪比为

$$
10\log_{10}\mathrm{SNR}^{\mathrm{foot}}_{\mathrm{tot}}=R_{\mathrm{dB}}-12.95+20\log_{10}x(R)+10\log_{10}(1-1/R).
$$

| $R_{\mathrm{dB}}$ | 20 | 26 | 30 | 35 | 40.53 | 45 | 50 | 60 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $x(R)$ | 3.748 | 4.126 | 4.357 | 4.628 | 4.908 | 5.123 | 5.351 | 5.780 |
| 脚 − $R_{\mathrm{dB}}$（总 SNR，dB） | −1.52 | −0.65 | −0.17 | +0.35 | +0.87 | +1.24 | +1.62 | +2.28 |
| ML 脚 − ZZB 脚（dB） | 1.06 | 0.84 | 0.74 | 0.64 | 0.56 | 0.51 | 0.46 | 0.39 |

〔VER〕

**〔已修〕** deep §4.2 和 §11 的“安全规则：总 SNR $\gtrsim R_{\mathrm{dB}}+1$ dB”只在 $R\approx40$–45 dB 时成立。一般情形应当用上表，或直接用上式 $R_{\mathrm{dB}}-12.95+20\log_{10}x(R)$；deep 与 README 现已加注指向本表。

*物理意义：* 脚就是“一次滑周的概率 × 代价 $\lambda^2$”追上细尺方差的地方。由于 $\lambda^2/\mathrm{HMSE}\propto R$，而 Q 的尾部极陡，$R$ 只能通过对数影响 $x$。从脚到中心是 $20\log_{10}x=13.8$ dB，而两条直线高差 40.5 dB，所以曲线看起来几乎垂直。

### 3.4 下界的脚：位置被“信息”钉住

由 2.4 节，ZZB 的滑周项在稀少区约为 $\tfrac12\lambda^2Q(x)$。令它等于 HMSE：

$$
2x_Z^2Q(x_Z)=\frac1R\qquad(\text{主例 } x_Z=4.6022).
$$

**推论（在格点 ZZB 近似下）。** 任何估计量都满足 $\mathrm{MSE}\ge\mathrm{ZZB}$，所以任何估计量的脚（$\mathrm{MSE}=2\,\mathrm{HMSE}$）都不可能出现在比 $\sigma_{\mathrm{jump}}/x_Z$ 更大的噪声上。相干 ML 的脚只比这个极限早 $20\log_{10}(x/x_Z)=0.56$ dB（主例），在 $R$ 为 20–60 dB 时是 0.39–1.06 dB〔VER〕。

*物理意义：* 在脚附近，ML 与下界的滑周项差 4 倍（6 dB），可是 Q 太陡，这 6 dB 只对应不到 1 dB 的信噪比。**拐点的位置是信息决定的，过渡段的高度才是估计量和界的松紧决定的。**

### 3.5 平台与第二道门限

$s\gg1$ 时 $\mathbb E[\hat k^2]\to s^2+\tfrac1{12}$（Sheppard 修正）〔VER〕，于是

$$
\mathrm{MSE}\to\mathrm{MMSE}+\lambda_{\mathrm{eff}}^2/12 .
$$

*物理意义：* 取整等于给粗尺再加一份宽 $\lambda$ 的均匀量化噪声。

平台一直持续到第二道门限：单个相位开始被噪声折叠，粗尺自身也失效（Rife–Boorstyn）。坍缩实验中，第二道门限出现在总 SNR 8.0–11.8 dB〔N〕。在主例上，格点模型在 $u\approx6.8$ 处偏 1.09 dB，在 $u=8$ 处偏 5.48 dB〔N〕。

---

## 4 格点模型：把拐点接到整条曲线

### 4.1 推导要点（deep 附录 A）

记 $a_i=2\pi(f_i-\bar f)/\bar f$，即“错一个平均波长”时第 $i$ 个频点剩下的相位。在第二道门限之前，$\delta=k\lambda+\epsilon$ 处的代价为 $\lVert e-ka-\epsilon g\rVert_W^2$。对 $\epsilon$ 求最小，再把 $a$ 投影到 $g$ 的 $W$-正交补 $\tilde a$ 上，得到

$$
J(k)=\text{常数}+\tilde a^\top W\tilde a\,(k-z)^2,\qquad z=\frac{\tilde a^\top We}{\tilde a^\top W\tilde a}\sim\mathcal N(0,s^2),\qquad
s=\frac{\sigma}{2\sigma_{\mathrm{jump}}\sqrt{1-1/R}} .
$$

细误差与 $z$ 不相关（因为 $\tilde a\perp_W g$），所以

$$
\text{误差}=\hat k\lambda_{\mathrm{eff}}+(\text{独立的细误差，方差 HMSE}),\qquad \hat k=\mathrm{round}(z+\beta).
$$

$\beta=\tilde a^\top Wm/\tilde a^\top W\tilde a$ 是确定性相位偏差 $m$ 在“整周方向”上的投影。

### 4.2 一条公式，三个区段

$$
G_R(u)\equiv\frac{\mathrm{MSE}}{\mathrm{MMSE}}=\frac1R+\frac{4}{u^2}\Bigl(1-\frac1R\Bigr)^2\sum_k k^2\Bigl[\Phi\Bigl(\tfrac{k+\frac12}{s}\Bigr)-\Phi\Bigl(\tfrac{k-\frac12}{s}\Bigr)\Bigr].
$$

| 区段 | 极限形式 | 物理意义 |
| --- | --- | --- |
| $u\ll1$ | $\mathrm{HMSE}\,[1+8x^2Q(x)R]$ | 载频直线加稀少滑周 |
| $u\approx1$ | $\mathrm{MMSE}\cdot 10^{0.1145}$ | 中心与过冲 |
| $u\gg1$ | $\mathrm{MMSE}+\lambda_{\mathrm{eff}}^2/12$ | 带宽直线加量化 |

除了 $1/R$ 地板，$G_R$ 与设计无关，这就是“坍缩”：8 组等间距设计与各自模型的偏差，全窗不超过 0.94 dB，局部窗不超过 0.41 dB〔N〕。主例与蒙特卡洛的偏差：全窗 $u\le6$ 不超过 0.18 dB，局部窗不超过 0.04 dB，一直验证到滑周概率 $5\times10^{-5}$〔N〕。

### 4.3 偏差

公共相位 $\varepsilon$ 给出 $\beta=-\varepsilon/2\pi$。近侧的滑周概率变成 $Q((\tfrac12-|\beta|)/s)$，门限因此移动

$$
-20\log_{10}(1-2|\beta|)\ \mathrm{dB}.
$$

$|\varepsilon|=19.57^\circ$ 时移动 1 dB〔N〕。群时延 $m=gc\tau$ 给出 $\beta=0$：距离平移 $c\tau$，门限不动。

### 4.4 多层

两簇设计在第二层上再套用一次同一个模型：细尺换成全带宽，粗尺换成簇内带宽。主例的两簇设计中 $R_2=S/S_{\mathrm{within}}=15.89$ dB，$\sigma_{\mathrm{jump},2}=0.7225$ rad，合成波长 5.45 m〔N〕。

### 4.5 格点模型何时成立

条件是：所有起作用的谷都在 $k\lambda_{\mathrm{eff}}$ 上。第 6.1 节的定理给出一个**充分**条件：$\delta<c/(2\Delta_{\max})$。等间距时这个条件覆盖整个先验窗。条件不满足时（非格点谷起作用），用 4.6 节的多谷模型。

### 4.6 多谷线性化模型（新；取代谷混合作为主预测）

**推导。** 记 $J$ 在 $(\lambda/4,T)$ 内的谷底为 $\delta_v$，两侧都算（$\pm\delta_v$）。谷底处的无噪声残差是固定向量 $A_v=\mathrm{wrap}(g\delta_v)$。在第二道门限之前，在每个谷底把包裹代价线性化：$\delta=\delta_v+\epsilon$ 处代价 $\approx\lVert e-A_v-\epsilon g\rVert^2$。对 $\epsilon$ 求最小，就是把 $e-A_v$ 投影到 $g$ 的正交补上。真齿（$v=0$）的代价是 $\lVert P_\perp e\rVert^2$，所以谷 $v$ 赢过真齿当且仅当

$$
\Delta J_v=\lVert\tilde a_v\rVert^2-2\,\tilde a_v^\top e<0,\qquad \tilde a_v=A_v-\frac{g^\top A_v}{g^\top g}\,g .
$$

所有竞争者都是同一个高斯向量 $e$ 的线性统计量。ML 选 $\hat v=\arg\min_v\Delta J_v$（$\Delta J_0=0$），误差为 $D_{\hat v}+g^\top e/g^\top g$，其中 $D_v=\delta_v-g^\top A_v/g^\top g$。细误差 $g^\top e/g^\top g\sim\mathcal N(0,\mathrm{HMSE})$ 只依赖 $e$ 沿 $g$ 的分量，而 $\hat v$ 只依赖 $P_\perp e$，两者精确独立，所以

$$
\boxed{\ \mathrm{MSE}=\mathrm{HMSE}+\sum_v P(\hat v=v)\,D_v^2\ }
$$

$P(\hat v=v)$ 是高斯多面体概率。谷底处 $\partial J/\partial\delta=2g^\top A_v=0$，所以 $D_v\approx\delta_v$，$\lVert\tilde a_v\rVert^2\approx J_v$。

**三个性质。**

1. **共线时退化为格点模型。** 若所有 $\tilde a_v=k_v\tilde a$（等间距在近区内），则 $\Delta J_k=\lVert\tilde a\rVert^2[(k-z)^2-z^2]$，$z=\tilde a^\top e/\lVert\tilde a\rVert^2$，赢家就是 $\mathrm{round}(z)$（4.1 节）。数值：等间距上在 $u\in\{0.3,0.5,1,2\}$ 与格点模型差不超过 0.025 dB〔VER〕。
2. **滑周稀少时退化为谷和。** 此时 $P(\hat v=v)\approx P(\Delta J_v<0)=Q(\sqrt{J_v}/2\sigma)$，回到 $\mathrm{HMSE}+\sum_v2Q(\sqrt{J_v}/2\sigma)\delta_v^2$。在 $u=0.25$ 上相对差 0.6%〔VER〕。所以脚 $\sigma_{\mathrm{leave}}$ 的定义（第 6 节）不受影响。
3. **中间区段是联合竞争。** 谷 $v$ 与 $w$ 的相关系数是 $\tilde a_v^\top\tilde a_w/(\lVert\tilde a_v\rVert\lVert\tilde a_w\rVert)$。谷混合 $\sum_vp_vD_v^2/(1+\sum p)$ 用的是成对概率 $p_v=2Q(\sqrt{J_v}/2\sigma)$：它只问 $v$ 是否赢过真齿，不问 $v$ 是否也赢过其他谷。远谷的 $D_v^2$ 大，这部分被高估，所以谷混合在过渡区偏悲观。

**计算。** `multivalley.mse(f, σ)`。先剪掉 $Q(\lVert\tilde a_v\rVert/2\sigma)\max D^2<10^{-6}\,\mathrm{HMSE}$ 的竞争者，然后在 $\mathbb R^N$ 中做混合重要性抽样：以 0.2 的概率抽 $\mathcal N(0,\sigma^2I)$，其余按权重 $\propto Q(\lVert\tilde a_v\rVert/2\sigma)D_v^2$ 抽以边界点 $\tilde a_v/2$ 为中心的 $\mathcal N(\tilde a_v/2,\sigma^2I)$。$10^5$ 个样本时，$\sigma=0.2$ 上滑周项的相对标准误为 0.2–0.6%。模型给出的单次平方误差满足 $\mathbb E[X^2]=\sum P_vD_v^4+6\,\mathrm{HMSE}\sum P_vD_v^2+3\,\mathrm{HMSE}^2$，由此得到 $n$ 次蒙特卡洛样本均值围绕预测的 95% 带（`mc_band_db`）。

**适用范围。** 线性化要求竞争谷底处每个频点的残差离 $\pm\pi$ 足够远，噪声不会把它绕过去，也就是第二道门限之前。本文验证到 $\sigma=0.45$ rad（每频点 3.9 dB）。它是相干 ML 的**模型**，不是界；PM 不适用。

**与 e07 蒙特卡洛对比**〔V2:e07〕〔VER〕（dB；“允许区间”是多谷 − MC 在 95% 带内应落的范围，“—”表示带的一侧超过 20 dB，即 MC 次数太少、未分辨）：

| $\sigma$ (rad) | 0.051 | 0.067 | 0.088 | 0.115 | 0.151 | 0.198 | 0.260 | 0.342 | 0.448 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 几何：MC (dB·m²) | −69.05 | −66.64 | −64.25 | −18.40 | −5.46 | +1.95 | +6.62 | +9.63 | +11.68 |
| 几何：谷混合 − MC | −0.00 | +1.13 | +26.61 | +0.99 | +0.29 | +0.57 | +0.92 | +1.56 | +2.28 |
| 几何：多谷 − MC | −0.00 | +1.13 | +26.59 | +0.96 | +0.15 | +0.14 | +0.07 | +0.03 | +0.13 |
| 几何：允许区间 | −1.34/+1.94 | — | — | — | −0.97/+1.25 | −0.43/+0.48 | −0.25/+0.26 | −0.17/+0.18 | ±0.14 |
| 两簇：MC (dB·m²) | −31.67 | −26.93 | −24.28 | −22.13 | −20.15 | −17.28 | −6.37 | +0.29 | +5.17 |
| 两簇：谷混合 − MC | +0.02 | −0.36 | +0.31 | +1.14 | +2.04 | +4.22 | +3.01 | +3.60 | +3.06 |
| 两簇：多谷 − MC | +0.17 | −0.14 | +0.03 | −0.03 | +0.01 | +1.25 | −0.60 | +0.07 | −0.06 |
| 两簇：允许区间 | −0.41/+0.45 | −0.24/+0.25 | −0.16/+0.17 | ±0.14 | −0.39/+0.43 | −1.61/+2.58 | −0.94/+1.20 | −0.44/+0.49 | −0.26/+0.27 |

多谷模型在全部 18 个点上都落在 MC 的 95% 带内〔VER〕。在带的两侧都不超过 1.5 dB 的点上（几何 5 个、两簇 8 个），多谷的最大误差是几何 0.15 dB、两簇 0.60 dB；谷混合是 2.28 dB 和 3.60 dB，而且在过渡区系统偏悲观。几何设计在 0.067–0.115 rad 上两者都比 MC 高 1–27 dB：那里的滑周来自 6.2 m 远谷的稀有事件，归档的 MC 次数看不到。

**e09 的四个设计**〔V2:e09〕〔VER〕（MC 每点 6000 次）：

| 设计 | MC @0.183 | 谷混合 − MC | 多谷 − MC | MC @0.234 | 谷混合 − MC | 多谷 − MC | $\sigma_{\mathrm{ref}}=0.2$：谷混合 / 多谷 (dB·m²) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 等间距 | −17.11 | +2.80 | −0.02 | −15.04 | +3.12 | −0.04 | −13.44 / −16.38 |
| front_lo | −19.62 | +2.22 | +0.03 | −16.76 | +2.71 | −0.16 | −16.44 / −18.84 |
| front_mid | −5.60 | +1.66 | −1.40† | −0.45 | +2.98 | −0.09 | −1.23 / −4.31 |
| front_hi | +7.30 | +2.44 | +0.04 | +9.20 | +3.08 | −0.00 | +10.70 / +8.05 |

† 重跑 8×6000 次得 −7.10 dB（8 批范围 −8.05 至 −6.24），多谷比重跑值高 0.11 dB〔VER〕。归档的 −5.60 是重尾（远谷）估计的一次涨落。

*结论：* 在 $\sigma\approx0.2$ 上，谷混合对四个设计都悲观 2.2–3.2 dB（front_mid 在 0.183 rad 处用重跑值），在 $\sigma_{\mathrm{ref}}$ 上比多谷高 2.4–3.1 dB，但排序不变（front_lo < 等间距 < front_mid < front_hi）。所以 e08/e09 选出的设计和前沿结论成立，只是归档的 $\mathrm{MSE}_{\mathrm{ref}}$ 绝对值偏高约 2.5–3 dB。**从本文起，谷混合只作为历史搜索代理**：e08/e09 的归档值保持原样，不重写；新预测一律用多谷模型。

---

## 5 分账：ML − PM ≈ 1.2 dB，PM − 最紧已知界 ≈ 0.7 dB（门限上）

### 5.1 门限上的数字（$u=1$，同一均匀先验）

| | 蒙特卡洛 / 数值界〔N〕〔VER〕 | 格点闭式〔VER〕 |
| --- | --- | --- |
| 相干 ML | −23.0859 dB·m²（MMSE +1.13） | MMSE +1.145 |
| 后验均值 PM | −24.3104 | MMSE −0.004 |
| 非相干 | −24.2951 | （＝MMSE，见 5.2） |
| WWB（单测试点，2.5 节） | −24.957 | MMSE −0.743 |
| ZZB | −27.502 | MMSE −3.30 |
| MMSE | −24.2142 | — |
| **ML − PM** | **1.2245 dB** | 1.149 dB |
| **PM − WWB** | **0.647 dB** | 0.738 dB |
| WWB − ZZB | 2.545 dB | 2.556 dB |
| PM − ZZB（旧账） | 3.1916 dB | 3.295 dB |

### 5.2 为什么 PM 就是带宽直线

整数的后验均值（整数上取平坦先验）可以写成

$$
\mathbb E[k\mid z]=z+s^2\partial_z\log\vartheta_s(z),\qquad \vartheta_s(z)=\sum_k e^{-(z-k)^2/2s^2}.
$$

由 Poisson 求和，$\vartheta_s\propto1+2e^{-2\pi^2s^2}\cos2\pi z+\dots$，所以

$$
\mathbb E[k\mid z]=z-4\pi s^2e^{-2\pi^2s^2}\sin2\pi z+\dots
$$

在 $s=\tfrac12$ 时修正项的幅度只有 0.0226〔VER〕。于是 $\mathrm{MSE}_{\mathrm{PM}}\approx\mathrm{HMSE}+\lambda_{\mathrm{eff}}^2s^2=\mathrm{MMSE}$（1.3 节的恒等式）。

*物理意义：* 在 $u\gtrsim1$ 时，最优估计量“不取整”，它就是粗尺。这正是“门限上非相干估计与最优估计一样好”（−24.30 对 −24.31）的数学原因。只有在 $u<1$ 时，PM 才通过整数结构低于带宽直线。

### 5.3 1.2 dB：硬判决的代价

$$
\mathrm{ML}-\mathrm{PM}\approx10\log_{10}\frac{\mathbb E[\mathrm{round}(z)^2]}{\mathbb E[z^2]}\quad(u\gtrsim1).
$$

它就是 3.2 节的过冲。格点模型中，这个差在 $u\approx0.75$ 处达到峰值 1.33 dB〔VER〕；蒙特卡洛的峰值为 1.43 dB，出现在 $u\approx0.80$〔V2:e04〕。

### 5.4 0.7 dB：界的松弛（原 3.2 dB 中 2.5 dB 是 ZZB 自己的）

在这个先验下，PM 是均方误差意义上的最优估计，所以 PM 与任何下界之差全部是**界的松弛**，不是任何估计量能拿回的。旧账用 ZZB，由 2.4 节的闭式，

$$
\frac{\mathrm{ZZB}}{\mathrm{MMSE}}\Big|_{u=1}=\frac1R+\Bigl(1-\frac1R\Bigr)\frac{\sum(2k-1)Q(k/2s)}{2s^2}\Big|_{s=\frac12}\Rightarrow-3.30\ \mathrm{dB}.
$$

单测试点 WWB 在同一点上，由 2.5 节的格点闭式，

$$
\frac{\mathrm{WWB}}{\mathrm{MMSE}}\Big|_{u=1}\approx\Bigl(1-\frac1R\Bigr)\frac{1}{\sinh 1}\Rightarrow-0.70\ \mathrm{dB}\quad(T\to\infty;\ T=14.99\ \mathrm m\ \text{时}-0.743\ \mathrm{dB}).
$$

所以门限上 PM 与最紧已知界只差 0.65 dB（蒙特卡洛 PM）或 0.74 dB（格点 PM）。原来的 3.2 dB 中有 2.5 dB 是 ZZB 自己的松弛：二元检验在 $(0,\lambda]$ 上只积到 $\lambda^2/2$（5.5 节第 2 条），而 WWB 直接以 $h^2=\lambda_{\mathrm{eff}}^2$ 计价。

### 5.5 滑周项的逐段比例

〔VER〕

| $u$ | ML/ZZB | PM/ZZB | ML/PM |
| --- | --- | --- | --- |
| 0.2（脚附近） | 4.00 | 3.12 | 1.28 |
| 0.5 | 3.98 | 3.00 | 1.33 |
| 1.0 | 2.78 | 2.14 | 1.30 |
| 2.0 | 1.60 | 1.48 | 1.08 |

脚附近的因子 4 可以精确拆成两个因子 2：

1. ML 在 $\pm$ 两侧都会滑周（概率 $2Q$），而 ZZB 只做一次二元检验（概率 $Q$）；
2. ML 每次滑周付 $\lambda^2$，而 ZZB 在 $(0,\lambda]$ 上只积到 $\int_0^\lambda h\,dh=\lambda^2/2$。

改用最紧已知界 $\max(\mathrm{ZZB},\mathrm{WWB})$ 时，PM 与它的差（全 MSE，dB；2.5 节表）〔VER〕：

| $u$ | 0.2 | 0.3 | 0.5 | 0.7 | 1.0 | 2.0 |
| --- | --- | --- | --- | --- | --- | --- |
| PM − max(ZZB, WWB) | 1.12 | 4.87 | 4.76 | 2.37 | 0.74 | 0.14 |
| 最紧的界 | ZZB | ZZB | ZZB | WWB | WWB | WWB |

*结论：* ML 与 PM 的差是**估计量**的账，最大约 1.3–1.4 dB。PM 与最紧已知界的差是**界**的账：门限上约 0.7 dB（WWB），$u\ge2$ 时不到 0.3 dB；但在 $u\approx0.3$–0.5 的下坡段仍有约 4.8 dB。那里最紧的仍是 ZZB，滑周项因子约 3。拐点位置不受这两笔账影响，最多差 0.56 dB（3.4 节）。

---

## 6 频点设计的多目标形式化

### 6.0 问题

**决策变量：** $\mathcal F\subset\{2400,2401,\dots,2470\}$ MHz，$|\mathcal F|=N=8$，频点互不相同（整数 MHz，$B\le70$ MHz）。

**目标：**

| 目标 | 量 | 方向 | 由谁决定 |
| --- | --- | --- | --- |
| O1 门限早晚（精度） | $\sigma_{\mathrm{leave}}$ | ↑ | 主要是 $S$；有上界（6.2） |
| O2 门限下平台 | $\mathrm{MSE}(\sigma_{\mathrm{ref}}=0.2)$：归档值用谷混合（搜索代理），预测用多谷模型（4.6） | ↓ | 非格点谷 $\rho_{\mathrm{ng}}$、$D_{\mathrm{ng}}$（6.3） |
| O3 非格点模糊 | $\rho_{\mathrm{ng}}$ ↑；另报近区长度 $L$ ↑ 和 $\Delta_{\mathrm{slip}}$ | ↑ | 间隔结构（6.1、6.3） |
| O4 精确周期（可选） | $P=c/\gcd(\mathcal F)$ | ↑ | 数论（6.5） |

HMSE 不是目标：在固定频段内它几乎不随排布变化（1.2 节）。

**定义。**

- **格点谷与非格点谷。** 因为 $\mathrm{wrap}$ 对每个频点取最短代表，总有 $J\le J_2$。谷 $v$ 称为**格点谷**，若 $J(\delta_v)=J_2(\delta_v)$（实现中 $J_2-J\le10^{-8}$ rad²）：没有任何频点绕离“载波梳子 + 包络抛物线”，于是谷底在 $k\lambda_{\mathrm{eff}}$，深度 $k^2J(\lambda)(1-1/R)$，只由 $(S,\bar f)$ 决定。$J_2-J>0$ 的谷称为**非格点谷**（`design.valley_ratios`）。
- **O3：** $\rho_{\mathrm{ng}}=\min_{v\ \text{非格点},\ \delta_v\in(\lambda/4,T)}J_v/J(\lambda)$，对应位置 $D_{\mathrm{ng}}$。
- **近区长度** $L=c/(2\Delta_{\max})$（6.1 节引理的保证范围）。
- **旧 O3（仅作对比）：** 远区起点 $\delta_B=c/(2B)$，$\rho_{\mathrm{far}}=\min_{v:\,\delta_v>\delta_B}J_v/J(\lambda)$（格点谷与非格点谷都算），位置 $D_{\mathrm{far}}$。
- 第一跳模糊尺度 $\Delta_{\mathrm{slip}}(\sigma)=\arg\max_{v}\,2Q\bigl(\sqrt{J_v}/2\sigma\bigr)\delta_v^2$，即主导滑周预算的那个谷的距离。

**为什么换掉 $\rho_{\mathrm{far}}$。** 格点谷不是排布的性质：它的深度只看 $(S,\bar f)$，它的作用已经计入 O1（脚由 $k=\pm1$ 格点谷决定，6.2 节）和 O2 的格点部分。把它算进 O3，$\rho_{\mathrm{far}}$ 就被 $\delta_B$ 之外的第一个格点齿封顶，约为 $(\delta_B/\lambda)^2$：带宽越窄，$\delta_B$ 越远，$\rho_{\mathrm{far}}$ 越大，与模糊无关。等间距的 $\rho_{\mathrm{far}}=324$ 就是 $k=18$ 的格点齿；它最浅的非格点谷是 $672\,J(\lambda)$（6.1 节）。$\rho_{\mathrm{ng}}$ 量的恰好是格点模型解释不了的那部分模糊。

### 6.1 近区定理：等间距使深谷落在 $k\lambda$

> **引理（近区精确性，任意排布）。** 若 $0\le\delta<c/(2\Delta_{\max})$ 且 $J(\delta)<\pi^2/4$，则 $J(\delta)=J_2(\delta)=N\,\mathrm{wrap}(2\pi\delta/\lambda)^2+(2\pi/c)^2S\delta^2$。
>
> *证明。* 记 $r_i=\mathrm{wrap}(2\pi f_i\delta/c)$。由 $J<\pi^2/4$ 得每个 $|r_i|<\pi/2$，所以 $|r_{i+1}-r_i|<\pi$。另一方面 $r_{i+1}-r_i\equiv2\pi(f_{i+1}-f_i)\delta/c\pmod{2\pi}$，而右边落在 $[0,\pi)$ 内。两个都在 $(-\pi,\pi)$ 内且模 $2\pi$ 同余的数必须相等，所以 $r_i=r_1+2\pi(f_i-f_1)\delta/c$ 对 $f_i$ 是线性的。它的均值 $\bar r\equiv2\pi\bar f\delta/c\pmod{2\pi}$ 且 $|\bar r|<\pi/2$，因此 $\bar r=\mathrm{wrap}(2\pi\delta/\lambda)$。于是 $J=N\bar r^2+\sum_i\bigl(2\pi(f_i-\bar f)\delta/c\bigr)^2=J_2(\delta)$。∎

**推论 1（格点谷）。** 在近区内，每一个 $J<\pi^2/4$ 的谷都位于 $k\lambda_{\mathrm{eff}}$，深度为 $k^2J(\lambda)(1-1/R)$，与排布无关，只看 $(S,\bar f)$。

**推论 2（等间距唯一最小化 $\Delta_{\max}$）。** 因为 $\sum(f_{i+1}-f_i)=B$，有 $\Delta_{\max}\ge B/(N-1)$，等号当且仅当等间距。所以**在给定 $B$ 时**，等间距使“有保证的近区” $L=c(N-1)/(2B)$ 最长。主例中 $L=c/(20\ \mathrm{MHz})=14.99\ \mathrm m=T$：**整个先验窗都是近区**。

**推论 3（近区内非格点谷的深度下限，任意排布）。** 由引理，$\delta<L$ 内的非格点谷都有 $J\ge\pi^2/4$，所以只看近区时 $\rho_{\mathrm{ng}}\ge\pi^2/(4J(\lambda))=\bar f^2/(16S)$。等间距 $L=T$，于是在整个窗上

$$
\rho_{\mathrm{ng}}^{\mathrm{equal}}\ \ge\ \frac{\pi^2/4}{J(\lambda)}=\frac{\bar f^2}{16S}=\frac{R-1}{16N}=88.2 .
$$

数值核验〔VER〕：153 组设计、$1.55\times10^6$ 个点上，$|J-J_2|\le1.8\times10^{-13}$。等间距网格的所有浅谷与 $k\lambda_{\mathrm{eff}}$ 的偏差不超过 0.05 mm，深度相对误差不超过 $4\times10^{-5}$。窗内有 653 个非格点谷，最浅的是 $\rho_{\mathrm{ng}}=672.3$（5.725 m），满足上式。远区（$\delta>2.14$ m）最浅的谷是格点谷 $k=18$，$\rho_{\mathrm{far}}=324=18^2$。

**〔更正〕** 本文初稿说“非格点谷中最浅的是 $100\,J(\lambda)$”，有误：1.231 m 处的 $100\,J(\lambda)$ 是 $k=10$ 的**格点**谷（$J\ge\pi^2/4$，超出引理的条件，但 $J_2-J=7\times10^{-15}$）〔VER〕。

*物理意义：* 等间距时，每一对相邻频点都构成同一把合成尺 $c/\Delta=2T$。没有哪一部分频点能在全体之前“先对齐”。不等间距时，最大间隔 $\Delta_{\max}$ 两侧的频点会以更短的合成波长 $c/\Delta_{\max}$ 各自重新对齐，远谷就从那里长出来。例如纯两端设计 $\Delta_{\max}=64$ MHz，近区只到 2.34 m，远谷出现在 4.494 m ≈ $c/67\ \mathrm{MHz}$（67 MHz 是两簇中心之差）。

*注意：* 这是**保证**，不是必要条件。矩匹配设计 $\Delta_{\max}=28$ MHz，近区只有 5.35 m，但在 15 m 内也没有浅的非格点谷：$\rho_{\mathrm{ng}}=498.9$（11.02 m）。它的 $\rho_{\mathrm{far}}=256$ 是 1.97 m 处 $k=16$ 的格点齿〔VER〕。

### 6.2 精度旋钮：拉开 S

格点区内 $\sigma_{\mathrm{leave}}=\sigma_{\mathrm{jump}}/x(R)=\pi\sqrt S/(\bar f\,x(R))$。它随 $S$ 单调增；$x$ 只对数地依赖 $R$，而 $R$ 随 $S$ 增大而变小，所以还有一点额外收益。

- **上限：** 给定 $B$ 和 $N$（偶数），$S\le NB^2/4$（一半频点放在每一端）。与等间距相比 $S_{\max}/S_{\mathrm{eq}}=3(N-1)/(N+1)=7/3$，即 $\sigma_{\mathrm{jump}}$ 最多 +3.68 dB〔VER〕。
- **整数 MHz 下的两端：** 2400–2403 / 2467–2470 的 $S=8988\ \mathrm{MHz}^2$（2.14 倍），$\sigma_{\mathrm{leave}}$ +3.17 dB；e08 搜到的最优设计 +3.56 dB〔V2:e08〕〔N〕。

> **命题 6.2（O1 上界）。** 令 $\sigma_{\mathrm{cap}}(\mathcal F)$ 为方程 $2Q\bigl(\sqrt{J(\lambda_{\mathrm{eff}})}/2\sigma\bigr)\lambda_{\mathrm{eff}}^2=\mathrm{HMSE}(\sigma)$ 的根，即只由 $k=\pm1$ 格点谷决定的脚（`design.leave_cap`）。则对频段内任意排布，
> $$\sigma_{\mathrm{leave}}(\mathcal F)\ \le\ \sigma_{\mathrm{cap}}(\mathcal F)\ \approx\ \frac{\sigma_{\mathrm{jump}}}{x(R)}\ \le\ \sigma_{\mathrm{cap}}(\text{4+4 两端})=\sigma_{\mathrm{leave}}^{\mathrm{equal}}\cdot10^{3.60/20}.$$
>
> *证明。* (i) $B\le70$ MHz，所以 $L\ge c/(2B)\ge2.14\ \mathrm m\gg\lambda$；又 $J(\lambda_{\mathrm{eff}})\approx4\sigma_{\mathrm{jump}}^2<\pi^2/4$。由 6.1 节引理，每个排布都有位于 $\pm\lambda_{\mathrm{eff}}$、深度 $J(\lambda)(1-1/R)$ 的格点谷。$\sigma_{\mathrm{leave}}$ 按谷和定义（第一个使谷和 $\ge2\,\mathrm{HMSE}$ 的 $\sigma$），其他谷只会让谷和更大，所以在 $\sigma_{\mathrm{cap}}$ 处谷和已经 $\ge2\,\mathrm{HMSE}$，第一次越过只能更早。(ii) $\sigma_{\mathrm{cap}}$ 的方程与 3.3 节相同，$\sigma_{\mathrm{cap}}\approx\pi\sqrt S/(\bar f\,x(R))$ 随 $S$ 增大。8 个互不相同的整数 MHz 落在 2400–2470 时，$S=\sum(f_i-\bar f)^2\le\sum(f_i-2435)^2\le2(35^2+34^2+33^2+32^2)=8988$，等号只在 4+4 两端成立。把 $\bar f$ 的变化也算进去（$S\le8988-8(\bar f-2435)^2$），上确界只比 4+4 两端高 0.0009 dB〔VER〕。∎

数值〔VER〕：

- 4+4 两端的上界是 +3.602 dB（数值求根与 $\sigma_{\mathrm{jump}}/x(R)$ 一致）。两端自己只有 +3.17 dB，因为 4.49 m 的非格点谷（$\rho_{\mathrm{ng}}=1.51$）把脚拉早了。
- e08 最优设计 +3.565 dB，它自己的上界 +3.568 dB，离全局上界只差 **0.037 dB**。
- 表 6.6 的四个前沿点都恰好落在各自的上界上（+2.95、+2.98、+3.33、+3.50 dB）。
- 归档中全部 4527 组设计都满足 $\sigma_{\mathrm{leave}}\le\sigma_{\mathrm{cap}}$（最大比值 $1-6\times10^{-8}$），也都不超过全局上界。

*物理意义：* O1 已经被搜到几乎封顶。在这个频段、$N=8$、整数 MHz 的约束下，任何进一步的搜索最多再推迟门限 0.04 dB；剩下能交换的只有 O2 和 O3。

1500 组随机设计的 $\sigma_{\mathrm{jump}}/\sigma_{\mathrm{leave}}$ 落在 4.80–5.32 之间（中位数 4.98）〔N〕，这是推论 1 的直接后果。

**可检验判据 C1：** 若 $x_{\mathrm{eff}}=\sigma_{\mathrm{jump}}/\sigma_{\mathrm{leave}}\approx x(R)$，说明脚由格点谷决定，O1 是纯精度量；若 $x_{\mathrm{eff}}>x(R)$，说明远谷已经提前把脚拉早。表 6.6 中只有纯两端一组属于后者（4.983 对 4.743）。

### 6.3 模糊旋钮：远谷判据

比较远谷与 $1\lambda$ 谷对滑周预算的贡献。只保留指数项（Chernoff 精度）：

$$
\frac{p_{\mathrm{far}}D^2}{p_\lambda\lambda^2}\approx\Bigl(\frac D\lambda\Bigr)^2\exp\!\left(-\frac{J_{\mathrm{far}}-J(\lambda)}{8\sigma^2}\right).
$$

令它等于 1，并用 $J(\lambda)=4\sigma_{\mathrm{jump}}^2$：

$$
\boxed{\ u_{\mathrm{far}}=\sqrt{\frac{\rho_{\mathrm{far}}-1}{4\ln(D_{\mathrm{far}}/\lambda)}}\ }\qquad\Longleftrightarrow\qquad \rho_{\mathrm{far}}\ \ge\ 1+4u_*^2\ln\frac{D_{\mathrm{far}}}{\lambda}.
$$

*物理意义：* 当 $u<u_{\mathrm{far}}$ 时，第一跳仍是 $1\lambda$（$\Delta_{\mathrm{slip}}=\lambda$）；当 $u>u_{\mathrm{far}}$ 时，远谷接管误差预算。要求“远谷不早于噪声 $u_*$ 接管”，就要求上式右边的 $\rho$。

数值核验：$u_{\mathrm{far}}$ 与逐谷数值交点 $u_{\mathrm{cross}}$ 的对比〔VER〕：

| 设计 | 两端 | e08 最优 | front_hi | front_mid | front_lo | 等间距 |
| --- | --- | --- | --- | --- | --- | --- |
| $u_{\mathrm{far}}$ | 0.19 | 0.27 | 0.36 | 1.02 | 2.05 | 5.29 |
| $u_{\mathrm{cross}}$ | 0.19 | 0.28 | 0.38 | 1.10 | 2.26 | 6.0 |

表中除等间距外，最浅远谷都是非格点谷，即 $\rho_{\mathrm{far}}=\rho_{\mathrm{ng}}$、$D_{\mathrm{far}}=D_{\mathrm{ng}}$〔VER〕，所以本节的 $\rho$ 都可以读作 $\rho_{\mathrm{ng}}$。等间距的 $\rho_{\mathrm{far}}=324$ 是 $k=18$ 的格点齿，格点模型已经计入；它的 $u_{\mathrm{far}}=5.29$ 只是这个齿与 $1\lambda$ 齿的交点，与非格点模糊无关。

**可检验判据 C2：**

| 要求 | 条件（$D/\lambda\approx37$，即 4.5 m） |
| --- | --- |
| 远谷不拉早脚（$u_*=1/x(R)\approx0.21$） | $\rho\gtrsim1.6$ |
| 远谷在 $u_*=0.8$ 之前不接管 | $\rho\gtrsim10$ |
| 远谷在门限中心 $u_*=1$ 之前不接管 | $\rho\gtrsim15$ |
| 在 $\sigma_{\mathrm{ref}}=0.2$ rad（$u_*\approx1.75$）保住平台 | $\rho\gtrsim45$ |

deep 中的工程判据 $J_{\mathrm{far}}/J(\lambda)\gtrsim10$ 对应第二行：远谷在中心之后才接管。第一行解释了 C1：纯两端 $\rho=1.51<1.6$，脚被拉早；e08 最优 $\rho=2.05>1.6$，脚仍由格点决定（$x_{\mathrm{eff}}=4.747\approx x(R)=4.745$）〔VER〕。

### 6.4 锚点：几何条件、混叠与 $K\propto N_{\mathrm{edge}}$

**设定。** 两端各 $n$ 个频点，中心为 $f_{\mathrm{lo}},f_{\mathrm{hi}}$，$\Delta F=f_{\mathrm{hi}}-f_{\mathrm{lo}}$。远谷出现在合成波长的各阶 $\delta_m\approx mc/\Delta F$，$m=1,\dots,M$，$M=\lfloor T\Delta F/c\rfloor$（70 MHz、15 m 时 $M=3$）。它们都是非格点谷，所以本节的 $\rho_{\mathrm{far}}$ 就是 $\rho_{\mathrm{ng}}$。在 $\delta_m$ 处两端同相，稍微挪一下 $\delta$ 就能把公共相位抵消，所以只靠两端时 $\rho_{0,m}$ 很小：$m=1,2,3$ 分别为 0.082、0.184、0.0205〔VER〕。

在相对位置 $\alpha=(f_a-f_{\mathrm{lo}})/\Delta F$ 放一个锚点。两端抵消公共相位之后，锚点在第 $m$ 阶剩下的残差是 $2\pi\lVert m\alpha\rVert$（$\lVert\cdot\rVert$ 表示到最近整数的距离）。若两端足够“硬”，谷底挪不动，则

$$
\boxed{\ \rho_{\mathrm{far}}\approx\min_{1\le m\le M}\Bigl[\rho_{0,m}+\frac{\sum_j(2\pi\lVert m\alpha_j\rVert)^2}{J(\lambda)}\Bigr],\qquad J(\lambda)\approx2\pi^2n\Bigl(\frac{\Delta F}{\bar f}\Bigr)^2 .}
$$

在 $n=64,512$、$\alpha\in\{0.2,0.25,0.5\}$ 上，与数值的差在 2% 以内〔VER〕。由此可以读出三件事。

**(a) 混叠。** 锚点在 $m\alpha\in\mathbb Z$ 的那些阶上是“瞎”的。中点锚点 $\alpha=\tfrac12$ 在 $m=2$ 阶看不见远谷。以 4+4 整数两端加 2435 MHz 为例：$m=1$ 的谷被填到 $\rho=148$（4.487 m），但全窗最浅谷移到 8.988 m，$\rho=6.03$〔VER〕。扫描单个锚点，最好的位置在 $\alpha\approx\tfrac14,\tfrac34$：2419 和 2451 MHz 给 $\rho\approx39$，2418 和 2452 MHz 给 $\approx35$〔VER〕。原因是 $\alpha=\tfrac14$ 时 $\lVert\alpha\rVert,\lVert2\alpha\rVert,\lVert3\alpha\rVert\ge\tfrac14$，三阶都照顾到。

**(b) 均匀插 $K$ 个锚点。** 位置 $\alpha_j=j/(K+1)$ 的锚点在 $m=K+1$ 阶上全体同相，所以必须 $K+1>M$。$M=3$ 时：

| 每端点数 $n$ | $K=1$ | $K=2$ | $K=3$ | $K=4$ | $K=8$ |
| --- | --- | --- | --- | --- | --- |
| 4 | 0.184 | 0.021 | 196 | 210 | 318 |
| 64 | 0.184 | 0.021 | 14.1 | 15.0 | 24.8 |
| 256 | 0.184 | 0.021 | 3.57 | 3.81 | 6.29 |
| 512 | 0.184 | 0.021 | 1.80 | 1.92 | 3.16 |

〔VER〕 $K=1,2$ 两列与 $n$ 无关，就是混叠的直接体现。

**(c) $K\propto N_{\mathrm{edge}}$。** 一个“一般位置”的锚点，平均贡献 $\mathbb E(2\pi\lVert U\rVert)^2=\pi^2/3$（$U$ 均匀分布）。代入上式：

$$
\rho_{\mathrm{far}}\approx\rho_0+\frac{K}{6n}\Bigl(\frac{\bar f}{\Delta F}\Bigr)^2,\qquad
K_{\mathrm{req}}\approx6n(\rho_*-\rho_0)\Bigl(\frac{\Delta F}{\bar f}\Bigr)^2\approx\frac{n}{20}\quad(\rho_*=10,\ \bar f/\Delta F=34.8).
$$

用黄金比例散开的锚点核验〔VER〕：

| $(n,K)$ | (64, 4) | (256, 16) | (512, 16) |
| --- | --- | --- | --- |
| 数值 $\rho$ | 13.75 | 11.93 | 6.01 |
| 公式 | 12.6 | 12.6 | 6.3 |

*物理意义：* 锚点有没有用，是**几何**问题（它在每一阶远谷上的相对相位），与它自己的信噪比无关。两端越厚，$J(\lambda)\propto n$ 越大，锚点的相对作用就越小，所以锚点数要按 $n$ 成比例增加。

**判据 C3：** 对 $m=1,\dots,\lfloor T\Delta F/c\rfloor$，要求 $\sum_j(2\pi\lVert m\alpha_j\rVert)^2\gtrsim\rho_*J(\lambda)$。避免分母 $\le M$ 的有理位置 $\alpha$（例如 $\tfrac12$、$\tfrac13$）。

**〔已修〕deep §12.2 的锚点表**（每端 4/64/256/512，$K=0,1,2,4,8$）在仓库里没有生成脚本。deep 原值与复算对比如下；deep 现已换成复算值，并注明以本文和核验脚本为准：

- $K=8$ 列一致：deep 为 388 / 25 / 6.4 / 3.2，复算为 318 / 24.8 / 6.3 / 3.2；
- $K=4$ 列偏大约 1.5 倍：deep 为 156 / 10 / 2.5 / 1.3，复算为 210 / 15.0 / 3.8 / 1.9；
- $K=1,2$ 两列不一致：deep 中随 $n$ 下降（例如 $K=1$ 为 2.7 / 0.21 / 0.07 / 0.04），复算中因混叠恒为 0.184 和 0.021；
- deep 原来说“只在 2435 MHz 加 1 个锚点，$J$ 升到约 $150J(\lambda)$”，这只对 4.49 m 那个谷成立，全窗最浅谷是 $6.0\,J(\lambda)$（8.99 m）。

定性结论“$K\propto N_{\mathrm{edge}}$，判据看 $J_{\mathrm{far}}/J(\lambda)$”不变。

### 6.5 CRT、GCD 与精确周期

$P=c/\gcd(\mathcal F)$ 只决定 $J$ 在哪里**精确**回到 0。等间距 $P=29.98$ m；表中其他整数设计 $\gcd=1$ MHz，$P=299.8$ m〔VER〕。

- $\sigma_{\mathrm{jump}}$ 和近区（推论 1）只看 $(S,\bar f)$，所以 CRT 推迟不了第一道门限：矩匹配设计与等间距的 $\sigma_{\mathrm{leave}}$ 相同〔V2:e08〕。
- 窗口内的风险来自**近似**重合（丢番图逼近），由 $\rho_{\mathrm{ng}}$ 度量，不由 gcd 度量。例如两个频点 2400 和 2470 MHz：$\gcd=10$ MHz，$P=30$ m，但 12.87 m 处 $\rho_0=0.02$〔VER〕。

**判据 C4：** 只要 $P>2T$，gcd 对 O1–O3 都无影响；CRT 只服务 O4。

### 6.6 帕累托与具体频点组

〔N〕〔V2:e08/e09〕〔VER〕；$T=15$ m，$\sigma_{\mathrm{ref}}=0.2$ rad。“上界”是命题 6.2 的 $\sigma_{\mathrm{cap}}$ 相对等间距 $\sigma_{\mathrm{leave}}$ 的增益。$\mathrm{MSE}_{\mathrm{ref}}$ 列为“归档谷混合 / 多谷模型”，“—”表示归档里没有。

| 设计（MHz） | $S$ (MHz²) | $\sigma_{\mathrm{jump}}$ | $\sigma_{\mathrm{leave}}$（增益） | 上界 | 旧 $\rho_{\mathrm{far}}$ @ $D$ | $\rho_{\mathrm{ng}}$ @ $D_{\mathrm{ng}}$ | $L$ (m) | $u_{\mathrm{far}}$ | $u_{\mathrm{ref}}$ | $\mathrm{MSE}_{\mathrm{ref}}$ (dB·m²) | $\Delta_{\mathrm{slip}}$ @ 脚 / @0.2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 等间距 2400:10:2470 | 4200 | 0.0836 | 0.01703（0） | 0 | 324 @ 2.22 m（格点 $k=18$） | **672 @ 5.73 m** | **14.99** | 5.29 | 2.39 | −13.44 / −16.38 | 0.123 / 0.369 |
| 矩匹配 | 4200 | 0.0836 | 0.01703（0） | 0 | 256 @ 1.97 m（格点 $k=16$） | 499 @ 11.02 m | 5.35 | 4.79 | 2.39 | — / −16.39 | 0.123 / 0.369 |
| front_lo 2400,01,02,**09,50**,68,69,70 | 7826 | 0.1142 | 0.02392（+2.95） | +2.95 | 62.6 @ 4.74 m | 62.6 @ 4.74 m | 3.66 | 2.05 | 1.75 | **−16.44 / −18.84** | 0.123 / 0.246 |
| 2400,01,03,**19,63**,68,69,70 | 7894 | 0.1146 | 0.02400（+2.98） | +2.98 | 51.2 @ 4.74 m | 51.2 @ 4.74 m | 3.41 | 1.85 | 1.75 | −16.11 / −18.58 | 0.123 / 0.246 |
| front_mid 2400–03,**58**,68–70 | 8483 | 0.1189 | 0.02500（+3.33） | +3.33 | 16.2 @ 4.62 m | 16.2 @ 4.62 m | 2.73 | 1.02 | 1.68 | −1.23 / −4.31 | 0.123 / 4.62 |
| front_hi 2400,01,03,04,**66**,68–70 | 8797 | 0.1210 | 0.02548（+3.50） | +3.50 | 2.91 @ 4.49 m | 2.91 @ 4.49 m | 2.42 | 0.36 | 1.65 | +10.70 / +8.05 | 0.123 / 4.49 |
| e08 最优 2400–02,04,67–70 | 8925 | 0.1219 | 0.02568（+3.56） | +3.57 | 2.05 @ 4.49 m | 2.05 @ 4.49 m | 2.38 | 0.27 | 1.64 | — / +9.17 | 0.123 / 4.49 |
| 纯两端 2400–03,67–70 | 8988 | 0.1223 | 0.02455（+3.17） | **+3.60** | 1.51 @ 4.49 m | 1.51 @ 4.49 m | 2.34 | 0.19 | 1.64 | — / +10.33 | **4.49** / 8.99 |

**读表：**

1. **(O1, O2) 平面。** 等间距被 front_lo 和 2419/2463 组严格支配；用多谷模型的 $\mathrm{MSE}_{\mathrm{ref}}$ 也一样（−18.84、−18.58 对 −16.38）。按归档的谷混合 O2，3000 组随机设计中有 12.2% 也支配等间距〔VER〕。随机设计 $\sigma_{\mathrm{leave}}$ 的中位数为 0.0141，$\mathrm{MSE}_{\mathrm{ref}}$ 的中位数为 −11.7 dB·m²〔VER〕。**〔已修〕** deep §12.2 原写 −11.6；本文按 e09 云图取中位数再换算成 dB。
2. **前沿就是 O1 与 $\rho_{\mathrm{ng}}$ 的交换。** 前沿设计的最浅远谷都是非格点谷（$\rho_{\mathrm{far}}=\rho_{\mathrm{ng}}$）。沿前沿，$\sigma_{\mathrm{leave}}$ 从 +2.95 dB 走到 +3.50 dB，每一点都贴着自己的上界；$\rho_{\mathrm{ng}}$ 则从 62.6 掉到 2.91。最后约 0.5 dB 的门限推迟，是把锚点挪到两端换来的，代价是非格点远谷。O2 的好坏完全由 $u_{\mathrm{far}}$ 与 $u_{\mathrm{ref}}$ 的大小关系预测：$u_{\mathrm{far}}>u_{\mathrm{ref}}$ 的组（等间距、front_lo、2419 组）平台好，其余平台差（6.3 节）。
3. **(O1, O2, O3) 空间：重判。** 在 e09 归档的 3017 组互不相同的设计上〔VER〕：
   - **旧 O3（$\rho_{\mathrm{far}}$）下，等间距被 61 组支配。** 2317 组的 $\rho_{\mathrm{far}}\ge324$，因为等间距的 324 只是 $k=18$ 格点齿（6.0 节）。例如 2407,14,15,26,38,52,69,70：$\sigma_{\mathrm{leave}}=0.01739$，归档 $\mathrm{MSE}_{\mathrm{ref}}=-13.63$，$\rho_{\mathrm{far}}=400$ 是 2.46 m 处 $k=20$ 的格点齿，而它的 $\rho_{\mathrm{ng}}=409<672$。**本文初稿“等间距的 $\rho_{\mathrm{far}}$ 最大、不被支配”的说法不成立。**
   - **新 O3（$\rho_{\mathrm{ng}}$）下，0 组支配等间距。** 更强的是：没有任何设计在 $(\sigma_{\mathrm{leave}},\rho_{\mathrm{ng}})$ 两项上同时不差于等间距。所以这个结论与 O2 用哪个代理（归档谷混合或多谷模型）无关。
   - **等间距不是 $\rho_{\mathrm{ng}}$ 的全局最大者。** 1005 组（33%）的 $\rho_{\mathrm{ng}}\ge672$，但它们都是窄带设计（$B$ 中位数 50 MHz），$\sigma_{\mathrm{leave}}$ 都比等间距低至少 0.6 dB。在 (O1, O2) 上胜过等间距的 380 组中，$\rho_{\mathrm{ng}}$ 最大只有 518（10.59 m）。
   - **近区长度 $L$。** 推论 2 是固定 $B$ 的结论：所有非等间距设计都有 $L<c(N-1)/(2B)$。不固定 $B$ 时，紧凑设计的 $L$ 可达 29.98 m；但在 $\sigma_{\mathrm{leave}}\ge$ 等间距的设计中，$L$ 最大只有 11.53 m，不到等间距的 14.99 m。

   所以 deep 说的“等间距不是帕累托最优”只在二维目标 (O1, O2) 上成立。加上 O3＝$\rho_{\mathrm{ng}}$，等间距在归档上不被支配，是前沿的一个端点：在 $\sigma_{\mathrm{leave}}$ 不低于它的设计里，它的 $\rho_{\mathrm{ng}}$ 和 $L$ 都最大。这是归档上的结论，不是对全部 $\binom{71}{8}$ 种组合的证明（第 7 节）。
4. **蒙特卡洛核验**〔V2:e09〕：

| 每频点 SNR | 等间距 | front_lo |
| --- | --- | --- |
| 29.0 dB | −52.9 dB·m² | −75.3 dB·m² |
| 18.2 dB | −22.8 dB·m² | −25.3 dB·m² |
| 7.4 dB | −13.0 dB·m² | −8.4 dB·m² |

在 7.4 dB 时 front_lo 输给等间距，这是它的两簇结构触发了第二层门限（4.4 节）。

**各类排布分别优化什么：**

| 排布 | 优化的目标 | 代价 / 风险 | 判据 |
| --- | --- | --- | --- |
| 等间距 | O3：$\rho_{\mathrm{ng}}$ 与 $L$ 有保证；O2 不错 | O1 中等（离上界 3.6 dB） | 推论 2、3，$\rho_{\mathrm{ng}}\ge(R-1)/16N$ |
| 两端 + 锚点 | O1（+3 dB，贴住各自上界）且保住 O2 | 锚点混叠；锚点数要 $\propto n$ | C2、C3 |
| 纯两端 | $S$ 最大（上界 +3.60 dB） | O3 崩溃（$\rho_{\mathrm{ng}}=1.51$，$\Delta_{\mathrm{slip}}=4.49$ m），脚被拉早到 +3.17 dB | $\rho<1.6$ |
| CRT / 互素 | O4 | 对 O1–O3 无保证 | C4 |
| 均匀随机 | 无 | 中位数在 O1、O2 上都不如等间距 | — |

---

## 7 仍开放（研究缺口）

1. **下坡段的界缺口。** 门限上已由单测试点 WWB 收到 0.7 dB（2.5 节），$u\ge1.5$ 时不到 0.3 dB。没收回的是 $u\approx0.3$–0.5 处的约 4.8 dB：那里 ZZB 最紧，单测试点 WWB 因指数翻倍而无用。候选是多测试点 WWB（以 $\{\pm k\lambda_{\mathrm{eff}}\}$ 为联合测试点的矩阵形式），以及 Bell 扩展 ZZB 或多元假设版本的 ZZB。5.5 节把因子 4 拆成“双侧”和“$\int h\,dh$”两个因子 2：后一个正是 WWB 在门限上修掉的，前一个需要双侧或多元检验。
2. **非格点设计的逐谷 ZZB 闭式。** 推广形式 $\mathrm{ZZB}\approx\mathrm{HMSE}+\sum_v\tfrac12(\delta_v^2-\delta_{v-1}^2)\,Q(\sqrt{J_v}/2\sigma)$（$v$ 取记录谷）尚未在几何间隔、两簇设计上验证。
3. **多谷模型的解析误差界与快速近似。** 多谷线性化模型（4.6 节）已取代谷混合，在 e07/e09 上都落在 MC 带内；但它的误差只有经验核验，没有解析界（线性化在残差接近 $\pm\pi$ 时和第二道门限附近失效）。它目前靠重要性抽样，每个点约 $10^5$ 个样本，当搜索目标太慢，需要确定性近似（例如成对相关修正的 Bonferroni 型界或 Genz 积分）。用它当 O2 重跑 e08/e09 前沿也还没做：因为 $\sigma_{\mathrm{ref}}$ 上的排序不变，预期前沿成员基本不变，只是 O2 整体下移约 2.5–3 dB。
4. **(O1, O3) 前沿的证明。** 在 $\rho_{\mathrm{ng}}$ 下，等间距在 3017 组归档设计上不被支配（6.6 节），但这只是归档上的结论。“$\sigma_{\mathrm{leave}}\ge\sigma_{\mathrm{leave}}^{\mathrm{equal}}\Rightarrow\rho_{\mathrm{ng}}<\rho_{\mathrm{ng}}^{\mathrm{equal}}$”对全部组合是否成立，没有证明。在给定 $(N,B,T)$ 下 $\max_{\mathcal F}\rho_{\mathrm{ng}}$ 的结构也未知；窄带设计说明它不由等间距达到。
5. **精确组合最优。** $\binom{71}{8}\approx1.06\times10^{10}$ 种组合。近区引理和 O1 上界（命题 6.2：$\sigma_{\mathrm{leave}}\le\sigma_{\mathrm{cap}}$，只依赖 $(S,\bar f)$）都可以用来剪枝。还要加入蓝牙信道约束、其他 $N$ 与 $B$，以及三目标或四目标的完整前沿。
6. **相关相位噪声与不等 SNR。** 共用本振时 $\Sigma$ 非对角，衰落时权重不等。`lattice.params` 已支持三明治形式，但没有蒙特卡洛验证；加权情形下近区引理的条件 $J_W<\pi^2/4\cdot\min w_i$ 也需要重新推导。
7. **多条多径。** 单反射下 $\beta(\Delta)$ 已有闭式；多条反射、时变多径下 $P(|\beta|>\tfrac12)$ 的统计，以及以此为目标的排频，都还没做。
8. **第二道门限的统一模型。** 用复观测把 Rife–Boorstyn 门限与格点模型接成一条曲线（目前相位模型在每频点低于约 6 dB 时失效，deep §7）。
9. **先验依赖。** 本文全部贝叶斯数字都用 $T=P/2$。当 $T>P/2$ 时会出现周期混叠；先验非均匀时 PM ≈ MMSE 的结论需要重新检查。
