# 多频相位测距的门限：数学专文

*本页只讲数学：符号、界、拐点、格点模型、分账、频点设计的多目标形式。叙事版见 [README.md](README.md)，推导与实验全文见 [deep.zh.md](deep.zh.md)。*

**数字来源标记。**
〔N〕＝[`numbers.json`](numbers.json)（由 `computation/article_numbers.py` 从归档算出）；
〔V2:eXX〕＝`computation/data/v2/eXX_*.json`；
〔VER〕＝本文新增的核验脚本 [`computation/verify_math_idents.py`](../computation/verify_math_idents.py) 的输出：它只读归档与库函数，不写任何文件；运行 `python3 computation/verify_math_idents.py`，全部 PASS 约 10 s。
与已有文字不一致之处在文中用 **〔不一致〕** 标出。

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

把真值 $d$ 与 $d+h$ 当作两个等概假设，在高斯相位向量下，最小错误概率为

$$
P_e(h)=Q\!\left(\frac{\sqrt{J(h)}}{2\sigma}\right).
$$

谷填充的 Ziv–Zakai 下界为

$$
\mathrm{ZZB}=\frac1T\int_0^T h\;\mathcal V\!\left[(T-h)P_e(h)\right]dh,\qquad \mathcal V g(h)=\max_{v\ge h}g(v).
$$

实现见 `mfpr/bounds.py`。计算时每个谷底都单独精细化；否则在低噪声端会错 30 dB 以上（deep 附录 C）。

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

**〔不一致〕** deep §4.2 和 §11 的“安全规则：总 SNR $\gtrsim R_{\mathrm{dB}}+1$ dB”只在 $R\approx40$–45 dB 时成立。一般情形应当用上表，或直接用上式 $R_{\mathrm{dB}}-12.95+20\log_{10}x(R)$。

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

条件是：所有起作用的谷都在 $k\lambda_{\mathrm{eff}}$ 上。第 6.1 节的定理给出一个**充分**条件：$\delta<c/(2\Delta_{\max})$。等间距时这个条件覆盖整个先验窗。

---

## 5 分账：ML − PM ≈ 1.2 dB，PM − ZZB ≈ 3.2 dB（门限上）

### 5.1 门限上的数字（$u=1$，同一均匀先验）

| | 蒙特卡洛〔N〕 | 格点闭式〔VER〕 |
| --- | --- | --- |
| 相干 ML | −23.0859 dB·m²（MMSE +1.13） | MMSE +1.145 |
| 后验均值 PM | −24.3104 | MMSE −0.004 |
| 非相干 | −24.2951 | （＝MMSE，见 5.2） |
| ZZB | −27.502 | MMSE −3.30 |
| MMSE | −24.2142 | — |
| **ML − PM** | **1.2245 dB** | 1.149 dB |
| **PM − ZZB** | **3.1916 dB** | 3.295 dB |

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

### 5.4 3.2 dB：界的松弛

在这个先验下，PM 是均方误差意义上的最优估计，所以 PM − ZZB 全部是**界的松弛**，不是任何估计量能拿回的。由 2.4 节的闭式，

$$
\frac{\mathrm{ZZB}}{\mathrm{MMSE}}\Big|_{u=1}=\frac1R+\Bigl(1-\frac1R\Bigr)\frac{\sum(2k-1)Q(k/2s)}{2s^2}\Big|_{s=\frac12}\Rightarrow-3.30\ \mathrm{dB}.
$$

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

*结论：* ML 与 PM 的差是**估计量**的账，最大约 1.3–1.4 dB；PM 与 ZZB 的差是**界**的账，门限上约 3.2–3.3 dB，脚附近约 5 dB（因子 3.1）。拐点位置不受这两笔账影响，最多差 0.56 dB（3.4 节）。

---

## 6 频点设计的多目标形式化

### 6.0 问题

**决策变量：** $\mathcal F\subset\{2400,2401,\dots,2470\}$ MHz，$|\mathcal F|=N=8$，频点互不相同（整数 MHz，$B\le70$ MHz）。

**目标：**

| 目标 | 量 | 方向 | 由谁决定 |
| --- | --- | --- | --- |
| O1 门限早晚（精度） | $\sigma_{\mathrm{leave}}$ | ↑ | 主要是 $S$（6.2） |
| O2 门限下平台 | 谷混合 $\mathrm{MSE}(\sigma_{\mathrm{ref}}=0.2)$ | ↓ | 远谷 $\rho_{\mathrm{far}}$、$D_{\mathrm{far}}$（6.3） |
| O3 第一跳模糊尺度 | $\Delta_{\mathrm{slip}}$ ↓，$\rho_{\mathrm{far}}$ ↑ | | 间隔结构（6.1、6.3） |
| O4 精确周期（可选） | $P=c/\gcd(\mathcal F)$ | ↑ | 数论（6.5） |

HMSE 不是目标：在固定频段内它几乎不随排布变化（1.2 节）。

**定义。**

- 远区起点 $\delta_B=c/(2B)$。
- 远谷比 $\rho_{\mathrm{far}}=\min_{v:\,\delta_v>\delta_B}J_v/J(\lambda)$，对应位置记为 $D_{\mathrm{far}}$。
- 第一跳模糊尺度 $\Delta_{\mathrm{slip}}(\sigma)=\arg\max_{v}\,2Q\bigl(\sqrt{J_v}/2\sigma\bigr)\delta_v^2$，即主导滑周预算的那个谷的距离。

### 6.1 近区定理：等间距使深谷落在 $k\lambda$

> **引理（近区精确性，任意排布）。** 若 $0\le\delta<c/(2\Delta_{\max})$ 且 $J(\delta)<\pi^2/4$，则 $J(\delta)=J_2(\delta)=N\,\mathrm{wrap}(2\pi\delta/\lambda)^2+(2\pi/c)^2S\delta^2$。
>
> *证明。* 记 $r_i=\mathrm{wrap}(2\pi f_i\delta/c)$。由 $J<\pi^2/4$ 得每个 $|r_i|<\pi/2$，所以 $|r_{i+1}-r_i|<\pi$。另一方面 $r_{i+1}-r_i\equiv2\pi(f_{i+1}-f_i)\delta/c\pmod{2\pi}$，而右边落在 $[0,\pi)$ 内。两个都在 $(-\pi,\pi)$ 内且模 $2\pi$ 同余的数必须相等，所以 $r_i=r_1+2\pi(f_i-f_1)\delta/c$ 对 $f_i$ 是线性的。它的均值 $\bar r\equiv2\pi\bar f\delta/c\pmod{2\pi}$ 且 $|\bar r|<\pi/2$，因此 $\bar r=\mathrm{wrap}(2\pi\delta/\lambda)$。于是 $J=N\bar r^2+\sum_i\bigl(2\pi(f_i-\bar f)\delta/c\bigr)^2=J_2(\delta)$。∎

**推论 1（格点谷）。** 在近区内，每一个 $J<\pi^2/4$ 的谷都位于 $k\lambda_{\mathrm{eff}}$，深度为 $k^2J(\lambda)(1-1/R)$，与排布无关，只看 $(S,\bar f)$。

**推论 2（等间距唯一最小化 $\Delta_{\max}$）。** 因为 $\sum(f_{i+1}-f_i)=B$，有 $\Delta_{\max}\ge B/(N-1)$，等号当且仅当等间距。所以等间距使“有保证的近区” $L=c(N-1)/(2B)$ 最长。主例中 $L=c/(20\ \mathrm{MHz})=14.99\ \mathrm m=T$：**整个先验窗都是近区**。在这个窗里，任何 $J<\pi^2/4$ 的点都在某个 $k\lambda$ 齿内，并且

$$
\rho_{\mathrm{far}}^{\mathrm{equal}}\ \ge\ \frac{\pi^2/4}{J(\lambda)}=\frac{\bar f^2}{16S}=\frac{R-1}{16N}=88.2 .
$$

数值核验〔VER〕：153 组设计、$1.55\times10^6$ 个点上，$|J-J_2|\le1.8\times10^{-13}$。等间距网格的所有浅谷与 $k\lambda_{\mathrm{eff}}$ 的偏差不超过 0.05 mm，深度相对误差不超过 $4\times10^{-5}$；非格点谷中最浅的是 $100\,J(\lambda)$；远区（$\delta>2.14$ m）最浅的仍是格点谷 $k=18$，$\rho=324=18^2$。

*物理意义：* 等间距时，每一对相邻频点都构成同一把合成尺 $c/\Delta=2T$。没有哪一部分频点能在全体之前“先对齐”。不等间距时，最大间隔 $\Delta_{\max}$ 两侧的频点会以更短的合成波长 $c/\Delta_{\max}$ 各自重新对齐，远谷就从那里长出来。例如纯两端设计 $\Delta_{\max}=64$ MHz，近区只到 2.34 m，远谷出现在 4.494 m ≈ $c/67\ \mathrm{MHz}$（67 MHz 是两簇中心之差）。

*注意：* 这是**保证**，不是必要条件。矩匹配设计 $\Delta_{\max}=28$ MHz，近区只有 5.35 m，但在 15 m 内也没有浅远谷（$\rho_{\mathrm{far}}=256$）〔VER〕。

### 6.2 精度旋钮：拉开 S

格点区内 $\sigma_{\mathrm{leave}}=\sigma_{\mathrm{jump}}/x(R)=\pi\sqrt S/(\bar f\,x(R))$。它随 $S$ 单调增；$x$ 只对数地依赖 $R$，而 $R$ 随 $S$ 增大而变小，所以还有一点额外收益。

- **上限：** 给定 $B$ 和 $N$（偶数），$S\le NB^2/4$（一半频点放在每一端）。与等间距相比 $S_{\max}/S_{\mathrm{eq}}=3(N-1)/(N+1)=7/3$，即 $\sigma_{\mathrm{jump}}$ 最多 +3.68 dB〔VER〕。
- **整数 MHz 下的两端：** 2400–2403 / 2467–2470 的 $S=8988\ \mathrm{MHz}^2$（2.14 倍），$\sigma_{\mathrm{leave}}$ +3.17 dB；e08 搜到的最优设计 +3.56 dB〔V2:e08〕〔N〕。

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

**可检验判据 C2：**

| 要求 | 条件（$D/\lambda\approx37$，即 4.5 m） |
| --- | --- |
| 远谷不拉早脚（$u_*=1/x(R)\approx0.21$） | $\rho\gtrsim1.6$ |
| 远谷在 $u_*=0.8$ 之前不接管 | $\rho\gtrsim10$ |
| 远谷在门限中心 $u_*=1$ 之前不接管 | $\rho\gtrsim15$ |
| 在 $\sigma_{\mathrm{ref}}=0.2$ rad（$u_*\approx1.75$）保住平台 | $\rho\gtrsim45$ |

deep 中的工程判据 $J_{\mathrm{far}}/J(\lambda)\gtrsim10$ 对应第二行：远谷在中心之后才接管。第一行解释了 C1：纯两端 $\rho=1.51<1.6$，脚被拉早；e08 最优 $\rho=2.05>1.6$，脚仍由格点决定（$x_{\mathrm{eff}}=4.747\approx x(R)=4.745$）〔VER〕。

### 6.4 锚点：几何条件、混叠与 $K\propto N_{\mathrm{edge}}$

**设定。** 两端各 $n$ 个频点，中心为 $f_{\mathrm{lo}},f_{\mathrm{hi}}$，$\Delta F=f_{\mathrm{hi}}-f_{\mathrm{lo}}$。远谷出现在合成波长的各阶 $\delta_m\approx mc/\Delta F$，$m=1,\dots,M$，$M=\lfloor T\Delta F/c\rfloor$（70 MHz、15 m 时 $M=3$）。在 $\delta_m$ 处两端同相，稍微挪一下 $\delta$ 就能把公共相位抵消，所以只靠两端时 $\rho_{0,m}$ 很小：$m=1,2,3$ 分别为 0.082、0.184、0.0205〔VER〕。

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

**〔不一致〕deep §12.2 的锚点表**（每端 4/64/256/512，$K=0,1,2,4,8$）在仓库里没有生成脚本。复算结果：

- $K=8$ 列一致：deep 为 388 / 25 / 6.4 / 3.2，复算为 318 / 24.8 / 6.3 / 3.2；
- $K=4$ 列偏大约 1.5 倍：deep 为 156 / 10 / 2.5 / 1.3，复算为 210 / 15.0 / 3.8 / 1.9；
- $K=1,2$ 两列不一致：deep 中随 $n$ 下降（例如 $K=1$ 为 2.7 / 0.21 / 0.07 / 0.04），复算中因混叠恒为 0.184 和 0.021；
- deep 所说“只在 2435 MHz 加 1 个锚点，$J$ 升到约 $150J(\lambda)$”只对 4.49 m 那个谷成立，全窗最浅谷是 $6.0\,J(\lambda)$（8.99 m）。

定性结论“$K\propto N_{\mathrm{edge}}$，判据看 $J_{\mathrm{far}}/J(\lambda)$”不变；具体数字以本文复算为准。

### 6.5 CRT、GCD 与精确周期

$P=c/\gcd(\mathcal F)$ 只决定 $J$ 在哪里**精确**回到 0。等间距 $P=29.98$ m；表中其他整数设计 $\gcd=1$ MHz，$P=299.8$ m〔VER〕。

- $\sigma_{\mathrm{jump}}$ 和近区（推论 1）只看 $(S,\bar f)$，所以 CRT 推迟不了第一道门限：矩匹配设计与等间距的 $\sigma_{\mathrm{leave}}$ 相同〔V2:e08〕。
- 窗口内的风险来自**近似**重合（丢番图逼近），由 $\rho_{\mathrm{far}}$ 度量，不由 gcd 度量。例如两个频点 2400 和 2470 MHz：$\gcd=10$ MHz，$P=30$ m，但 12.87 m 处 $\rho_0=0.02$〔VER〕。

**判据 C4：** 只要 $P>2T$，gcd 对 O1–O3 都无影响；CRT 只服务 O4。

### 6.6 帕累托与具体频点组

〔N〕〔V2:e08/e09〕〔VER〕；$T=15$ m，$\sigma_{\mathrm{ref}}=0.2$ rad。

| 设计（MHz） | $S$ (MHz²) | $\sigma_{\mathrm{jump}}$ | $\sigma_{\mathrm{leave}}$（增益） | $\rho_{\mathrm{far}}$ @ $D$ | $u_{\mathrm{far}}$ | $u_{\mathrm{ref}}$ | $\mathrm{MSE}_{\mathrm{ref}}$ (dB·m²) | $\Delta_{\mathrm{slip}}$ @ 脚 / @0.2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 等间距 2400:10:2470 | 4200 | 0.0836 | 0.01703（0） | 324 @ 2.22 m（格点） | 5.29 | 2.39 | −13.44 | 0.123 / 0.369 |
| 矩匹配 | 4200 | 0.0836 | 0.01703（0） | 256 @ 1.97 m | 4.79 | 2.39 | — | 0.123 / 0.369 |
| front_lo 2400,01,02,**09,50**,68,69,70 | 7826 | 0.1142 | 0.02392（+2.95） | 62.6 @ 4.74 m | 2.05 | 1.75 | **−16.44** | 0.123 / 0.246 |
| 2400,01,03,**19,63**,68,69,70 | 7894 | 0.1146 | 0.02400（+2.98） | 51.2 @ 4.74 m | 1.85 | 1.75 | −16.11 | 0.123 / 0.246 |
| front_mid 2400–03,**58**,68–70 | 8483 | 0.1189 | 0.02500（+3.33） | 16.2 @ 4.62 m | 1.02 | 1.68 | −1.23 | 0.123 / 4.62 |
| front_hi 2400,01,03,04,**66**,68–70 | 8797 | 0.1210 | 0.02548（+3.50） | 2.91 @ 4.49 m | 0.36 | 1.65 | +10.70 | 0.123 / 4.49 |
| e08 最优 2400–02,04,67–70 | 8925 | 0.1219 | 0.02568（+3.56） | 2.05 @ 4.49 m | 0.27 | 1.64 | — | 0.123 / 4.49 |
| 纯两端 2400–03,67–70 | 8988 | 0.1223 | 0.02455（+3.17） | 1.51 @ 4.49 m | 0.19 | 1.64 | — | **4.49** / 8.99 |

**读表：**

1. **(O1, O2) 平面。** 等间距被 front_lo 和 2419/2463 组严格支配。3000 组随机设计中有 12.2% 也支配等间距〔VER〕。随机设计 $\sigma_{\mathrm{leave}}$ 的中位数为 0.0141，$\mathrm{MSE}_{\mathrm{ref}}$ 的中位数为 −11.7 dB·m²〔VER〕。**〔不一致，轻微〕** deep §12.2 写的是 −11.6；本文按 e09 云图取中位数再换算成 dB。
2. **前沿就是 O1 与 $\rho_{\mathrm{far}}$ 的交换。** 沿前沿，$\sigma_{\mathrm{leave}}$ 从 +2.95 dB 走到 +3.50 dB，$\rho_{\mathrm{far}}$ 则从 62.6 掉到 2.91。最后约 0.5 dB 的门限推迟，是把锚点挪到两端换来的，代价是远谷。O2 的好坏完全由 $u_{\mathrm{far}}$ 与 $u_{\mathrm{ref}}$ 的大小关系预测：$u_{\mathrm{far}}>u_{\mathrm{ref}}$ 的组（等间距、front_lo、2419 组）平台好，其余平台差（6.3 节）。
3. **(O1, O2, O3) 空间。** 等间距**不被支配**：它的 $\rho_{\mathrm{far}}$ 最大，而且 $\ge88$ 有定理保证；它在全窗上 $\Delta_{\mathrm{slip}}\in\{k\lambda\}$。所以 deep 说的“等间距不是帕累托最优”只在二维目标 (O1, O2) 上成立；加上 O3，等间距是前沿的一个端点。
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
| 等间距 | O3（有保证）；O2 不错 | O1 中等 | 推论 2，$\rho\ge(R-1)/16N$ |
| 两端 + 锚点 | O1（+3 dB）且保住 O2 | 锚点混叠；锚点数要 $\propto n$ | C2、C3 |
| 纯两端 | O1 最大 | O3 崩溃（$\Delta_{\mathrm{slip}}=4.49$ m），脚也被拉早 | $\rho<1.6$ |
| CRT / 互素 | O4 | 对 O1–O3 无保证 | C4 |
| 均匀随机 | 无 | 中位数在 O1、O2 上都不如等间距 | — |

---

## 7 仍开放（研究缺口）

1. **更紧的全局界。** Weiss–Weinstein 界、Bell 扩展 ZZB，以及在 $h=k\lambda$ 处取测试点的格点特化版本，能否收回门限上的 3.3 dB 和脚附近的因子 3.1？5.5 节把因子 4 拆成了“双侧”和“$\int h\,dh$”两部分，给出了改进的方向。
2. **非格点设计的逐谷 ZZB 闭式。** 推广形式 $\mathrm{ZZB}\approx\mathrm{HMSE}+\sum_v\tfrac12(\delta_v^2-\delta_{v-1}^2)\,Q(\sqrt{J_v}/2\sigma)$（$v$ 取记录谷）尚未在几何间隔、两簇设计上验证。
3. **谷混合近似的误差界。** 在深过渡区它偏保守 1–4.5 dB（deep §9.1），而 e08/e09 的搜索目标正是它。需要一个带误差界的替代，或者直接用格点 ZZB 当目标。
4. **等间距是否唯一最大化 $\rho_{\mathrm{far}}$。** 推论 2 只证明等间距唯一最大化**有保证的**近区；矩匹配的例子说明其他设计也可能有大的 $\rho_{\mathrm{far}}$。在给定 $(N,B,T)$ 下，$\max_{\mathcal F}\rho_{\mathrm{far}}$ 的结构未知。
5. **精确组合最优。** $\binom{71}{8}\approx1.06\times10^{10}$ 种组合；近区引理可以用来剪枝。还要加入蓝牙信道约束、其他 $N$ 与 $B$，以及三目标或四目标的完整前沿。
6. **相关相位噪声与不等 SNR。** 共用本振时 $\Sigma$ 非对角，衰落时权重不等。`lattice.params` 已支持三明治形式，但没有蒙特卡洛验证；加权情形下近区引理的条件 $J_W<\pi^2/4\cdot\min w_i$ 也需要重新推导。
7. **多条多径。** 单反射下 $\beta(\Delta)$ 已有闭式；多条反射、时变多径下 $P(|\beta|>\tfrac12)$ 的统计，以及以此为目标的排频，都还没做。
8. **第二道门限的统一模型。** 用复观测把 Rife–Boorstyn 门限与格点模型接成一条曲线（目前相位模型在每频点低于约 6 dB 时失效，deep §7）。
9. **先验依赖。** 本文全部贝叶斯数字都用 $T=P/2$。当 $T>P/2$ 时会出现周期混叠；先验非均匀时 PM ≈ MMSE 的结论需要重新检查。
