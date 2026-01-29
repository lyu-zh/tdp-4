# SOCP近似约束分析：r=4, p=3 的情况

## 问题重述

代码中验证的是SOCP近似约束：
\begin{align}
\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} + \sqrt{\frac{1-\gamma_a}{\gamma_a}} \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}} &\leqslant 0, \label{eq:assign_dep_rel_balance_d1_adm_lower} \\
\boldsymbol{\mu}_w^\top \mathbf{v}_{j,U} + \sqrt{\frac{1-\gamma_b}{\gamma_b}} \sqrt{\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U}} &\leqslant 0. \label{eq:assign_dep_rel_balance_d1_adm_upper}
\end{align}

当 r=4, p=3 时，确定性约束总是满足。那么SOCP近似约束是否也总是满足？

## 关键理解：SOCP近似约束是保守的

### 1. SOCP近似约束的来源

SOCP近似约束是通过**分解方法（ADM）**得到的：

原始精确约束：
$$\frac{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L} + (\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L})^2} + \frac{\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U}}{\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U} + (\boldsymbol{\mu}_w^\top \mathbf{v}_{j,U})^2} \leqslant \gamma$$

被分解为两个独立的约束：
- $\frac{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L} + (\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L})^2} \leqslant \gamma_a$
- $\frac{\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U}}{\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U} + (\boldsymbol{\mu}_w^\top \mathbf{v}_{j,U})^2} \leqslant \gamma_b$

其中 $\gamma_a + \gamma_b = \gamma$。

然后每个约束被转换为SOCP形式（使用Cantelli不等式）。

### 2. SOCP近似约束的数学含义

SOCP约束：
$$\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} + \sqrt{\frac{1-\gamma_a}{\gamma_a}} \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}} \leqslant 0$$

等价于：
$$\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} \leqslant -\sqrt{\frac{1-\gamma_a}{\gamma_a}} \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}$$

这意味着：**均值项必须足够负，以抵消不确定性项的影响**。

### 3. 为什么即使确定性约束满足，SOCP近似约束仍可能不满足？

#### 关键点：SOCP近似约束是保守的

即使确定性约束总是满足：
- 对于任意分布，$\widetilde{\mathbf{d}}_w^\top \mathbf{v}_{j,L} \leqslant 0$ 总是成立
- 所以 $\mathbb{E}[\widetilde{\mathbf{d}}_w^\top \mathbf{v}_{j,L}] = \boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} \leqslant 0$

**但是**，SOCP近似约束要求：
$$\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} + \sqrt{\frac{1-\gamma_a}{\gamma_a}} \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}} \leqslant 0$$

即使 $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} \leqslant 0$，如果：
- $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L}$ 接近0（均值项很小）
- $\sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}$ 很大（方差项很大）
- $\sqrt{\frac{1-\gamma_a}{\gamma_a}}$ 很大（$\gamma_a$ 很小）

那么：
$$\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} + \sqrt{\frac{1-\gamma_a}{\gamma_a}} \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}} > 0$$

**约束不满足！**

### 4. 当 r=4, p=3 时的具体情况

当 r=4, p=3 时：
- $\mathbf{v}_{j,L} = -(\sum_{k \in P} \mathbf{z}_k + \mathbf{z}_j)$
- $\mathbf{v}_{j,U} = \mathbf{z}_j - \frac{5}{3} \times \sum_{k \in P} \mathbf{z}_k$

#### 下界约束分析

SOCP约束：
$$\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} + \sqrt{\frac{1-\gamma_a}{\gamma_a}} \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}} \leqslant 0$$

由于：
- $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} = -\boldsymbol{\mu}_w^\top (\sum_{k \in P} \mathbf{z}_k + \mathbf{z}_j)$
- 如果 $\boldsymbol{\mu}_w^\top (\sum_{k \in P} \mathbf{z}_k + \mathbf{z}_j)$ 很大（总工作量很大），那么 $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L}$ 是一个很大的负数
- 但如果 $\sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}$ 也很大（不确定性很大），而 $\sqrt{\frac{1-\gamma_a}{\gamma_a}}$ 很大（$\gamma_a$ 很小），那么不确定性项可能超过均值项的绝对值

**关键**：即使确定性约束总是满足，SOCP近似约束仍可能不满足，因为：
1. **SOCP近似约束是保守的**：它提供了一个更严格的上界
2. **不确定性项的影响**：即使均值项为负，如果方差项很大，不确定性项可能使总和为正

#### 上界约束分析

SOCP约束：
$$\boldsymbol{\mu}_w^\top \mathbf{v}_{j,U} + \sqrt{\frac{1-\gamma_b}{\gamma_b}} \sqrt{\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U}} \leqslant 0$$

由于：
- $\mathbf{v}_{j,U} = \mathbf{z}_j - \frac{5}{3} \times \sum_{k \in P} \mathbf{z}_k$
- 如果区域 $j$ 的工作量较小，而总工作量较大，$\mathbf{v}_{j,U}$ 的分量可能是负数
- $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,U}$ 可能是负数
- 但如果 $\sqrt{\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U}}$ 很大，而 $\sqrt{\frac{1-\gamma_b}{\gamma_b}}$ 很大，不确定性项可能使总和为正

### 5. SOCP近似约束 vs 精确约束

**重要区别**：

1. **精确约束**：如果确定性约束总是满足，那么精确约束也应该总是满足（如之前的分析）

2. **SOCP近似约束**：即使确定性约束总是满足，SOCP近似约束仍可能不满足，因为：
   - SOCP近似约束是**保守的**（提供了更严格的上界）
   - 它考虑了不确定性项的影响
   - 即使均值项满足，如果方差项很大，约束仍可能不满足

### 6. 为什么SOCP近似约束可能失败？

SOCP近似约束使用Cantelli不等式，它提供了一个**保守的上界**：

$$P \left\{ \widetilde{\mathbf{d}}_w^\top \mathbf{v}_{j,L} > 0 \right\} \leqslant \frac{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L} + (\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L})^2}$$

当 $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} \leqslant 0$ 时，这个上界是有效的。

但是，SOCP约束：
$$\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} + \sqrt{\frac{1-\gamma_a}{\gamma_a}} \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}} \leqslant 0$$

要求均值项必须足够负，以抵消不确定性项。

**即使确定性约束总是满足**（即对于任意分布，$\widetilde{\mathbf{d}}_w^\top \mathbf{v}_{j,L} \leqslant 0$ 总是成立），如果：
- 均值项 $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L}$ 接近0
- 方差项 $\sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}$ 很大
- 系数 $\sqrt{\frac{1-\gamma_a}{\gamma_a}}$ 很大（$\gamma_a$ 很小）

那么SOCP约束可能不满足。

## 总结

**关键理解**：

1. ✅ **确定性约束总是满足**：当 r=4, p=3 时，确定性约束总是满足
2. ✅ **精确约束也应该总是满足**：如果确定性约束总是满足，精确约束也应该总是满足
3. ⚠️ **SOCP近似约束可能不满足**：即使确定性约束总是满足，SOCP近似约束仍可能不满足，因为：
   - SOCP近似约束是**保守的**（提供了更严格的上界）
   - 它要求均值项必须足够负，以抵消不确定性项的影响
   - 即使均值项为负，如果方差项很大，不确定性项可能使总和为正

**因此，代码中验证SOCP近似约束可能失败，这是正常的，因为SOCP近似约束比精确约束更严格。**

## 建议

1. **理解SOCP近似约束的特性**：它是保守的，可能比精确约束更严格
2. **检查实际数值**：运行代码，查看 `meanTerm_vjL`、`sqrtVar_vjL`、`factorLower` 的值
3. **如果总是违反**：考虑使用精确方法（`useExactMethod = true`），或者调整参数（使用更小的 r 值或更大的 $\gamma_a$、$\gamma_b$）
