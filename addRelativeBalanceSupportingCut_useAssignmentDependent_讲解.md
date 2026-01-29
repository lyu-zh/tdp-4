# addRelativeBalanceSupportingCut 方法中 useAssignmentDependent = true 情况详解

本文档详细讲解 `addRelativeBalanceSupportingCut` 方法在处理 `useAssignmentDependent = true` 情况下的代码实现，包括每个计算步骤对应的数学公式。

## 方法概述

`addRelativeBalanceSupportingCut` 方法用于生成相对平衡性约束的支撑超平面cut。该方法位于 `DistributionallyRobustAlgo.java` 文件的第 4289 行。

**方法签名：**
```java
private void addRelativeBalanceSupportingCut(GRBModel model, GRBVar[][] x, int j, boolean isLower,
double[][] xVal, double[] v_j, double meanTerm, double sqrtVar,
double factor, double coeff, double coeffUpper) throws GRBException
```

## 一、方法参数说明

**代码位置：** 第 4289-4291 行

| 参数 | 类型 | 说明 | 对应数学符号 |
|------|------|------|-------------|
| `model` | `GRBModel` | Gurobi 优化模型 | - |
| `x` | `GRBVar[][]` | 决策变量矩阵 | $\mathbf{x}$ |
| `j` | `int` | 区域索引 | $j$ |
| `isLower` | `boolean` | 是否为下界约束 | - |
| `xVal` | `double[][]` | 当前解的值 | $\bar{\mathbf{x}}^k$ |
| `v_j` | `double[]` | 当前解下的 $\mathbf{v}_{j,L}$ 或 $\mathbf{v}_{j,U}$ 向量 | $\bar{\mathbf{v}}_{j,L}^k$ 或 $\bar{\mathbf{v}}_{j,U}^k$ |
| `meanTerm` | `double` | 均值项 | $\boldsymbol{\mu}_w^\top \mathbf{v}_j$ |
| `sqrtVar` | `double` | 标准差项 | $\sqrt{\mathbf{v}_j^\top \boldsymbol{\Sigma}_w \mathbf{v}_j}$ |
| `factor` | `double` | 系数 | $t_L$ 或 $t_U$ |
| `coeff` | `double` | 系数 | $\frac{1-\alpha}{p}$ |
| `coeffUpper` | `double` | 系数 | $\frac{1+\alpha}{p}$ |

## 二、计算 $\varphi(\bar{\mathbf{x}}^k)$

**代码位置：** 第 4292-4293 行

```java
// 计算 φ 的值
double phi = meanTerm + factor * sqrtVar;
```

**对应公式：**

**下界约束（isLower = true）：**
$$\varphi_L(\bar{\mathbf{x}}^k) = \boldsymbol{\mu}_w^\top \bar{\mathbf{v}}_{j,L}^k + t_L \sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}$$

**上界约束（isLower = false）：**
$$\varphi_U(\bar{\mathbf{x}}^k) = \boldsymbol{\mu}_w^\top \bar{\mathbf{v}}_{j,U}^k + t_U \sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}$$

其中：
- $\text{meanTerm} = \boldsymbol{\mu}_w^\top \mathbf{v}_j$
- $\text{sqrtVar} = \sqrt{\mathbf{v}_j^\top \boldsymbol{\Sigma}_w \mathbf{v}_j}$
- $\text{factor} = t_L$（下界）或 $t_U$（上界）

## 三、初始化 Cut 表达式

**代码位置：** 第 4295-4299 行

```java
// 构建cut表达式
// 公式：φ(x̄ᵏ) + (∂φ/∂x_j)ᵀ (x_j - x̄_jᵏ) + Σ_{k≠j} (∂φ/∂x_k)ᵀ (x_k - x̄_kᵏ) ≤ 0
// 展开：φ(x̄ᵏ) - (∂φ/∂x_j)ᵀ x̄_jᵏ - Σ_{k≠j} (∂φ/∂x_k)ᵀ x̄_kᵏ + (∂φ/∂x_j)ᵀ x_j + Σ_{k≠j} (∂φ/∂x_k)ᵀ x_k ≤ 0
GRBLinExpr cutExpr = new GRBLinExpr();
double constantTerm = phi; // 初始化为 φ(x̄ᵏ)
```

**对应公式：**

支撑超平面cut的标准形式：
$$\varphi(\bar{\mathbf{x}}^k) + \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top (\mathbf{x}_j - \bar{\mathbf{x}}_j^k) + \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top (\mathbf{x}_k - \bar{\mathbf{x}}_k^k) \leqslant 0$$

展开后：
$$\varphi(\bar{\mathbf{x}}^k) - \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top \bar{\mathbf{x}}_j^k - \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top \bar{\mathbf{x}}_k^k + \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top \mathbf{x}_j + \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top \mathbf{x}_k \leqslant 0$$

**说明：**
- `cutExpr`：用于存储线性表达式的变量项
- `constantTerm`：初始化为 $\varphi(\bar{\mathbf{x}}^k)$，后续会减去梯度与当前解的乘积

## 四、Assignment-Dependent 模型的核心计算

**代码位置：** 第 4301-4362 行

### 4.1 初始化变量

**代码位置：** 第 4304-4305 行

```java
int n = inst.getN();
int p = numRegionsForAssignmentDependent;
```

**说明：**
- $n$：基本单元数量
- $p$：区域数量（在 assignment-dependent 模型中等于 `numRegionsForAssignmentDependent`）

### 4.2 计算 $\boldsymbol{\Sigma}_w \mathbf{v}_j$

**代码位置：** 第 4307-4308 行

```java
// 计算 Σ_w * v_j（N*p维）
double[] sigmaV = computeCovarianceMatrixVectorProduct(v_j);
```

**对应公式：**
$$\boldsymbol{\Sigma}_w \mathbf{v}_j = \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k \quad \text{或} \quad \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k$$

**说明：**
- `sigmaV` 是一个 $N \times p$ 维向量
- `sigmaV[idx]` 表示 $\boldsymbol{\Sigma}_w \mathbf{v}_j$ 在位置 `idx` 的条目
- 这个向量用于后续梯度计算中的方差项

## 五、计算区域 $j$ 的梯度

**代码位置：** 第 4310-4334 行

### 5.1 位置索引计算

**代码位置：** 第 4311-4313 行

```java
for (int i = 0; i < n; i++) {
    // 位置索引：对应 d_ij 在向量中的位置
    int idx_ij = i * p + j;
```

**对应公式：**
$$\text{idx\_ij} = i \times p + j$$

**说明：**
- 在 assignment-dependent 模型中，向量化顺序为：$[d_{11}, d_{12}, \ldots, d_{1p}, d_{21}, \ldots, d_{Np}]$
- `idx_ij` 对应 $d_{ij}$ 在向量中的位置

### 5.2 下界约束的梯度（区域 $j$）

**代码位置：** 第 4316-4321 行

```java
if (isLower) {
    // 下界约束：∂φ_L/∂x_ij = ((1-α)/p - 1) * μ_{w,ij} + factor * ((1-α)/p - 1) * (Σ_w * v_{j,L}^k)_{ij} / sqrt(...)
    // 其中 ∂v_{j,L}/∂x_ij = ((1-α)/p - 1) * e_ij
    double coeff_grad = coeff - 1.0; // (1-α)/p - 1
    grad_ij = coeff_grad * meanVector[idx_ij] + 
              factor * coeff_grad * sigmaV[idx_ij] / sqrtVar;
}
```

**对应公式：**

根据论文公式，对于下界约束：
$$\frac{\partial \varphi_L}{\partial x_{ij}} = \left(\frac{1-\alpha}{p} - 1\right) \mu_{w,ij} + t_L \cdot \frac{\left(\frac{1-\alpha}{p} - 1\right) (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ij}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

**推导过程：**

1. **$\mathbf{v}_{j,L}$ 对 $x_{ij}$ 的偏导数：**
   $$\frac{\partial \mathbf{v}_{j,L}}{\partial x_{ij}} = \frac{\partial}{\partial x_{ij}} \left(\frac{1-\alpha}{p}\sum_{k \in P} \mathbf{z}_k - \mathbf{z}_j\right) = \left(\frac{1-\alpha}{p} - 1\right) \mathbf{e}_{ij}$$

   其中 $\mathbf{e}_{ij}$ 是 $N \times p$ 维单位向量，在位置对应 $d_{ij}$ 处为 1，其他位置为 0。

2. **$\varphi_L$ 对 $x_{ij}$ 的偏导数：**
   $$\frac{\partial \varphi_L}{\partial x_{ij}} = \frac{\partial}{\partial x_{ij}} \left(\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} + t_L \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}\right)$$

   展开后：
   $$\frac{\partial \varphi_L}{\partial x_{ij}} = \boldsymbol{\mu}_w^\top \frac{\partial \mathbf{v}_{j,L}}{\partial x_{ij}} + t_L \cdot \frac{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \frac{\partial \mathbf{v}_{j,L}}{\partial x_{ij}}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

   代入 $\frac{\partial \mathbf{v}_{j,L}}{\partial x_{ij}} = \left(\frac{1-\alpha}{p} - 1\right) \mathbf{e}_{ij}$：
   $$\frac{\partial \varphi_L}{\partial x_{ij}} = \left(\frac{1-\alpha}{p} - 1\right) \mu_{w,ij} + t_L \cdot \frac{\left(\frac{1-\alpha}{p} - 1\right) (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ij}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

**代码实现：**
- `coeff_grad = coeff - 1.0 = (1-α)/p - 1`
- `meanVector[idx_ij] = μ_{w,ij}`
- `sigmaV[idx_ij] = (Σ_w * v_{j,L}^k)_{ij}`
- `sqrtVar = sqrt((v_{j,L}^k)^T * Σ_w * v_{j,L}^k)`

### 5.3 上界约束的梯度（区域 $j$）

**代码位置：** 第 4322-4328 行

```java
else {
    // 上界约束：∂φ_U/∂x_ij = (1 - (1+α)/p) * μ_{w,ij} + factor * (1 - (1+α)/p) * (Σ_w * v_{j,U}^k)_{ij} / sqrt(...)
    // 其中 ∂v_{j,U}/∂x_ij = (1 - (1+α)/p) * e_ij
    double coeff_grad = 1.0 - coeffUpper; // 1 - (1+α)/p
    grad_ij = coeff_grad * meanVector[idx_ij] + 
              factor * coeff_grad * sigmaV[idx_ij] / sqrtVar;
}
```

**对应公式：**

对于上界约束：
$$\frac{\partial \varphi_U}{\partial x_{ij}} = \left(1 - \frac{1+\alpha}{p}\right) \mu_{w,ij} + t_U \cdot \frac{\left(1 - \frac{1+\alpha}{p}\right) (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k)_{ij}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

**推导过程：**

1. **$\mathbf{v}_{j,U}$ 对 $x_{ij}$ 的偏导数：**
   $$\frac{\partial \mathbf{v}_{j,U}}{\partial x_{ij}} = \frac{\partial}{\partial x_{ij}} \left(\mathbf{z}_j - \frac{1+\alpha}{p}\sum_{k \in P} \mathbf{z}_k\right) = \left(1 - \frac{1+\alpha}{p}\right) \mathbf{e}_{ij}$$

2. **$\varphi_U$ 对 $x_{ij}$ 的偏导数：**
   类似下界约束的推导，得到：
   $$\frac{\partial \varphi_U}{\partial x_{ij}} = \left(1 - \frac{1+\alpha}{p}\right) \mu_{w,ij} + t_U \cdot \frac{\left(1 - \frac{1+\alpha}{p}\right) (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k)_{ij}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

**代码实现：**
- `coeff_grad = 1.0 - coeffUpper = 1 - (1+α)/p`
- `meanVector[idx_ij] = μ_{w,ij}`
- `sigmaV[idx_ij] = (Σ_w * v_{j,U}^k)_{ij}`
- `sqrtVar = sqrt((v_{j,U}^k)^T * Σ_w * v_{j,U}^k)`

### 5.4 添加变量项和更新常数项（区域 $j$）

**代码位置：** 第 4330-4333 行

```java
// 添加变量项：(∂φ/∂x_j)ᵀ x_j
cutExpr.addTerm(grad_ij, x[i][j]);
// 减去常数项：(∂φ/∂x_j)ᵀ x̄_jᵏ
constantTerm -= grad_ij * xVal[i][j];
```

**对应公式：**

- **添加变量项：** $\frac{\partial \varphi}{\partial x_{ij}} \cdot x_{ij}$
- **更新常数项：** $\text{constantTerm} = \text{constantTerm} - \frac{\partial \varphi}{\partial x_{ij}} \cdot \bar{x}_{ij}^k$

**说明：**
- 变量项会在 cut 约束中作为决策变量的系数
- 常数项需要减去梯度与当前解的乘积，以构建完整的线性表达式

## 六、计算其他区域 $k \neq j$ 的梯度

**代码位置：** 第 4336-4362 行

### 6.1 位置索引计算

**代码位置：** 第 4340-4342 行

```java
for (int k = 0; k < centers.size(); k++) {
    if (k == j) continue;
    
    for (int i = 0; i < n; i++) {
        // 位置索引：对应 d_ik 在向量中的位置
        int idx_ik = i * p + k;
```

**对应公式：**
$$\text{idx\_ik} = i \times p + k$$

**说明：**
- `idx_ik` 对应 $d_{ik}$ 在向量中的位置
- 遍历所有区域 $k \neq j$ 和所有基本单元 $i$

### 6.2 下界约束的梯度（其他区域 $k \neq j$）

**代码位置：** 第 4344-4349 行

```java
if (isLower) {
    // 下界约束：∂φ_L/∂x_ik = (1-α)/p * μ_{w,ik} + factor * (1-α)/p * (Σ_w * v_{j,L}^k)_{ik} / sqrt(...)
    // 其中 ∂v_{j,L}/∂x_ik = (1-α)/p * e_ik
    grad_ik = coeff * meanVector[idx_ik] + 
              factor * coeff * sigmaV[idx_ik] / sqrtVar;
}
```

**对应公式：**

对于下界约束，其他区域 $k \neq j$ 的梯度：
$$\frac{\partial \varphi_L}{\partial x_{ik}} = \frac{1-\alpha}{p} \mu_{w,ik} + t_L \cdot \frac{\frac{1-\alpha}{p} (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ik}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

**推导过程：**

1. **$\mathbf{v}_{j,L}$ 对 $x_{ik}$ 的偏导数（$k \neq j$）：**
   $$\frac{\partial \mathbf{v}_{j,L}}{\partial x_{ik}} = \frac{\partial}{\partial x_{ik}} \left(\frac{1-\alpha}{p}\sum_{k \in P} \mathbf{z}_k - \mathbf{z}_j\right) = \frac{1-\alpha}{p} \mathbf{e}_{ik}$$

   因为只有 $\mathbf{z}_k$ 依赖于 $x_{ik}$，而 $\mathbf{z}_j$ 不依赖于 $x_{ik}$（当 $k \neq j$ 时）。

2. **$\varphi_L$ 对 $x_{ik}$ 的偏导数：**
   $$\frac{\partial \varphi_L}{\partial x_{ik}} = \boldsymbol{\mu}_w^\top \frac{\partial \mathbf{v}_{j,L}}{\partial x_{ik}} + t_L \cdot \frac{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \frac{\partial \mathbf{v}_{j,L}}{\partial x_{ik}}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

   代入 $\frac{\partial \mathbf{v}_{j,L}}{\partial x_{ik}} = \frac{1-\alpha}{p} \mathbf{e}_{ik}$：
   $$\frac{\partial \varphi_L}{\partial x_{ik}} = \frac{1-\alpha}{p} \mu_{w,ik} + t_L \cdot \frac{\frac{1-\alpha}{p} (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ik}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

**代码实现：**
- `coeff = (1-α)/p`
- `meanVector[idx_ik] = μ_{w,ik}`
- `sigmaV[idx_ik] = (Σ_w * v_{j,L}^k)_{ik}`
- `sqrtVar = sqrt((v_{j,L}^k)^T * Σ_w * v_{j,L}^k)`

### 6.3 上界约束的梯度（其他区域 $k \neq j$）

**代码位置：** 第 4350-4355 行

```java
else {
    // 上界约束：∂φ_U/∂x_ik = -(1+α)/p * μ_{w,ik} - factor * (1+α)/p * (Σ_w * v_{j,U}^k)_{ik} / sqrt(...)
    // 其中 ∂v_{j,U}/∂x_ik = -(1+α)/p * e_ik
    grad_ik = -coeffUpper * meanVector[idx_ik] - 
              factor * coeffUpper * sigmaV[idx_ik] / sqrtVar;
}
```

**对应公式：**

对于上界约束，其他区域 $k \neq j$ 的梯度：
$$\frac{\partial \varphi_U}{\partial x_{ik}} = -\frac{1+\alpha}{p} \mu_{w,ik} - t_U \cdot \frac{\frac{1+\alpha}{p} (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k)_{ik}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

**推导过程：**

1. **$\mathbf{v}_{j,U}$ 对 $x_{ik}$ 的偏导数（$k \neq j$）：**
   $$\frac{\partial \mathbf{v}_{j,U}}{\partial x_{ik}} = \frac{\partial}{\partial x_{ik}} \left(\mathbf{z}_j - \frac{1+\alpha}{p}\sum_{k \in P} \mathbf{z}_k\right) = -\frac{1+\alpha}{p} \mathbf{e}_{ik}$$

2. **$\varphi_U$ 对 $x_{ik}$ 的偏导数：**
   $$\frac{\partial \varphi_U}{\partial x_{ik}} = \boldsymbol{\mu}_w^\top \frac{\partial \mathbf{v}_{j,U}}{\partial x_{ik}} + t_U \cdot \frac{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \frac{\partial \mathbf{v}_{j,U}}{\partial x_{ik}}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

   代入 $\frac{\partial \mathbf{v}_{j,U}}{\partial x_{ik}} = -\frac{1+\alpha}{p} \mathbf{e}_{ik}$：
   $$\frac{\partial \varphi_U}{\partial x_{ik}} = -\frac{1+\alpha}{p} \mu_{w,ik} - t_U \cdot \frac{\frac{1+\alpha}{p} (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k)_{ik}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

**代码实现：**
- `coeffUpper = (1+α)/p`
- `meanVector[idx_ik] = μ_{w,ik}`
- `sigmaV[idx_ik] = (Σ_w * v_{j,U}^k)_{ik}`
- `sqrtVar = sqrt((v_{j,U}^k)^T * Σ_w * v_{j,U}^k)`

### 6.4 添加变量项和更新常数项（其他区域 $k \neq j$）

**代码位置：** 第 4357-4360 行

```java
// 添加变量项：(∂φ/∂x_k)ᵀ x_k
cutExpr.addTerm(grad_ik, x[i][k]);
// 减去常数项：(∂φ/∂x_k)ᵀ x̄_kᵏ
constantTerm -= grad_ik * xVal[i][k];
```

**对应公式：**

- **添加变量项：** $\frac{\partial \varphi}{\partial x_{ik}} \cdot x_{ik}$
- **更新常数项：** $\text{constantTerm} = \text{constantTerm} - \frac{\partial \varphi}{\partial x_{ik}} \cdot \bar{x}_{ik}^k$

## 七、构建最终的 Cut 约束

**代码位置：** 第 4473-4481 行

```java
// 设置常数项
// 常数项 = φ(x̄ᵏ) - (∂φ/∂x_j)ᵀ x̄_jᵏ - Σ_{k≠j} (∂φ/∂x_k)ᵀ x̄_kᵏ
// 对于assignment-dependent模型，梯度涉及所有区域的变量，因此需要从所有区域减去梯度与当前解的乘积
cutExpr.addConstant(constantTerm);

// 添加cut约束
String cutName = isLower ? "rel_balance_cut_L_" + j + "_" + System.currentTimeMillis() : 
                               "rel_balance_cut_U_" + j + "_" + System.currentTimeMillis();
model.addConstr(cutExpr, GRB.LESS_EQUAL, 0, cutName);
```

**对应公式：**

最终的支撑超平面cut为：
$$\varphi(\bar{\mathbf{x}}^k) + \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top (\mathbf{x}_j - \bar{\mathbf{x}}_j^k) + \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top (\mathbf{x}_k - \bar{\mathbf{x}}_k^k) \leqslant 0$$

**展开形式：**

$$\underbrace{\varphi(\bar{\mathbf{x}}^k) - \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top \bar{\mathbf{x}}_j^k - \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top \bar{\mathbf{x}}_k^k}_{\text{常数项}} + \underbrace{\left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top \mathbf{x}_j + \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top \mathbf{x}_k}_{\text{变量项}} \leqslant 0$$

**详细展开：**

对于下界约束：
$$\begin{aligned}
&\varphi_L(\bar{\mathbf{x}}^k) - \sum_{i=1}^{n} \frac{\partial \varphi_L}{\partial x_{ij}} \bar{x}_{ij}^k - \sum_{k \in P, k \neq j} \sum_{i=1}^{n} \frac{\partial \varphi_L}{\partial x_{ik}} \bar{x}_{ik}^k \\
&\quad + \sum_{i=1}^{n} \frac{\partial \varphi_L}{\partial x_{ij}} x_{ij} + \sum_{k \in P, k \neq j} \sum_{i=1}^{n} \frac{\partial \varphi_L}{\partial x_{ik}} x_{ik} \leqslant 0
\end{aligned}$$

对于上界约束：
$$\begin{aligned}
&\varphi_U(\bar{\mathbf{x}}^k) - \sum_{i=1}^{n} \frac{\partial \varphi_U}{\partial x_{ij}} \bar{x}_{ij}^k - \sum_{k \in P, k \neq j} \sum_{i=1}^{n} \frac{\partial \varphi_U}{\partial x_{ik}} \bar{x}_{ik}^k \\
&\quad + \sum_{i=1}^{n} \frac{\partial \varphi_U}{\partial x_{ij}} x_{ij} + \sum_{k \in P, k \neq j} \sum_{i=1}^{n} \frac{\partial \varphi_U}{\partial x_{ik}} x_{ik} \leqslant 0
\end{aligned}$$

**说明：**
- **常数项**：$\text{constantTerm} = \varphi(\bar{\mathbf{x}}^k) - \sum_{i=1}^{n} \frac{\partial \varphi}{\partial x_{ij}} \bar{x}_{ij}^k - \sum_{k \in P, k \neq j} \sum_{i=1}^{n} \frac{\partial \varphi}{\partial x_{ik}} \bar{x}_{ik}^k$
- **变量项**：$\sum_{i=1}^{n} \frac{\partial \varphi}{\partial x_{ij}} x_{ij} + \sum_{k \in P, k \neq j} \sum_{i=1}^{n} \frac{\partial \varphi}{\partial x_{ik}} x_{ik}$

## 八、梯度计算公式总结

### 8.1 下界约束的梯度

**区域 $j$：**
$$\frac{\partial \varphi_L}{\partial x_{ij}} = \left(\frac{1-\alpha}{p} - 1\right) \mu_{w,ij} + t_L \cdot \frac{\left(\frac{1-\alpha}{p} - 1\right) (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ij}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

**其他区域 $k \neq j$：**
$$\frac{\partial \varphi_L}{\partial x_{ik}} = \frac{1-\alpha}{p} \mu_{w,ik} + t_L \cdot \frac{\frac{1-\alpha}{p} (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ik}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

### 8.2 上界约束的梯度

**区域 $j$：**
$$\frac{\partial \varphi_U}{\partial x_{ij}} = \left(1 - \frac{1+\alpha}{p}\right) \mu_{w,ij} + t_U \cdot \frac{\left(1 - \frac{1+\alpha}{p}\right) (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k)_{ij}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

**其他区域 $k \neq j$：**
$$\frac{\partial \varphi_U}{\partial x_{ik}} = -\frac{1+\alpha}{p} \mu_{w,ik} - t_U \cdot \frac{\frac{1+\alpha}{p} (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k)_{ik}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

## 九、关键点总结

1. **梯度计算涉及所有区域**：在 assignment-dependent 模型中，$\mathbf{v}_{j,L}$ 和 $\mathbf{v}_{j,U}$ 依赖于所有区域的决策变量，因此需要计算所有区域的梯度。

2. **位置索引映射**：向量化顺序为 $[d_{11}, d_{12}, \ldots, d_{1p}, d_{21}, \ldots, d_{Np}]$，位置索引为 $\text{idx} = i \times p + j$。

3. **梯度公式的对称性**：
   - 区域 $j$ 的梯度系数为 $\frac{1-\alpha}{p} - 1$（下界）或 $1 - \frac{1+\alpha}{p}$（上界）
   - 其他区域 $k \neq j$ 的梯度系数为 $\frac{1-\alpha}{p}$（下界）或 $-\frac{1+\alpha}{p}$（上界）

4. **常数项的计算**：需要从初始值 $\varphi(\bar{\mathbf{x}}^k)$ 中减去所有区域的梯度与当前解的乘积。

5. **Cut 的有效性**：支撑超平面 cut 是凸函数 $\varphi(\mathbf{x}) \leqslant 0$ 在点 $\bar{\mathbf{x}}^k$ 处的线性外逼近，对于所有 $\mathbf{x}$ 都有效。

## 总结

本文档详细说明了 `addRelativeBalanceSupportingCut` 方法在处理 `useAssignmentDependent = true` 情况下的实现，包括：

1. **参数说明**：所有输入参数的含义和对应的数学符号
2. **$\varphi$ 值的计算**：基于均值项和方差项
3. **梯度计算**：区域 $j$ 和其他区域 $k \neq j$ 的梯度公式及推导
4. **Cut 构建**：如何将梯度信息组合成最终的线性约束

所有计算步骤都与论文中的数学公式一一对应，确保了实现的正确性。
