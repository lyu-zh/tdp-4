# checkAndAddRelativeBalanceCuts 方法中 useAssignmentDependent = true 情况详解

本文档详细讲解 `checkAndAddRelativeBalanceCuts` 方法在处理 `useAssignmentDependent = true` 情况下的代码实现，包括每个计算步骤对应的数学公式。

## 方法概述

`checkAndAddRelativeBalanceCuts` 方法用于检查相对平衡性约束是否被违反，如果违反则生成支撑超平面cut。该方法位于 `DistributionallyRobustAlgo.java` 文件的第 3987 行。

## 一、获取当前解的值

**代码位置：** 第 3990-3996 行

```java
double[][] currentXVal = new double[inst.getN()][centers.size()];
for (int i = 0; i < inst.getN(); i++) {
    for (int k = 0; k < centers.size(); k++) {
        currentXVal[i][k] = x[i][k].get(GRB.DoubleAttr.X);
    }
}
```

**说明：** 获取当前解的值 $\bar{\mathbf{x}}^k$，其中：
- $x_{ik}$ 表示基本单元 $i$ 是否分配给区域 $k$
- $\bar{x}_{ik}^k$ 表示在当前解 $\bar{\mathbf{x}}^k$ 下变量 $x_{ik}$ 的值

## 二、遍历相对平衡性约束信息

**代码位置：** 第 4001-4007 行

```java
for (RelativeBalanceConstraintInfo info : relativeBalanceConstraintInfos) {
    int j = info.j;
    int p = centers.size();
    double coeff = info.coeff;
    double coeffUpper = info.coeffUpper;
    double factorLower = info.factorLower;
    double factorUpper = info.factorUpper;
```

**对应公式：**
- $j$：当前区域索引
- $p$：区域总数
- $\text{coeff} = \frac{1-\alpha}{p}$
- $\text{coeffUpper} = \frac{1+\alpha}{p}$
- $\text{factorLower} = t_L$：下界约束的系数（根据 D1 或 D2 模糊集计算）
- $\text{factorUpper} = t_U$：上界约束的系数（根据 D1 或 D2 模糊集计算）

## 三、Assignment-Dependent 模型的核心计算

**代码位置：** 第 4020-4077 行

### 3.1 初始化向量维度

**代码位置：** 第 4021-4027 行

```java
if (useAssignmentDependent) {
    int n = inst.getN();
    int p = numRegionsForAssignmentDependent;
    int nTimesP = n * p;
    v_jL = new double[nTimesP];
    v_jU = new double[nTimesP];
```

**说明：**
- $n$：基本单元数量
- $p$：区域数量
- $n \times p$：$\mathbf{v}_{j,L}$ 和 $\mathbf{v}_{j,U}$ 的维度

### 3.2 构建 $\mathbf{z}_k$ 向量

**代码位置：** 第 4029-4046 行

```java
double[][] z_k = new double[p][nTimesP];
for (int k = 0; k < p; k++) {
    for (int i = 0; i < n; i++) {
        for (int j2 = 0; j2 < p; j2++) {
            int idx = i * p + j2;
            if (j2 == k) {
                z_k[k][idx] = xVal[i][k];
            } else {
                z_k[k][idx] = 0.0;
            }
        }
    }
}
```

**对应公式：**

$\mathbf{z}_k$ 是一个 $N \times p$ 维向量，其中：
- 向量化顺序：$[d_{11}, d_{12}, \ldots, d_{1p}, d_{21}, \ldots, d_{Np}]$
- 位置索引：$\text{idx} = i \times p + j_2$
- $\mathbf{z}_k[\text{idx}] = x_{ik}$ 当 $j_2 = k$，否则为 $0$

**数学表示：**
$$\mathbf{z}_k[i \times p + j_2] = \begin{cases}
x_{ik} & \text{if } j_2 = k \\
0 & \text{otherwise}
\end{cases}$$

### 3.3 计算 $\sum_{k \in P} \mathbf{z}_k$

**代码位置：** 第 4048-4055 行

```java
double[] sum_z_k = new double[nTimesP];
for (int idx = 0; idx < nTimesP; idx++) {
    sum_z_k[idx] = 0.0;
    for (int k = 0; k < centers.size(); k++) {
        sum_z_k[idx] += z_k[k][idx];
    }
}
```

**对应公式：**
$$\sum_{k \in P} \mathbf{z}_k = \sum_{k=0}^{p-1} \mathbf{z}_k$$

### 3.4 计算 $\mathbf{v}_{j,L}$ 和 $\mathbf{v}_{j,U}$

**代码位置：** 第 4057-4065 行

```java
double[] z_j = z_k[j];

for (int idx = 0; idx < nTimesP; idx++) {
    v_jL[idx] = coeff * sum_z_k[idx] - z_j[idx];
    v_jU[idx] = z_j[idx] - coeffUpper * sum_z_k[idx];
}
```

**对应公式：**

**下界约束：**
$$\mathbf{v}_{j,L} = \frac{1-\alpha}{p} \sum_{k \in P} \mathbf{z}_k - \mathbf{z}_j$$

**上界约束：**
$$\mathbf{v}_{j,U} = \mathbf{z}_j - \frac{1+\alpha}{p} \sum_{k \in P} \mathbf{z}_k$$

### 3.5 计算均值项 $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,L}$ 和 $\boldsymbol{\mu}_w^\top \mathbf{v}_{j,U}$

**代码位置：** 第 4067-4073 行

```java
meanTerm_vjL = 0.0;
meanTerm_vjU = 0.0;
for (int idx = 0; idx < nTimesP; idx++) {
    meanTerm_vjL += meanVector[idx] * v_jL[idx];
    meanTerm_vjU += meanVector[idx] * v_jU[idx];
}
```

**对应公式：**

**下界约束：**
$$\text{meanTerm\_vjL} = \boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} = \sum_{\text{idx}=0}^{n \times p - 1} \mu_{w,\text{idx}} \cdot v_{j,L}[\text{idx}]$$

**上界约束：**
$$\text{meanTerm\_vjU} = \boldsymbol{\mu}_w^\top \mathbf{v}_{j,U} = \sum_{\text{idx}=0}^{n \times p - 1} \mu_{w,\text{idx}} \cdot v_{j,U}[\text{idx}]$$

### 3.6 计算方差项 $\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}$ 和 $\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U}$

**代码位置：** 第 4075-4077 行

```java
varTerm_vjL = computeQuadraticForm(v_jL);
varTerm_vjU = computeQuadraticForm(v_jU);
```

**对应公式：**

**下界约束：**
$$\text{varTerm\_vjL} = \mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}$$

**上界约束：**
$$\text{varTerm\_vjU} = \mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U}$$

## 四、计算 $\varphi_L$ 和 $\varphi_U$

**代码位置：** 第 4147-4152 行

```java
double sqrtVar_vjL = Math.sqrt(Math.max(0.0, varTerm_vjL));
double sqrtVar_vjU = Math.sqrt(Math.max(0.0, varTerm_vjU));

double phi_L = meanTerm_vjL + factorLower * sqrtVar_vjL;
double phi_U = meanTerm_vjU + factorUpper * sqrtVar_vjU;
```

**对应公式：**

**下界约束：**
$$\varphi_L(\bar{\mathbf{x}}^k) = \boldsymbol{\mu}_w^\top \mathbf{v}_{j,L} + t_L \sqrt{\mathbf{v}_{j,L}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,L}}$$

其中：
- $t_L = \sqrt{\delta_1} + \sqrt{\frac{1 - \gamma_a}{\gamma_a} (\delta_2 - \delta_1)}$（Case 1，D2 模糊集）
- 或 $t_L = \sqrt{\frac{1-\gamma_a}{\gamma_a}}$（D1 模糊集）

**上界约束：**
$$\varphi_U(\bar{\mathbf{x}}^k) = \boldsymbol{\mu}_w^\top \mathbf{v}_{j,U} + t_U \sqrt{\mathbf{v}_{j,U}^\top \boldsymbol{\Sigma}_w \mathbf{v}_{j,U}}$$

其中：
- $t_U = \sqrt{\delta_1} + \sqrt{\frac{1 - \gamma_b}{\gamma_b} (\delta_2 - \delta_1)}$（Case 1，D2 模糊集）
- 或 $t_U = \sqrt{\frac{1-\gamma_b}{\gamma_b}}$（D1 模糊集）

## 五、检查约束违反并记录信息

**代码位置：** 第 4158-4249 行

如果 $\varphi_L > \text{violationTolerance}$ 或 $\varphi_U > \text{violationTolerance}$，则记录违反信息到 `ViolationInfo` 对象中，包括：
- `j`：区域索引
- `isLower`：是否为下界约束
- `violation`：违反值（$\varphi_L$ 或 $\varphi_U$）
- `xVal`：当前解的值（深拷贝）
- `v_j`：$\mathbf{v}_{j,L}$ 或 $\mathbf{v}_{j,U}$（深拷贝）
- `meanTerm`：均值项
- `sqrtVar`：标准差项
- `factor`：系数 $t_L$ 或 $t_U$
- `coeff`：$\frac{1-\alpha}{p}$
- `coeffUpper`：$\frac{1+\alpha}{p}$

## 六、生成支撑超平面 Cut

**代码位置：** 第 4262-4266 行

```java
for (ViolationInfo info : violations) {
    addRelativeBalanceSupportingCut(model, x, info.j, info.isLower, info.xVal, info.v_j, 
        info.meanTerm, info.sqrtVar, info.factor, info.coeff, info.coeffUpper);
    cutsAdded = true;
}
```

调用 `addRelativeBalanceSupportingCut` 方法生成支撑超平面cut。

## 七、addRelativeBalanceSupportingCut 方法详解

**代码位置：** 第 4289-4482 行

### 7.1 计算 $\varphi(\bar{\mathbf{x}}^k)$

**代码位置：** 第 4292-4293 行

```java
double phi = meanTerm + factor * sqrtVar;
```

**对应公式：**
$$\varphi(\bar{\mathbf{x}}^k) = \text{meanTerm} + \text{factor} \times \text{sqrtVar}$$

### 7.2 初始化 Cut 表达式

**代码位置：** 第 4298-4299 行

```java
GRBLinExpr cutExpr = new GRBLinExpr();
double constantTerm = phi; // 初始化为 φ(x̄ᵏ)
```

**对应公式：**
$$\text{constantTerm} = \varphi(\bar{\mathbf{x}}^k)$$

### 7.3 计算 $\boldsymbol{\Sigma}_w \mathbf{v}_j$

**代码位置：** 第 4307-4308 行

```java
double[] sigmaV = computeCovarianceMatrixVectorProduct(v_j);
```

**对应公式：**
$$\boldsymbol{\Sigma}_w \mathbf{v}_j = \text{computeCovarianceMatrixVectorProduct}(\mathbf{v}_j)$$

### 7.4 计算区域 $j$ 的梯度

**代码位置：** 第 4310-4334 行

#### 7.4.1 下界约束的梯度

**代码位置：** 第 4316-4321 行

```java
if (isLower) {
    double coeff_grad = coeff - 1.0; // (1-α)/p - 1
    grad_ij = coeff_grad * meanVector[idx_ij] + 
              factor * coeff_grad * sigmaV[idx_ij] / sqrtVar;
}
```

**对应公式：**
$$\frac{\partial \varphi_L}{\partial x_{ij}} = \left(\frac{1-\alpha}{p} - 1\right) \mu_{w,ij} + t_L \cdot \frac{\left(\frac{1-\alpha}{p} - 1\right) (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ij}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

其中：
- $\text{idx\_ij} = i \times p + j$：对应 $d_{ij}$ 在向量中的位置
- $\mu_{w,ij}$：$\boldsymbol{\mu}_w$ 在位置 $\text{idx\_ij}$ 的条目
- $(\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ij}$：$\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k$ 在位置 $\text{idx\_ij}$ 的条目

#### 7.4.2 上界约束的梯度

**代码位置：** 第 4322-4328 行

```java
else {
    double coeff_grad = 1.0 - coeffUpper; // 1 - (1+α)/p
    grad_ij = coeff_grad * meanVector[idx_ij] + 
              factor * coeff_grad * sigmaV[idx_ij] / sqrtVar;
}
```

**对应公式：**
$$\frac{\partial \varphi_U}{\partial x_{ij}} = \left(1 - \frac{1+\alpha}{p}\right) \mu_{w,ij} + t_U \cdot \frac{\left(1 - \frac{1+\alpha}{p}\right) (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k)_{ij}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

#### 7.4.3 添加变量项和更新常数项

**代码位置：** 第 4330-4333 行

```java
cutExpr.addTerm(grad_ij, x[i][j]);
constantTerm -= grad_ij * xVal[i][j];
```

**对应公式：**
- 添加变量项：$\frac{\partial \varphi}{\partial x_{ij}} \cdot x_{ij}$
- 更新常数项：$\text{constantTerm} = \text{constantTerm} - \frac{\partial \varphi}{\partial x_{ij}} \cdot \bar{x}_{ij}^k$

### 7.5 计算其他区域 $k \neq j$ 的梯度

**代码位置：** 第 4336-4362 行

#### 7.5.1 下界约束的梯度

**代码位置：** 第 4344-4349 行

```java
if (isLower) {
    grad_ik = coeff * meanVector[idx_ik] + 
              factor * coeff * sigmaV[idx_ik] / sqrtVar;
}
```

**对应公式：**
$$\frac{\partial \varphi_L}{\partial x_{ik}} = \frac{1-\alpha}{p} \mu_{w,ik} + t_L \cdot \frac{\frac{1-\alpha}{p} (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k)_{ik}}{\sqrt{(\bar{\mathbf{v}}_{j,L}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,L}^k}}$$

其中：
- $\text{idx\_ik} = i \times p + k$：对应 $d_{ik}$ 在向量中的位置

#### 7.5.2 上界约束的梯度

**代码位置：** 第 4350-4355 行

```java
else {
    grad_ik = -coeffUpper * meanVector[idx_ik] - 
              factor * coeffUpper * sigmaV[idx_ik] / sqrtVar;
}
```

**对应公式：**
$$\frac{\partial \varphi_U}{\partial x_{ik}} = -\frac{1+\alpha}{p} \mu_{w,ik} - t_U \cdot \frac{\frac{1+\alpha}{p} (\boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k)_{ik}}{\sqrt{(\bar{\mathbf{v}}_{j,U}^k)^\top \boldsymbol{\Sigma}_w \bar{\mathbf{v}}_{j,U}^k}}$$

#### 7.5.3 添加变量项和更新常数项

**代码位置：** 第 4357-4360 行

```java
cutExpr.addTerm(grad_ik, x[i][k]);
constantTerm -= grad_ik * xVal[i][k];
```

**对应公式：**
- 添加变量项：$\frac{\partial \varphi}{\partial x_{ik}} \cdot x_{ik}$
- 更新常数项：$\text{constantTerm} = \text{constantTerm} - \frac{\partial \varphi}{\partial x_{ik}} \cdot \bar{x}_{ik}^k$

### 7.6 最终 Cut 约束

**代码位置：** 第 4473-4481 行

```java
cutExpr.addConstant(constantTerm);

String cutName = isLower ? "rel_balance_cut_L_" + j + "_" + System.currentTimeMillis() : 
                               "rel_balance_cut_U_" + j + "_" + System.currentTimeMillis();
model.addConstr(cutExpr, GRB.LESS_EQUAL, 0, cutName);
```

**对应公式：**

最终的支撑超平面cut为：
$$\varphi(\bar{\mathbf{x}}^k) + \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top (\mathbf{x}_j - \bar{\mathbf{x}}_j^k) + \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top (\mathbf{x}_k - \bar{\mathbf{x}}_k^k) \leqslant 0$$

展开后：
$$\varphi(\bar{\mathbf{x}}^k) - \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top \bar{\mathbf{x}}_j^k - \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top \bar{\mathbf{x}}_k^k + \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top \mathbf{x}_j + \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top \mathbf{x}_k \leqslant 0$$

其中：
- 常数项：$\text{constantTerm} = \varphi(\bar{\mathbf{x}}^k) - \left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top \bar{\mathbf{x}}_j^k - \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top \bar{\mathbf{x}}_k^k$
- 变量项：$\left(\frac{\partial \varphi}{\partial \mathbf{x}_j}\right)^\top \mathbf{x}_j + \sum_{k \in P, k \neq j} \left(\frac{\partial \varphi}{\partial \mathbf{x}_k}\right)^\top \mathbf{x}_k$

## 总结

本文档详细说明了 `checkAndAddRelativeBalanceCuts` 方法在处理 `useAssignmentDependent = true` 情况下的实现，包括：

1. **向量构建**：如何从决策变量构建 $\mathbf{z}_k$ 向量
2. **约束向量计算**：如何计算 $\mathbf{v}_{j,L}$ 和 $\mathbf{v}_{j,U}$
3. **约束值计算**：如何计算 $\varphi_L$ 和 $\varphi_U$
4. **梯度计算**：如何计算所有区域的梯度
5. **Cut 生成**：如何构建支撑超平面cut

所有计算步骤都与论文中的数学公式一一对应，确保了实现的正确性。
