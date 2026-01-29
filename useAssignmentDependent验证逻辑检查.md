# useAssignmentDependent 情况下验证逻辑检查

## 前提条件

- `centers.size() == 3`
- `numRegionsForAssignmentDependent == 3`
- 维度匹配，没有问题

## 检查验证逻辑

### 1. z_k 的构建逻辑

**代码（第4039-4052行）：**
```java
double[][] z_k = new double[p_assignment][nTimesP];
for (int k = 0; k < p_assignment; k++) {
    for (int i = 0; i < n; i++) {
        for (int j2 = 0; j2 < p_assignment; j2++) {
            int idx = i * p_assignment + j2;
            if (j2 == k) {
                z_k[k][idx] = xVal[i][k];
            } else {
                z_k[k][idx] = 0.0;
            }
        }
    }
}
```

**逻辑分析：**
- `z_k[k]` 是一个 `nTimesP` 维向量
- 对于每个 `i` 和 `j2`：
  - 如果 `j2 == k`：`z_k[k][i * p_assignment + j2] = xVal[i][k]`
  - 否则：`z_k[k][i * p_assignment + j2] = 0`

**结果：**
- `z_k[k]` 在位置 `i * p_assignment + k` 处的值为 `xVal[i][k]`
- 其他位置为0

**✅ 逻辑正确**

### 2. sum_z_k 的计算

**代码（第4056-4062行）：**
```java
double[] sum_z_k = new double[nTimesP];
for (int idx = 0; idx < nTimesP; idx++) {
    sum_z_k[idx] = 0.0;
    for (int k = 0; k < p_assignment; k++) {
        sum_z_k[idx] += z_k[k][idx];
    }
}
```

**逻辑分析：**
- 对于 `idx = i * p_assignment + j2`：
  - `sum_z_k[idx] = sum_k z_k[k][idx]`
  - 只有当 `k == j2` 时，`z_k[k][idx]` 才非零
  - 所以 `sum_z_k[idx] = z_k[j2][idx] = xVal[i][j2]`

**结果：**
- `sum_z_k[i * p_assignment + j2] = xVal[i][j2]`

**✅ 逻辑正确**

### 3. v_{j,L} 和 v_{j,U} 的计算

**代码（第4074-4077行）：**
```java
for (int idx = 0; idx < nTimesP; idx++) {
    v_jL[idx] = coeff * sum_z_k[idx] - z_j[idx];
    v_jU[idx] = z_j[idx] - coeffUpper * sum_z_k[idx];
}
```

**逻辑分析：**
- 对于 `idx = i * p_assignment + j2`：
  - `sum_z_k[idx] = xVal[i][j2]`
  - `z_j[idx] = z_k[j][idx] = xVal[i][j]`（当 `j2 == j` 时），否则为0
  - 所以：
    - `v_jL[idx] = coeff * xVal[i][j2] - (j2 == j ? xVal[i][j] : 0)`
    - `v_jU[idx] = (j2 == j ? xVal[i][j] : 0) - coeffUpper * xVal[i][j2]`

**关键问题：**
- 当 `j2 == j` 时：`v_jL[idx] = coeff * xVal[i][j] - xVal[i][j] = (coeff - 1) * xVal[i][j]`
- 当 `j2 != j` 时：`v_jL[idx] = coeff * xVal[i][j2]`

**数学上应该是什么？**

根据约束定义：
- `v_{j,L} = \frac{1-r}{p} \sum_{k \in P} \mathbf{z}_k - \mathbf{z}_j`
- `v_{j,U} = \mathbf{z}_j - \frac{1+r}{p} \sum_{k \in P} \mathbf{z}_k`

对于 `idx = i * p_assignment + j2`：
- `sum_z_k[idx] = xVal[i][j2]`（所有区域的总和）
- `z_j[idx] = xVal[i][j]`（当 `j2 == j` 时），否则为0

所以：
- `v_jL[idx] = coeff * xVal[i][j2] - (j2 == j ? xVal[i][j] : 0)`
- `v_jU[idx] = (j2 == j ? xVal[i][j] : 0) - coeffUpper * xVal[i][j2]`

**✅ 逻辑正确**

### 4. 均值项的计算

**代码（第4078-4083行）：**
```java
meanTerm_vjL = 0.0;
meanTerm_vjU = 0.0;
for (int idx = 0; idx < nTimesP; idx++) {
    meanTerm_vjL += meanVector[idx] * v_jL[idx];
    meanTerm_vjU += meanVector[idx] * v_jU[idx];
}
```

**逻辑分析：**
- `meanVector[idx]` 对应 `d_{ij}` 的均值，其中 `idx = i * p_assignment + j2`
- `v_jL[idx]` 和 `v_jU[idx]` 的计算如上

**✅ 逻辑正确**

### 5. 方差项的计算

**代码（第4085-4086行）：**
```java
varTerm_vjL = computeQuadraticForm(v_jL);
varTerm_vjU = computeQuadraticForm(v_jU);
```

**逻辑分析：**
- `computeQuadraticForm` 计算 `v^T * Σ * v`
- 对于 `useAssignmentDependent`，使用 `nTimesPForCovariance = n * numRegionsForAssignmentDependent`
- 如果 `useDiagonalCovariance = true`，只计算对角线项

**✅ 逻辑正确**

## 可能的问题

### 1. **coeff 的计算**

`coeff` 和 `coeffUpper` 是基于 `centers.size()` 计算的（第2471-2473行）：
```java
int p = centers.size();
double coeff = (1.0 - r) / p;
double coeffUpper = (1.0 + r) / p;
```

如果 `centers.size() == numRegionsForAssignmentDependent`，这没问题。

但如果将来这两个值不一致，可能会有问题。

### 2. **sum_z_k 的含义**

`sum_z_k[idx]` 表示所有区域的总和，其中 `idx = i * p_assignment + j2`。

但根据约束的定义：
- `\sum_{k \in P} \mathbf{z}_k` 应该是一个向量，其中每个位置对应所有区域的总和
- 对于 `idx = i * p_assignment + j2`，`sum_z_k[idx] = xVal[i][j2]`

**问题**：`sum_z_k[idx]` 的值是 `xVal[i][j2]`，但这似乎不对！

**应该是什么？**

根据约束定义：
- `\sum_{k \in P} \mathbf{z}_k` 在位置对应 `d_{ij}` 处的值应该是 `\sum_{k \in P} x_{ik}`（所有区域的总和）
- 但代码中 `sum_z_k[i * p_assignment + j2] = xVal[i][j2]`，这似乎只考虑了区域 `j2`

**让我重新理解：**

在 assignment-dependent 模型中：
- `\mathbf{z}_k` 是一个 `N \times p` 维向量
- `\mathbf{z}_k[i * p + j] = x_{ik}` 当 `j == k`，否则为0
- `\sum_{k \in P} \mathbf{z}_k[i * p + j] = \sum_{k \in P} (j == k ? x_{ik} : 0) = x_{ij}`（当 `j == k` 时）

**等等，这里有问题！**

`\sum_{k \in P} \mathbf{z}_k[i * p + j]` 应该等于 `x_{ij}`（当 `j == k` 时），但代码中：
- `sum_z_k[i * p_assignment + j2] = xVal[i][j2]`

这看起来是对的，因为：
- `z_k[k][i * p_assignment + j2] = xVal[i][k]`（当 `j2 == k` 时）
- `sum_z_k[i * p_assignment + j2] = sum_k z_k[k][i * p_assignment + j2] = z_k[j2][i * p_assignment + j2] = xVal[i][j2]`

**但是**，根据约束的定义，`\sum_{k \in P} \mathbf{z}_k` 在位置 `i * p + j` 处的值应该是所有区域的总和，而不是单个区域的值。

让我重新检查约束的定义...
