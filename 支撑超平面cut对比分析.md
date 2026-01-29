# Assignment-Dependent模型 vs 原始模型的支撑超平面Cut对比分析

## 一、数量对比

### 1. Cut数量（每个区域）
两种模型在每个区域的cut数量相同：
- **每个区域**：最多2个cut（下界约束1个 + 上界约束1个）
- **p个区域**：最多2p个cut
- **示例**（p=30）：最多60个cut

### 2. 每次迭代可能添加的cut数量
- **Assignment-dependent模型**：限制每次最多添加**10个cut**
- **原始模型**：限制每次最多添加**20个cut**

原因：assignment-dependent模型的每个cut更复杂，需要更严格的限制。

---

## 二、复杂度对比

### 1. 向量维度

#### Assignment-Dependent模型
- **v_{j,L} 和 v_{j,U} 的维度**：`N × p` 维
- **示例**（N=671, p=30）：**20,130维**
- **协方差矩阵维度**：`(N×p) × (N×p) = 20,130 × 20,130`
- **协方差矩阵元素数量**：405,216,900个元素
- **内存需求**：约3 GB

#### 原始模型
- **v_{j,L} 和 v_{j,U} 的维度**：`N` 维
- **示例**（N=671）：**671维**
- **协方差矩阵维度**：`N × N = 671 × 671`
- **协方差矩阵元素数量**：450,241个元素
- **内存需求**：约3.4 MB

**维度比例**：Assignment-dependent模型是原始模型的 **p倍**（30倍）

---

### 2. 每个Cut包含的变量数量

#### Assignment-Dependent模型
每个cut包含**所有区域的所有变量**：
- **变量数量**：`N × p = 671 × 30 = 20,130` 个变量
- **计算逻辑**：
  ```java
  // 对于区域 j：N个变量
  for (int i = 0; i < n; i++) {
      cutExpr.addTerm(grad_ij, x[i][j]);  // x[i][j]
  }
  
  // 对于其他区域 k ≠ j：(p-1) × N 个变量
  for (int k = 0; k < centers.size(); k++) {
      if (k == j) continue;
      for (int i = 0; i < n; i++) {
          cutExpr.addTerm(grad_ik, x[i][k]);  // x[i][k]
      }
  }
  ```
- **总变量数**：`N + (p-1) × N = N × p = 20,130`

#### 原始模型
每个cut**只包含区域j的变量**：
- **变量数量**：`N = 671` 个变量
- **计算逻辑**：
  ```java
  // 只对于区域 j：N个变量
  for (int i = 0; i < inst.getN(); i++) {
      cutExpr.addTerm(grad_i, x[i][j]);  // 只有 x[i][j]
  }
  
  // 对于其他区域 k ≠ j：梯度为0，不添加任何项
  // （论文明确说明：v_{j,L} 和 v_{j,U} 只依赖于区域 j 的决策变量）
  ```
- **总变量数**：`N = 671`

**变量数量比例**：Assignment-dependent模型是原始模型的 **p倍**（30倍）

---

### 3. 梯度计算复杂度

#### Assignment-Dependent模型
- **协方差矩阵-向量乘积**：`Σ_w * v_j`（20,130维向量）
- **计算复杂度**：O((N×p)²) = O(20,130²) ≈ O(405,000,000)
- **梯度计算**：
  - 区域j：需要访问`meanVector[idx_ij]`和`sigmaV[idx_ij]`（N×p维）
  - 其他区域k：需要访问`meanVector[idx_ik]`和`sigmaV[idx_ik]`（N×p维）
  - **总计算量**：N × p 次梯度计算

#### 原始模型
- **协方差矩阵-向量乘积**：`Σ * v_j`（671维向量）
- **计算复杂度**：O(N²) = O(671²) ≈ O(450,000)
- **梯度计算**：
  - 区域j：需要访问`meanVector[i]`和`sigmaV[i]`（N维）
  - 其他区域k：梯度为0，无需计算
  - **总计算量**：N 次梯度计算

**计算复杂度比例**：Assignment-dependent模型是原始模型的 **p²倍**（900倍）

---

### 4. Cut的稀疏性

#### Assignment-Dependent模型
- **稀疏性**：**密集**（dense）
- **原因**：每个cut包含所有N×p个变量，几乎所有变量都有非零系数
- **约束矩阵密度**：接近100%

#### 原始模型
- **稀疏性**：**稀疏**（sparse）
- **原因**：每个cut只包含区域j的N个变量，其他变量系数为0
- **约束矩阵密度**：约 `N / (N×p) = 1/p = 3.3%`

---

## 三、内存和计算资源对比

### 1. 协方差矩阵存储

| 模型类型 | 矩阵维度 | 元素数量 | 内存需求（double） |
|---------|---------|---------|------------------|
| Assignment-dependent | 20,130 × 20,130 | 405,216,900 | ~3.1 GB |
| 原始模型 | 671 × 671 | 450,241 | ~3.4 MB |

**内存比例**：Assignment-dependent模型是原始模型的 **~900倍**

### 2. 每次迭代的计算时间

假设添加相同数量的cut（例如10个）：

| 模型类型 | 每个cut的变量数 | 10个cut的总变量数 | 协方差计算复杂度 |
|---------|---------------|-----------------|----------------|
| Assignment-dependent | 20,130 | 201,300 | O(405,000,000) |
| 原始模型 | 671 | 6,710 | O(450,000) |

**计算时间比例**：Assignment-dependent模型预计是原始模型的 **数百倍**

---

## 四、代码实现差异

### 1. Assignment-Dependent模型的Cut生成

```java
// v_j是N×p维向量
int n = inst.getN();  // 671
int p = numRegionsForAssignmentDependent;  // 30
int nTimesP = n * p;  // 20,130

// 计算 Σ_w * v_j（N×p维）
double[] sigmaV = computeCovarianceMatrixVectorProduct(v_j);

// 对于区域 j：N个变量
for (int i = 0; i < n; i++) {
    int idx_ij = i * p + j;
    grad_ij = (coeff - 1.0) * meanVector[idx_ij] + 
              factor * (coeff - 1.0) * sigmaV[idx_ij] / sqrtVar;
    cutExpr.addTerm(grad_ij, x[i][j]);
}

// 对于其他区域 k ≠ j：(p-1) × N 个变量
for (int k = 0; k < centers.size(); k++) {
    if (k == j) continue;
    for (int i = 0; i < n; i++) {
        int idx_ik = i * p + k;
        grad_ik = coeff * meanVector[idx_ik] + 
                  factor * coeff * sigmaV[idx_ik] / sqrtVar;
        cutExpr.addTerm(grad_ik, x[i][k]);
    }
}
```

### 2. 原始模型的Cut生成

```java
// v_j是N维向量
int n = inst.getN();  // 671

// 计算 Σ * v_j（N维）
double[] sigmaV = new double[inst.getN()];
for (int i = 0; i < inst.getN(); i++) {
    for (int k = 0; k < inst.getN(); k++) {
        sigmaV[i] += covarianceMatrix[i][k] * v_j[k];
    }
}

// 只对于区域 j：N个变量
for (int i = 0; i < inst.getN(); i++) {
    grad_i = -meanVector[i] + factor * (-sigmaV[i]) / sqrtVar;
    cutExpr.addTerm(grad_i, x[i][j]);
}

// 对于其他区域 k ≠ j：梯度为0，不添加任何项
// （论文明确说明：v_{j,L} 和 v_{j,U} 只依赖于区域 j 的决策变量）
```

---

## 五、总结

### 关键差异

| 维度 | Assignment-Dependent模型 | 原始模型 | 比例 |
|-----|------------------------|---------|------|
| **向量维度** | N×p (20,130) | N (671) | **30倍** |
| **协方差矩阵维度** | (N×p)×(N×p) | N×N | **900倍** |
| **每个cut的变量数** | N×p (20,130) | N (671) | **30倍** |
| **计算复杂度** | O((N×p)²) | O(N²) | **900倍** |
| **内存需求** | ~3.1 GB | ~3.4 MB | **~900倍** |
| **Cut稀疏性** | 密集（~100%） | 稀疏（~3.3%） | - |
| **每次迭代cut限制** | 10个 | 20个 | - |

### 为什么Assignment-Dependent模型更复杂？

1. **高维向量空间**：v_j从N维扩展到N×p维，导致所有相关计算都按p倍增长
2. **全区域耦合**：每个cut涉及所有区域的变量，而原始模型只涉及单个区域
3. **密集约束矩阵**：所有变量都有非零系数，导致约束矩阵非常密集
4. **大规模协方差矩阵**：需要存储和计算(N×p)×(N×p)的协方差矩阵

### 优化建议

1. **限制cut数量**：已实现（每次最多10个）
2. **优先添加最违反的cut**：已实现（按违反程度排序）
3. **时间限制**：已实现（每次迭代独立时间限制）
4. **可能的进一步优化**：
   - 考虑使用稀疏矩阵存储协方差矩阵
   - 考虑延迟计算（lazy evaluation）某些项
   - 考虑使用近似方法计算协方差矩阵-向量乘积
