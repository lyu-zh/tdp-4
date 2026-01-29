import com.gurobi.gurobi.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * TSP对偶值验证器
 * 对比 LP 对偶值 和 MIP 边际成本 (Cost_with_i - Cost_without_i)
 */
public class TSPDualValidator {
    
    // 复用参数
    private static final double EARTH_RADIUS_METERS = 6371000.0;
    private static GRBEnv env;
    // 同步锁：确保 Gurobi 操作的线程安全
    private static final Object gurobiLock = new Object();

    // ----------------------------------------------------------------
    // 基础数据结构与工具方法 (复用自 TSPDualExtractor)
    // ----------------------------------------------------------------
    
    static class Point {
        int id; double longitude; double latitude;
        Point(int id, double lon, double lat) { this.id = id; this.longitude = lon; this.latitude = lat; }
    }
    
    private static double haversineDistance(double lon1, double lat1, double lon2, double lat2) {
        double lat1Rad = Math.toRadians(lat1);
        double lat2Rad = Math.toRadians(lat2);
        double deltaLatRad = Math.toRadians(lat2 - lat1);
        double deltaLonRad = Math.toRadians(lon2 - lon1);
        double a = Math.sin(deltaLatRad/2)*Math.sin(deltaLatRad/2) +
                   Math.cos(lat1Rad)*Math.cos(lat2Rad)*Math.sin(deltaLonRad/2)*Math.sin(deltaLonRad/2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        return EARTH_RADIUS_METERS * c;
    }
    
    private static List<Point> readCoordinates(String filename) throws IOException {
        List<Point> points = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = reader.readLine();
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            String[] parts = line.split(",");
            if (parts.length >= 3) {
                try {
                    points.add(new Point(Integer.parseInt(parts[0].trim()), 
                                       Double.parseDouble(parts[1].trim()), 
                                       Double.parseDouble(parts[2].trim())));
                } catch (Exception e) {}
            }
        }
        reader.close();
        return points;
    }

    private static double[][] computeDistanceMatrix(List<Point> points) {
        int n = points.size();
        double[][] dist = new double[n][n];
        System.out.println("Calculating distance matrix...");
        IntStream.range(0, n).parallel().forEach(i -> {
            for (int j = i + 1; j < n; j++) {
                double d = haversineDistance(points.get(i).longitude, points.get(i).latitude,
                                           points.get(j).longitude, points.get(j).latitude);
                dist[i][j] = d; dist[j][i] = d;
            }
        });
        return dist;
    }

    private static int findPointIndexById(List<Point> points, int id) {
        for (int i=0; i<points.size(); i++) if(points.get(i).id == id) return i;
        return -1;
    }

    static class TSPPointSet {
        Set<Integer> pointSet; int extraPointIndex; int[][] xMatrix;
        TSPPointSet(Set<Integer> p, int e, int[][] x) { pointSet = p; extraPointIndex = e; xMatrix = x; }
    }

    /**
     * 预计算x_{i,j}数组：对于每个区域j，找到最近的n/k个点，这些点的x_{i,j}=1
     * 然后对于遍历的点i，如果x_{i,j}=0，则将其设为1
     * @return x[i][j] = 1表示点i被分配给区域j，否则为0
     */
    private static int[][] computeXMatrix(List<Point> allPoints, List<Point> centers, double[][] dist, int n, int k) {
        int[][] x = new int[n][k];
        
        // 对于每个区域j，计算最近的n/k个点
        for (int j = 0; j < k; j++) {
            int centerId = centers.get(j).id;
            // 直接计算索引：因为centerId就是原数据中的序号（从1开始），索引=id-1
            int centerIdx = centerId - 1;
            if (centerIdx < 0 || centerIdx >= allPoints.size()) continue;
            
            // 查找所有点到中心点的距离（包括中心点本身，距离为0）
            List<AbstractMap.SimpleEntry<Integer, Double>> dists = new ArrayList<>();
            for (int idx = 0; idx < n; idx++) {
                dists.add(new AbstractMap.SimpleEntry<>(idx, dist[centerIdx][idx]));
            }
            dists.sort(Map.Entry.comparingByValue());
            
            // 取最近的n/k个点，设置x_{i,j}=1（包括中心点本身，因为中心点到自己的距离为0，会排在第一位）
            int numBasePoints = n / k;
            for (int idx = 0; idx < numBasePoints && idx < dists.size(); idx++) {
                int pointIdx = dists.get(idx).getKey();
                x[pointIdx][j] = 1;
            }
        }
        
        return x;
    }

    /**
     * 确定TSP点集：基于x_{i,j}数组，对于区域j，包含所有x_{i,j}=1的点
     * 另外，对于遍历的点i，如果x_{i,j}=0，则将其设为1并加入点集
     */
    private static TSPPointSet determineTSPPointSet(int j, int i, List<Point> allPoints, List<Point> centers, 
                                                     double[][] dist, int[][] x, int n, int k) {
        Set<Integer> pointSet = new HashSet<>();
        int centerId = centers.get(j).id;
        // 直接计算索引：因为centerId就是原数据中的序号（从1开始），索引=id-1
        int centerIdx = centerId - 1;
        if (centerIdx < 0 || centerIdx >= allPoints.size()) {
            return new TSPPointSet(pointSet, i, x);
        }
        
        // 创建x矩阵的副本，用于本次计算
        int[][] xForThis = new int[n][k];
        for (int p = 0; p < n; p++) {
            System.arraycopy(x[p], 0, xForThis[p], 0, k);
        }
        
        // 对于遍历的点i，如果x_{i,j}=0，则将其设为1
        if (xForThis[i][j] == 0) {
            xForThis[i][j] = 1;
        }
        
        // 添加所有x_{i,j}=1的点到点集
        // 注意：中心点已经包含在内，因为x[centerIdx][j]=1（中心点到自己的距离为0，在computeXMatrix中会被设置为1）
        for (int idx = 0; idx < n; idx++) {
            if (xForThis[idx][j] == 1) {
                pointSet.add(idx);
            }
        }
        
        return new TSPPointSet(pointSet, i, xForThis);
    }

    // ----------------------------------------------------------------
    // 方法 A: LP 目标函数值计算
    // 模型形式：
    // min sum_{(k,i) in A} d_{k,i} y_{k,i}
    // s.t. sum_{k: (k,i) in A} y_{k,i} = sum_{k: (i,k) in A} y_{i,k},  forall i in C ∪ {d}  (流量平衡)
    //      sum_{k: (k,i) in A} y_{k,i} = 1,  forall i in C ∪ {d}  (访问约束)
    //      sum_{(k,i) in A(S)} y_{k,i} <= |S| - 1,  forall S ⊆ C, S ≠ ∅  (子回路消除)
    //      y_{k,i} in [0,1]  (LP松弛)
    // 对偶值通过计算"with i"和"without i"两种情况下的LP目标函数差值得到
    // ----------------------------------------------------------------
    
    /**
     * LP求解结果类，包含目标函数值和对偶值
     */
    static class LPResult {
        double objVal;
        double dualValue;  // 对应约束的对偶值（如果targetPointIdx >= 0）
        LPResult(double obj, double dual) { objVal = obj; dualValue = dual; }
    }
    
    /**
     * 求解TSP LP松弛的目标函数值
     * 模型形式：
     * min sum_{(k,i) in A} d_{k,i} y_{k,i}
     * s.t. sum_{k: (k,i) in A} y_{k,i} = sum_{k: (i,k) in A} y_{i,k},  forall i in C ∪ {d}  (流量平衡)
     *      sum_{k: (k,i) in A} y_{k,i} = 1,  forall i in C ∪ {d}  (访问约束)
     *      sum_{(k,i) in A(S)} y_{k,i} <= |S| - 1,  forall S ⊆ C, S ≠ ∅  (子回路消除)
     *      y_{k,i} in [0,1]  (LP松弛)
     * @param pointSet 点集（包含中心点和客户点）
     * @param distMatrix 距离矩阵
     * @param xMatrix x_{i,j}数组（此参数在新模型中不再使用，但保留以保持接口一致）
     * @param regionIdx 区域索引j
     * @param centerIdx 中心点索引（depot）
     * @param targetPointIdx 目标点i在全局索引中的位置（用于获取对偶值，如果<0则不获取对偶值）
     * @return LPResult包含目标函数值和对偶值
     */
    private static LPResult solveLPCost(Set<Integer> pointSet, double[][] distMatrix, int[][] xMatrix, int regionIdx, int centerIdx, int targetPointIdx) throws GRBException {
        List<Integer> nodes = new ArrayList<>(pointSet);
        int size = nodes.size();
        if(size < 2) return new LPResult(0.0, Double.NaN);
        
        // 找到中心点在nodes列表中的索引
        int depotIdx = -1;
        for(int i=0; i<size; i++) {
            if(nodes.get(i) == centerIdx) {
                depotIdx = i;
                break;
            }
        }
        if(depotIdx == -1) {
            // 如果中心点不在点集中，返回NaN
            return new LPResult(Double.NaN, Double.NaN);
        }
        
        // 计算客户点集合C的大小（不包括中心点）
        int numCustomers = size - 1;
        if(numCustomers == 0) return new LPResult(0.0, Double.NaN);
        
        // 找到目标点i在nodes列表中的索引（用于获取对偶值）
        int targetLocalIdx = -1;
        if(targetPointIdx >= 0) {
            for(int i=0; i<size; i++) {
                if(nodes.get(i) == targetPointIdx) {
                    targetLocalIdx = i;
                    break;
                }
            }
        }
        
        // 同步 Gurobi 操作以确保线程安全
        synchronized (gurobiLock) {
            GRBModel model = new GRBModel(env);
            try {
                // 变量：y_{k,i} (连续变量，0 <= y <= 1)
                GRBVar[][] y = new GRBVar[size][size];
                for(int u=0; u<size; u++) {
                    for(int v=0; v<size; v++) {
                        if(u != v) {
                            int idxU = nodes.get(u);
                            int idxV = nodes.get(v);
                            y[u][v] = model.addVar(0.0, 1.0, distMatrix[idxU][idxV], GRB.CONTINUOUS, null);
                        }
                    }
                }
                
                // 约束1：流量平衡（对所有点 i ∈ C ∪ {d}）
                // sum_{k} y_{k,i} = sum_{k} y_{i,k}, forall i ∈ C ∪ {d}
                for(int i=0; i<size; i++) {
                    GRBLinExpr inFlow = new GRBLinExpr();
                    GRBLinExpr outFlow = new GRBLinExpr();
                    for(int k=0; k<size; k++) {
                        if(k != i) {
                            inFlow.addTerm(1.0, y[k][i]);
                            outFlow.addTerm(1.0, y[i][k]);
                        }
                    }
                    // 所有点的流量平衡（包括depot）
                    model.addConstr(inFlow, GRB.EQUAL, outFlow, "flow_bal_"+i);
                }
                
                // 约束2：访问约束 sum_{k} y_{k,i} >= 1 (对所有点 i ∈ C ∪ {d})
                // 保存目标点对应的约束引用，用于获取对偶值
                GRBConstr targetVisitConstr = null;
                for(int i=0; i<size; i++) {
                    GRBLinExpr inFlow = new GRBLinExpr();
                    for(int k=0; k<size; k++) {
                        if(k != i) inFlow.addTerm(1.0, y[k][i]);
                    }
                    GRBConstr constr = model.addConstr(inFlow, GRB.GREATER_EQUAL, 1.0, "visit_"+i);
                    // 如果是目标点，保存约束引用
                    if(i == targetLocalIdx) {
                        targetVisitConstr = constr;
                    }
                }
                
                // 约束3：子回路消除（迭代添加）
                // sum_{(k,i) in A(S)} y_{k,i} <= |S| - 1, forall S ⊆ C, S ≠ ∅
                // 注意：子回路消除约束只针对客户点集合C的子集，不包括中心点（depot）
                // for(int iter=0; iter<50; iter++) {
                //     model.optimize();
                //     int status = model.get(GRB.IntAttr.Status);
                //     if(status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                //         // 如果求解失败，返回NaN
                //         return new LPResult(Double.NaN, Double.NaN);
                //     }
                    
                //     double[][] sol = new double[size][size];
                //     for(int u=0; u<size; u++) {
                //         for(int v=0; v<size; v++) {
                //             if(u != v) sol[u][v] = y[u][v].get(GRB.DoubleAttr.X);
                //         }
                //     }
                    
                //     // 查找子回路（只针对客户点集合C的子集，不包含depot）
                //     List<List<Integer>> subtours = findSubtours(sol, size);
                //     if(subtours.isEmpty()) break;
                    
                //     boolean added = false;
                //     for(List<Integer> sub : subtours) {
                //         // 子回路消除约束只针对客户点集合的子集 S ⊂ C（不包含depot）
                //         if(sub.size() >= 2 && sub.size() < size) {
                //             // 检查子回路是否包含中心点
                //             boolean containsDepot = sub.contains(depotIdx);
                //             if(!containsDepot) {
                //                 // 只对不包含depot的子回路（即只包含客户点的子回路）添加消除约束
                //                 GRBLinExpr expr = new GRBLinExpr();
                //                 for(int u : sub) {
                //                     for(int v : sub) {
                //                         if(u != v) expr.addTerm(1.0, y[u][v]);
                //                     }
                //                 }
                //                 model.addConstr(expr, GRB.LESS_EQUAL, sub.size()-1, "subtour_"+iter+"_"+sub.size());
                //                 added = true;
                //             }
                //         }
                //     }
                //     if(!added) break;
                // }
                
                // 最后再次优化，确保得到最优解
                model.optimize();
                int finalStatus = model.get(GRB.IntAttr.Status);
                if(finalStatus != GRB.OPTIMAL && finalStatus != GRB.SUBOPTIMAL) {
                    return new LPResult(Double.NaN, Double.NaN);
                }
                
                // 验证解的质量：检查是否还有明显的子回路
                double[][] finalSol = new double[size][size];
                for(int u=0; u<size; u++) {
                    for(int v=0; v<size; v++) {
                        if(u != v) finalSol[u][v] = y[u][v].get(GRB.DoubleAttr.X);
                    }
                }
                List<List<Integer>> finalSubtours = findSubtours(finalSol, size);
                // 如果还有明显的子回路（大小>=2且<size），说明求解可能不准确
                boolean hasSubtour = false;
                int subtourCount = 0;
                for(List<Integer> sub : finalSubtours) {
                    if(sub.size() >= 2 && sub.size() < size) {
                        hasSubtour = true;
                        subtourCount++;
                        // 输出警告信息：显示子回路的局部索引
                        // System.out.println(String.format(
                        //     "[警告] Region=%d, TargetPoint=%d: 检测到明显子回路，大小=%d，局部节点索引=%s。这可能是数值精度问题，求解结果可能不准确。",
                        //     regionIdx, targetPointIdx, sub.size(), sub.toString()));
                    }
                }
                if(hasSubtour) {
                    System.out.println(String.format(
                        "[警告] Region=%d, TargetPoint=%d: LP求解后仍存在 %d 个明显子回路（大小>=2且<总点数%d）。",
                        regionIdx, targetPointIdx, subtourCount, size));
                }
                
                // 获取LP目标函数值
                double objVal = model.get(GRB.DoubleAttr.ObjVal);
                // 确保目标函数值非负
                objVal = (objVal < 0) ? 0.0 : objVal;
                
                // 获取目标点对应约束的对偶值
                double dualVal = Double.NaN;
                if(targetVisitConstr != null && targetLocalIdx >= 0) {
                    try {
                        dualVal = targetVisitConstr.get(GRB.DoubleAttr.Pi);
                    } catch(GRBException e) {
                        dualVal = Double.NaN;
                    }
                }
                
                return new LPResult(objVal, dualVal);
            } finally { model.dispose(); }
        }
    }

    // ----------------------------------------------------------------
    // 方法 B: MIP 求解
    // 模型形式与LP相同，但 y_{k,i} 是二进制变量 {0,1}
    // min sum_{(k,i) in A} d_{k,i} y_{k,i}
    // s.t. sum_{k: (k,i) in A} y_{k,i} = sum_{k: (i,k) in A} y_{i,k},  forall i in C ∪ {d}  (流量平衡)
    //      sum_{k: (k,i) in A} y_{k,i} >= 1,  forall i in C ∪ {d}  (访问约束)
    //      sum_{(k,i) in A(S)} y_{k,i} <= |S| - 1,  forall S ⊆ C, S ≠ ∅  (子回路消除)
    //      y_{k,i} in {0,1}  (二进制变量)
    // ----------------------------------------------------------------
    
    // Callback for MIP Subtour Elimination
    static class MIPSubtourCallback extends GRBCallback {
        private GRBVar[][] y; int n; int depotIdx;
        MIPSubtourCallback(GRBVar[][] y, int n, int depotIdx) { 
            this.y = y; 
            this.n = n; 
            this.depotIdx = depotIdx;
        }
        @Override protected void callback() {
            try {
                if(where == GRB.CB_MIPSOL) {
                    double[][] sol = new double[n][n];
                    for(int i=0; i<n; i++) for(int j=0; j<n; j++) if(i!=j) sol[i][j] = getSolution(y[i][j]);
                    List<List<Integer>> subs = findSubtours(sol, n);
                    for(List<Integer> sub : subs) {
                        // 子回路消除约束只针对客户点集合的子集 S ⊆ C（不包含depot）
                        if(sub.size() >= 2 && sub.size() < n) {
                            // 检查子回路是否包含中心点
                            boolean containsDepot = sub.contains(depotIdx);
                            if(!containsDepot) {
                                // 只对不包含depot的子回路（即只包含客户点的子回路）添加消除约束
                                GRBLinExpr expr = new GRBLinExpr();
                                for(int u : sub) for(int v : sub) if(u!=v) expr.addTerm(1.0, y[u][v]);
                                addLazy(expr, GRB.LESS_EQUAL, sub.size()-1);
                            }
                        }
                    }
                }
            } catch(GRBException e) {}
        }
    }

    private static double solveMIPTsp(Set<Integer> pointSet, double[][] distMatrix, int[][] xMatrix, int regionIdx, int centerIdx) throws GRBException {
        List<Integer> nodes = new ArrayList<>(pointSet);
        int size = nodes.size();
        if(size < 2) return 0.0;
        
        // 找到中心点在nodes列表中的索引
        int depotIdx = -1;
        for(int i=0; i<size; i++) {
            if(nodes.get(i) == centerIdx) {
                depotIdx = i;
                break;
            }
        }
        if(depotIdx == -1) {
            // 如果中心点不在点集中，返回NaN
            return Double.NaN;
        }
        
        // 计算客户点集合C的大小（不包括中心点）
        int numCustomers = size - 1;
        if(numCustomers == 0) return 0.0;
        
        // 同步 Gurobi 操作以确保线程安全
        synchronized (gurobiLock) {
            GRBModel model = new GRBModel(env);
            try {
                // 变量：y_{k,i} (二进制变量)
                GRBVar[][] y = new GRBVar[size][size];
                for(int u=0; u<size; u++) {
                    for(int v=0; v<size; v++) {
                        if(u != v) {
                            int idxU = nodes.get(u);
                            int idxV = nodes.get(v);
                            y[u][v] = model.addVar(0.0, 1.0, distMatrix[idxU][idxV], GRB.BINARY, null);
                        }
                    }
                }
                
                // 约束1：流量平衡（对所有点 i ∈ C ∪ {d}）
                // sum_{k} y_{k,i} = sum_{k} y_{i,k}, forall i ∈ C ∪ {d}
                for(int i=0; i<size; i++) {
                    GRBLinExpr inFlow = new GRBLinExpr();
                    GRBLinExpr outFlow = new GRBLinExpr();
                    for(int k=0; k<size; k++) {
                        if(k != i) {
                            inFlow.addTerm(1.0, y[k][i]);
                            outFlow.addTerm(1.0, y[i][k]);
                        }
                    }
                    // 所有点的流量平衡（包括depot）
                    model.addConstr(inFlow, GRB.EQUAL, outFlow, "flow_bal_"+i);
                }
                
                // 约束2：访问约束 sum_{k} y_{k,i} >= 1 (对所有点 i ∈ C ∪ {d})
                for(int i=0; i<size; i++) {
                    GRBLinExpr inFlow = new GRBLinExpr();
                    for(int k=0; k<size; k++) {
                        if(k != i) inFlow.addTerm(1.0, y[k][i]);
                    }
                    model.addConstr(inFlow, GRB.GREATER_EQUAL, 1.0, "visit_"+i);
                }
                
                // 约束3：子回路消除（使用lazy constraint）
                // 注意：子回路消除约束只针对客户点集合C，不包括中心点（depot）
                // model.setCallback(new MIPSubtourCallback(y, size, depotIdx));
                // model.set(GRB.IntParam.LazyConstraints, 1);
                model.optimize();
                
                if(model.get(GRB.IntAttr.Status) == GRB.OPTIMAL) return model.get(GRB.DoubleAttr.ObjVal);
                return Double.NaN;
            } finally { model.dispose(); }
        }
    }

    // Shared Helper: Find Subtours
    private static List<List<Integer>> findSubtours(double[][] sol, int n) {
        List<List<Integer>> subs = new ArrayList<>();
        boolean[] vis = new boolean[n];
        for(int i=0; i<n; i++) {
            if(!vis[i]) {
                List<Integer> sub = new ArrayList<>();
                dfs(i, sol, vis, sub, n);
                if(sub.size()>1) subs.add(sub); // Note: for LP heuristic, we take all >1
            }
        }
        return subs;
    }
    private static void dfs(int u, double[][] sol, boolean[] vis, List<Integer> sub, int n) {
        vis[u] = true; sub.add(u);
        for(int v=0; v<n; v++) if(!vis[v] && sol[u][v] > 0.5) dfs(v, sol, vis, sub, n);
    }

    // ----------------------------------------------------------------
    // 主程序
    // ----------------------------------------------------------------
    
    public static void main(String[] args) {
        try {
            System.out.println("=== TSP对偶值验证器 ===");
            env = new GRBEnv();
            env.set(GRB.IntParam.OutputFlag, 0);
            env.set(GRB.IntParam.Threads, 1);
            env.start();
            
            // 小数据集路径
            String pFile = "data/unique_coordinates_list_100.csv";
            String cFile = "data/selected_centers_100points_p5.csv";
            String outFile = "output/validation_results_4.csv";
            
            System.out.println("Reading data...");
            List<Point> allPoints = readCoordinates(pFile);
            List<Point> centers = readCoordinates(cFile); 
            double[][] distMatrix = computeDistanceMatrix(allPoints);
            
            int n = allPoints.size();
            int k = centers.size();
            System.out.println("Validation Start: N=" + n + ", K=" + k + ", Tasks=" + (n*k));
            
            // 预先计算x_{i,j}矩阵
            System.out.println("Computing x_{i,j} matrix...");
            int[][] xMatrix = computeXMatrix(allPoints, centers, distMatrix, n, k);
            
            BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
            writer.write("Region,Point,Dual_Value,MIP_With_i,MIP_Without_i,MIP_Diff,Gap,LP_With_i,LP_Without_i,Dual_From_Model\n");
            
            AtomicInteger processed = new AtomicInteger(0);
            long startTime = System.currentTimeMillis();
            
            // 并行验证
            IntStream.range(0, k).parallel().forEach(j -> {
                for(int i=0; i<n; i++) {
                    try {
                        // 1. 构建"without i"的点集和x矩阵（使用原始的xMatrix）
                        int[][] xWithoutI = new int[n][k];
                        for (int p = 0; p < n; p++) {
                            System.arraycopy(xMatrix[p], 0, xWithoutI[p], 0, k);
                        }
                        // 确保x_{i,j}=0（不包含点i）
                        xWithoutI[i][j] = 0;
                        
                        // 重新构建点集，只包含xWithoutI[idx][j]=1的点
                        Set<Integer> setWithoutI = new HashSet<>();
                        int centerId = centers.get(j).id;
                        // 直接计算索引：因为centerId就是原数据中的序号（从1开始），索引=id-1
                        int centerIdx = centerId - 1;
                        if (centerIdx >= 0 && centerIdx < allPoints.size()) {
                            setWithoutI.add(centerIdx);
                        }
                        for (int idx = 0; idx < n; idx++) {
                            if (xWithoutI[idx][j] == 1) {
                                setWithoutI.add(idx);
                            }
                        }
                        
                        // 2. 构建"with i"的点集和x矩阵
                        TSPPointSet setWithI = determineTSPPointSet(j, i, allPoints, centers, distMatrix, xMatrix, n, k);
                        
                        // 3. 计算LP Dual (LP目标函数差值)
                        // 方法：比较"without i"和"with i"两种情况下的LP目标函数差值
                        // 这个差值反映"加入点i"的边际成本（LP松弛下的）
                        // 注意：如果点集大小小于2，LP成本为0
                        // 注意：中心点作为depot，需要传递给求解方法
                        LPResult lpResultWithoutI = (setWithoutI.size() < 2) 
                            ? new LPResult(0.0, Double.NaN) 
                            : solveLPCost(setWithoutI, distMatrix, xWithoutI, j, centerIdx, -1);
                        LPResult lpResultWithI = (setWithI.pointSet.size() < 2) 
                            ? new LPResult(0.0, Double.NaN) 
                            : solveLPCost(setWithI.pointSet, distMatrix, setWithI.xMatrix, j, centerIdx, i);
                        
                        double lpCostWithoutI = lpResultWithoutI.objVal;
                        double lpCostWithI = lpResultWithI.objVal;
                        double dual = (Double.isNaN(lpCostWithoutI) || Double.isNaN(lpCostWithI)) 
                            ? Double.NaN : (lpCostWithI - lpCostWithoutI);
                        double dualFromModel = lpResultWithI.dualValue;
                        
                        // 计算MIP差值
                        double costWithI = solveMIPTsp(setWithI.pointSet, distMatrix, setWithI.xMatrix, j, centerIdx);
                        double costWithoutI = (setWithoutI.size() < 2) ? 0.0 : solveMIPTsp(setWithoutI, distMatrix, xWithoutI, j, centerIdx);
                        double mipDiff = (Double.isNaN(costWithI) || Double.isNaN(costWithoutI)) 
                            ? Double.NaN : (costWithI - costWithoutI);
                        
                        // 4. Gap = (MIP_Diff - Dual_Value) / MIP_Diff
                        double gap = (Double.isNaN(dual) || Double.isNaN(mipDiff) || mipDiff == 0.0) 
                            ? Double.NaN : ((mipDiff - dual) / mipDiff);
                        
                        synchronized(writer) {
                            writer.write(String.format("%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", 
                                j, i, dual, costWithI, costWithoutI, mipDiff, gap, 
                                lpCostWithI, lpCostWithoutI, dualFromModel));
                        }
                        
                        int c = processed.incrementAndGet();
                        if(c % 50 == 0) {
                            long elapsed = System.currentTimeMillis() - startTime;
                             System.out.println("Processed " + c + "/" + (n*k) + ", Time: " + (elapsed/1000.0) + "s");
                        }
                        
                    } catch(Exception e) {
                        e.printStackTrace();
                    }
                }
            });
            
            writer.close();
            env.dispose();
            System.out.println("Validation Done. Results in " + outFile);
            
        } catch(Exception e) { e.printStackTrace(); }
    }
}

