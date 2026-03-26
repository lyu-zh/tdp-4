import com.gurobi.gurobi.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * TSP对偶值提取器 (MIP差值版本 - 优化版)
 * 1. 使用并行流加速 N*K 次求解
 * 2. 使用两次MIP求解计算目标函数值差值（包含点i vs 不包含点i）
 * 3. 共享 Gurobi Environment 提升性能
 * 4. 【优化】对于每个区域 j，基础点集的 TSP（不包含点 i）只计算一次，避免重复计算
 *    对于固定的区域 j，不包含点 i 的 TSP 对所有点 i 都是一样的（中心点 + 最近的 n/k 个点）
 */
public class TSPDualExtractor {
    
    // 地球半径（米）
    private static final double EARTH_RADIUS_METERS = 6371000.0;
    
    // 特殊值：用于标记不在5n/k范围内的点（在CSV中显示为"Null"）
    private static final double OUT_OF_RANGE_MARKER = Double.NEGATIVE_INFINITY;
    
    // 全局 Gurobi 环境（线程安全，避免重复创建）
    private static GRBEnv env;

    /**
     * 点类
     */
    public static class Point {
        int id;
        double longitude;
        double latitude;
        
        public Point(int id, double longitude, double latitude) {
            this.id = id;
            this.longitude = longitude;
            this.latitude = latitude;
        }
    }
    
    /**
     * 使用 Haversine 公式计算两点间距离（米）
     */
    private static double haversineDistance(double lon1, double lat1, double lon2, double lat2) {
        double lat1Rad = Math.toRadians(lat1);
        double lat2Rad = Math.toRadians(lat2);
        double deltaLatRad = Math.toRadians(lat2 - lat1);
        double deltaLonRad = Math.toRadians(lon2 - lon1);
        
        double a = Math.sin(deltaLatRad / 2) * Math.sin(deltaLatRad / 2) +
                   Math.cos(lat1Rad) * Math.cos(lat2Rad) *
                   Math.sin(deltaLonRad / 2) * Math.sin(deltaLonRad / 2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        
        return EARTH_RADIUS_METERS * c;
    }
    
    // ---------------- 读取数据部分 ----------------
    
    public static List<Point> readCoordinates(String filename) throws IOException {
        List<Point> points = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = reader.readLine(); // skip header
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            String[] parts = line.split(",");
            if (parts.length >= 3) {
                try {
                    // 第一列支持整数或浮点（如 70 或 70.0），统一转为 int
                    int id = (int) Double.parseDouble(parts[0].trim());
                    points.add(new Point(id,
                                       Double.parseDouble(parts[1].trim()),
                                       Double.parseDouble(parts[2].trim())));
                } catch (Exception e) {}
            }
        }
        reader.close();
        return points;
    }

    public static List<Point> readCenters(String filename) throws IOException {
        return readCoordinates(filename);
    }

    public static int findPointIndexById(List<Point> points, int id) {
        for (int i = 0; i < points.size(); i++) {
            if (points.get(i).id == id) return i;
        }
        return -1;
    }
    
    /**
     * 从 demand_matrix.csv 中读取指定日期的列作为 xi 向量
     * 文件格式：第一行是表头（经度,纬度,日期1,日期2,...），后续行是数据
     * 根据经度和纬度匹配点，提取指定日期列的值作为 xi 向量
     * @param filename demand_matrix.csv 文件路径
     * @param targetDate 目标日期，格式：yyyy-MM-dd（例如："2023-09-03"）
     * @param allPoints 所有点的列表（用于匹配经度和纬度）
     * @return xi 向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     *         如果文件不存在或读取失败，返回 null（表示不使用 xi 过滤）
     */
    public static int[] readXiVectorFromDemandMatrix(String filename, String targetDate, List<Point> allPoints) {
        try {
            int n = allPoints.size();
            int[] xi = new int[n];
            // 初始化所有点为 0
            Arrays.fill(xi, 0);
            
            BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename), "UTF-8"));
            
            // 读取表头，找到目标日期对应的列索引
            String headerLine = reader.readLine();
            if (headerLine == null) {
                reader.close();
                System.out.println("警告: demand_matrix.csv 文件为空，将不使用 xi 过滤");
                return null;
            }
            
            // 解析表头，找到目标日期列的索引
            String[] headers = headerLine.split(",");
            int targetDateColIndex = -1;
            for (int i = 0; i < headers.length; i++) {
                if (headers[i].trim().equals(targetDate)) {
                    targetDateColIndex = i;
                    break;
                }
            }
            
            if (targetDateColIndex == -1) {
                reader.close();
                System.out.println("警告: 在 demand_matrix.csv 中未找到日期 " + targetDate + "，将不使用 xi 过滤");
                return null;
            }
            
            System.out.println("找到目标日期列: " + targetDate + " (列索引: " + targetDateColIndex + ")");
            
            // 读取数据行，根据经度和纬度匹配点
            String line;
            int matchedCount = 0;
            double tolerance = 1e-8; // 经度纬度匹配的容差
            
            // 【修复】检查表头格式：如果第一列是"序号"或"序号"相关的列，则经度纬度列需要偏移
            boolean hasIndexColumn = false;
            if (headers.length > 0) {
                String firstHeader = headers[0].trim().toLowerCase();
                // 检查第一列是否是序号列（可能是"序号"、"id"、"index"等）
                if (firstHeader.contains("序号") || firstHeader.contains("id") || 
                    firstHeader.contains("index") || firstHeader.equals("")) {
                    hasIndexColumn = true;
                }
            }
            
            // 【修复】如果第一列是序号列，则经度在parts[1]，纬度在parts[2]
            // targetDateColIndex已经是相对于完整表头的索引，不需要调整
            int lonColIndex = hasIndexColumn ? 1 : 0;
            int latColIndex = hasIndexColumn ? 2 : 1;
            int dateColIndex = targetDateColIndex; // targetDateColIndex已经是正确的索引
            
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                
                String[] parts = line.split(",");
                if (parts.length <= dateColIndex) continue;
                
                try {
                    double lon = Double.parseDouble(parts[lonColIndex].trim());
                    double lat = Double.parseDouble(parts[latColIndex].trim());
                    int value = Integer.parseInt(parts[dateColIndex].trim());
                    
                    // 将值转换为 0 或 1
                    if (value != 0 && value != 1) {
                        value = (value > 0) ? 1 : 0;
                    }
                    
                    // 在 allPoints 中查找匹配的点（根据经度和纬度）
                    for (int i = 0; i < n; i++) {
                        Point p = allPoints.get(i);
                        if (Math.abs(p.longitude - lon) < tolerance && 
                            Math.abs(p.latitude - lat) < tolerance) {
                            xi[i] = value;
                            matchedCount++;
                            break;
                        }
                    }
                } catch (NumberFormatException e) {
                    // 忽略无效行
                }
            }
            
            reader.close();
            
            int accessibleCount = Arrays.stream(xi).sum();
            System.out.println("成功从 demand_matrix.csv 读取 xi 向量:");
            System.out.println("  - 目标日期: " + targetDate);
            System.out.println("  - 匹配的点数: " + matchedCount + "/" + n);
            System.out.println("  - 可访问的点数: " + accessibleCount + "/" + n);
            
            return xi;
        } catch (IOException e) {
            System.out.println("读取 demand_matrix.csv 失败: " + e.getMessage() + "，将不使用 xi 过滤");
            return null;
        }
    }
    
    /**
     * 读取 xi 向量（0-1向量）- 旧版本，保留用于兼容性
     * 文件格式：每行一个整数，0 或 1
     * 如果文件不存在或读取失败，返回 null（表示不使用 xi 过滤）
     * @param filename 文件名
     * @param n 点的数量
     * @return xi 向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     */
    private static int[] readXiVector(String filename, int n) {
        try {
            int[] xi = new int[n];
            BufferedReader reader = new BufferedReader(new FileReader(filename));
            String line;
            int idx = 0;
            while ((line = reader.readLine()) != null && idx < n) {
                line = line.trim();
                if (line.isEmpty()) continue;
                try {
                    xi[idx] = Integer.parseInt(line);
                    if (xi[idx] != 0 && xi[idx] != 1) {
                        System.out.println("警告: xi[" + idx + "] = " + xi[idx] + "，应该是 0 或 1，已设置为 0");
                        xi[idx] = 0;
                    }
                    idx++;
                } catch (NumberFormatException e) {
                    // 忽略无效行
                }
            }
            reader.close();
            
            // 如果读取的点数不足 n，剩余的点设置为 0
            for (int i = idx; i < n; i++) {
                xi[i] = 0;
            }
            
            System.out.println("成功读取 xi 向量，共 " + idx + " 个点，其中 " + 
                Arrays.stream(xi).sum() + " 个点可以被访问");
            return xi;
        } catch (IOException e) {
            System.out.println("未找到 xi 向量文件: " + filename + "，将不使用 xi 过滤");
            return null;
        }
    }
    
    // ---------------- 核心计算部分 ----------------

    /**
     * 计算距离矩阵（并行优化，带进度提醒）
     */
    // 全局日志记录器（用于距离矩阵计算等非并行计算部分）
    private static Logger globalLogger = null;
    
    public static void setGlobalLogger(Logger logger) {
        globalLogger = logger;
    }
    
    public static double[][] computeDistanceMatrix(List<Point> points) {
        int n = points.size();
        double[][] dist = new double[n][n];
        
        final Logger loggerRef = globalLogger; // 创建 final 引用
        String msgStart = "正在使用 " + Runtime.getRuntime().availableProcessors() + " 个线程并行计算距离矩阵...\n";
        if (loggerRef != null) loggerRef.log(msgStart);
        else System.out.print(msgStart);
        AtomicInteger completedRows = new AtomicInteger(0);
        long startTime = System.currentTimeMillis();
        
        IntStream.range(0, n).parallel().forEach(i -> {
            Point p1 = points.get(i);
            for (int j = i + 1; j < n; j++) {
                Point p2 = points.get(j);
                double distance = haversineDistance(p1.longitude, p1.latitude,
                                                   p2.longitude, p2.latitude);
                dist[i][j] = distance;
                dist[j][i] = distance;
            }
            dist[i][i] = 0.0;
            
            // 进度提醒
            int c = completedRows.incrementAndGet();
            if (n >= 10 && c % (n / 10) == 0) { // 每完成10%提示一次
                String msgProgress = String.format("距离矩阵计算进度: %d%% (%d/%d), 耗时: %.1fs%n", 
                    (c * 100 / n), c, n, (System.currentTimeMillis() - startTime) / 1000.0);
                if (loggerRef != null) loggerRef.log(msgProgress);
                else System.out.print(msgProgress);
            }
        });
        
        String msgEnd = "距离矩阵计算完成！\n";
        if (loggerRef != null) loggerRef.log(msgEnd);
        else System.out.print(msgEnd);
        return dist;
    }
    
    static class DistancePair {
        int index; double distance;
        DistancePair(int index, double distance) { this.index = index; this.distance = distance; }
    }
    
    /**
     * 获取离中心点所有点的排序后的距离列表（只排序一次，可复用）
     * 【优化】避免重复排序：对同一个中心点的所有点只排序一次，然后可以多次使用
     * @param centerIdx 中心点索引
     * @param distMatrix 距离矩阵
     * @param n 总点数
     * @return 按距离排序后的距离对列表（不包括中心点本身）
     */
    private static List<DistancePair> getSortedDistances(int centerIdx, double[][] distMatrix, int n) {
        List<DistancePair> distances = new ArrayList<>(n);
        
        // 收集所有点（除了中心点）的距离信息
        for (int idx = 0; idx < n; idx++) {
            if (idx != centerIdx) {
                double d = distMatrix[centerIdx][idx];
                distances.add(new DistancePair(idx, d));
            }
        }
        
        // 按距离排序（只排序一次）
        distances.sort(Comparator.comparingDouble(p -> p.distance));
        
        return distances;
    }
    
    /**
     * 从排序后的距离列表中获取最近的 numPoints 个点的索引列表
     * 【优化】使用预排序的距离列表，避免重复排序
     * @param sortedDistances 已排序的距离列表
     * @param numPoints 要获取的点数（例如：5*n/k）
     * @return 最近的 numPoints 个点的索引列表
     */
    private static List<Integer> getNearestPointsFromSorted(List<DistancePair> sortedDistances, int numPoints) {
        List<Integer> nearestPoints = new ArrayList<>();
        for (int idx = 0; idx < numPoints && idx < sortedDistances.size(); idx++) {
            nearestPoints.add(sortedDistances.get(idx).index);
        }
        return nearestPoints;
    }
    
    /**
     * 获取离区域 j 中心点最近的 numPoints 个点的索引列表（保留用于兼容性）
     * 【已废弃】建议使用 getSortedDistances + getNearestPointsFromSorted 以避免重复排序
     */
    private static List<Integer> getNearestPoints(int centerIdx, double[][] distMatrix, int n, int numPoints) {
        List<DistancePair> sortedDistances = getSortedDistances(centerIdx, distMatrix, n);
        return getNearestPointsFromSorted(sortedDistances, numPoints);
    }
    
    /**
     * 从排序后的距离列表中获取区域 j 的基础点集：中心点 + 最近的 n/k 个点（不包含任何额外的点 i）
     * 【优化】使用预排序的距离列表，避免重复排序
     * 【修复】根据 xi 向量过滤点集：先从所有点中选择最近的 n/k 个点，然后只保留其中 xi_i = 1 的点
     * @param centerIdx 中心点索引
     * @param sortedDistances 已排序的距离列表（从 getSortedDistances 获取）
     * @param n 总点数
     * @param k 区域数
     * @param xi 0-1向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     *           如果为 null，则不进行过滤（保持原有行为）
     * @return 基础点集
     */
    private static Set<Integer> getBasePointSetFromSorted(int centerIdx,
                                                           List<DistancePair> sortedDistances,
                                                           int n, int k,
                                                           int[] xi) {
        Set<Integer> pointSet = new HashSet<>();
        
        // 【修改】无论中心点是否可以被访问，都要将其加入点集（中心点充当模型中的depot d）
        // 中心点总是加入点集
        pointSet.add(centerIdx);
        
        int numBasePoints = n / k;
        // 从排序后的距离列表中，取前 n/k 个点，只保留其中 xi_i = 1 的点（如果提供了 xi 向量）
        for (int idx = 0; idx < numBasePoints && idx < sortedDistances.size(); idx++) {
            int pointIdx = sortedDistances.get(idx).index;
            // 如果提供了 xi 向量，只保留 xi[pointIdx] = 1 的点
            if (xi == null || xi[pointIdx] == 1) {
                pointSet.add(pointIdx);
            }
        }
        
        // 注意：如果点集中只有中心点（大小为1），solveTSPMIP 会自动返回 0.0（不需要做TSP）
        
        return pointSet;
    }
    
    /**
     * 获取区域 j 的基础点集：中心点 + 最近的 n/k 个点（不包含任何额外的点 i）
     * 【优化】直接使用预计算的 distMatrix，避免重复计算
     * 【修复】根据 xi 向量过滤点集：先从所有点中选择最近的 n/k 个点，然后只保留其中 xi_i = 1 的点
     * 【已废弃】建议使用 getSortedDistances + getBasePointSetFromSorted 以避免重复排序
     * @param xi 0-1向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     *           如果为 null，则不进行过滤（保持原有行为）
     */
    private static Set<Integer> getBasePointSetForRegion(int j,
                                                          List<Point> allPoints,
                                                          List<Point> centers,
                                                          double[][] distMatrix,
                                                          int n, int k,
                                                          int[] xi) {
        int centerId = centers.get(j).id;
        int centerIdx = findPointIndexById(allPoints, centerId);
        
        if (centerIdx == -1) return new HashSet<>();
        
        List<DistancePair> sortedDistances = getSortedDistances(centerIdx, distMatrix, n);
        return getBasePointSetFromSorted(centerIdx, sortedDistances, n, k, xi);
    }
    
    /**
     * 重载方法：不提供 xi 向量时，保持原有行为
     */
    private static Set<Integer> getBasePointSetForRegion(int j,
                                                          List<Point> allPoints,
                                                          List<Point> centers,
                                                          double[][] distMatrix,
                                                          int n, int k) {
        return getBasePointSetForRegion(j, allPoints, centers, distMatrix, n, k, null);
    }
    
    /**
     * 确定包含点 i 的 TSP 点集：基础点集 + 额外点 i
     * 【优化】直接使用预计算的基础点集，避免重复计算
     * 【新增】根据 xi 向量过滤点集，只保留 xi_i = 1 的点
     * 【修改】中心点（depot）总是被包含，无论其 xi 值如何
     * @param basePointSet 基础点集
     * @param i 要添加的点索引
     * @param centerIdx 中心点（depot）的索引，总是被包含在点集中
     * @param xi 0-1向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     *           如果为 null，则不进行过滤（保持原有行为）
     */
    private static List<Integer> getPointSetWithI(Set<Integer> basePointSet, int i, int centerIdx, int[] xi) {
        Set<Integer> pointSet = new HashSet<>();
        
        // 【修改】中心点（depot）总是被包含，无论其 xi 值如何
        if (basePointSet.contains(centerIdx)) {
            pointSet.add(centerIdx);
        }
        
        // 添加基础点集中 xi = 1 的点（不包括中心点，因为已经添加）
        for (int idx : basePointSet) {
            if (idx == centerIdx) {
                continue; // 中心点已经添加
            }
            if (xi == null || xi[idx] == 1) {
                pointSet.add(idx);
            }
        }
        
        // 添加点 i（如果 xi[i] = 1，且 i 不是中心点）
        if (i != centerIdx && (xi == null || xi[i] == 1)) {
            pointSet.add(i);
        }
        
        return new ArrayList<>(pointSet);
    }
    
    /**
     * 重载方法：不提供 xi 向量时，保持原有行为
     */
    private static List<Integer> getPointSetWithI(Set<Integer> basePointSet, int i, int centerIdx) {
        return getPointSetWithI(basePointSet, i, centerIdx, null);
    }
    
    /**
     * 构建并求解TSP MIP模型
     * @param pointList 点集索引列表
     * @param distMatrix 距离矩阵
     * @return 目标函数值，如果求解失败返回Double.NaN
     */
    private static double solveTSPMIP(List<Integer> pointList, double[][] distMatrix) throws GRBException {
        int tspSize = pointList.size();
        if (tspSize < 2) return 0.0;
        
        GRBModel model = new GRBModel(env);
        
        try {
            // 【优化】设置超时时间，避免某些难解问题占用过多时间
            // 根据问题规模设置超时：小问题(<20点)5秒，中等问题(20-50点)10秒，大问题(>50点)30秒
            double timeLimit = tspSize < 20 ? 5.0 : (tspSize < 50 ? 10.0 : 30.0);
            model.set(GRB.DoubleParam.TimeLimit, timeLimit);
            
            // 【优化】设置 MIP 求解参数，加速求解
            model.set(GRB.IntParam.MIPFocus, 1); // 专注于找到可行解
            model.set(GRB.DoubleParam.MIPGap, 0.01); // 允许1%的gap，加速求解
            GRBVar[][] y = new GRBVar[tspSize][tspSize];
            
            // 1. 定义变量 (MIP: 0 <= y <= 1, BINARY)
            for (int u = 0; u < tspSize; u++) {
                for (int v = 0; v < tspSize; v++) {
                    if (u != v) {
                        int idxU = pointList.get(u);
                        int idxV = pointList.get(v);
                        y[u][v] = model.addVar(0.0, 1.0, distMatrix[idxU][idxV], 
                                             GRB.BINARY, "y_" + u + "_" + v);
                    }
                }
            }
            
            // 2. Flow Balance Constraints (流量平衡: 入度 = 出度)
            for (int idx = 0; idx < tspSize; idx++) {
                GRBLinExpr balanceExpr = new GRBLinExpr();
                // Add in-degree terms (+1)
                for (int u = 0; u < tspSize; u++) {
                    if (u != idx) balanceExpr.addTerm(1.0, y[u][idx]);
                }
                // Add out-degree terms (-1)
                for (int v = 0; v < tspSize; v++) {
                    if (v != idx) balanceExpr.addTerm(-1.0, y[idx][v]);
                }
                model.addConstr(balanceExpr, GRB.EQUAL, 0.0, "balance_" + idx);
            }
            
            // 3. Visit Constraints (访问约束: 入度 >= 1)
            for (int idx = 0; idx < tspSize; idx++) {
                GRBLinExpr inExpr = new GRBLinExpr();
                for (int u = 0; u < tspSize; u++) {
                    if (u != idx) inExpr.addTerm(1.0, y[u][idx]);
                }
                model.addConstr(inExpr, GRB.GREATER_EQUAL, 1.0, "visit_" + idx);
            }
            
            // 4. 使用MTZ约束消除子回路 (Miller-Tucker-Zemlin formulation)
            // 对于MIP，MTZ约束比迭代添加约束更高效
            if (tspSize > 2) {
                GRBVar[] u = new GRBVar[tspSize];
                for (int i = 1; i < tspSize; i++) {
                    u[i] = model.addVar(1.0, tspSize - 1.0, 0.0, GRB.CONTINUOUS, "u_" + i);
                }
                
                for (int i = 1; i < tspSize; i++) {
                    for (int j = 1; j < tspSize; j++) {
                        if (i != j) {
                            GRBLinExpr mtzExpr = new GRBLinExpr();
                            mtzExpr.addTerm(1.0, u[i]);
                            mtzExpr.addTerm(-1.0, u[j]);
                            mtzExpr.addTerm((double)(tspSize - 1), y[i][j]);
                            model.addConstr(mtzExpr, GRB.LESS_EQUAL, (double)(tspSize - 2), "mtz_" + i + "_" + j);
                        }
                    }
                }
            }
            
            // 5. 求解MIP
            model.optimize();
            
            int status = model.get(GRB.IntAttr.Status);
            // 【优化】接受最优解、可行解或超时后的最佳解
            if (status == GRB.OPTIMAL || status == GRB.TIME_LIMIT) {
                // 如果超时，尝试获取当前最佳解
                if (status == GRB.TIME_LIMIT) {
                    try {
                        double objVal = model.get(GRB.DoubleAttr.ObjVal);
                        if (!Double.isNaN(objVal) && !Double.isInfinite(objVal)) {
                            return objVal;
                        }
                    } catch (GRBException e) {
                        // 如果无法获取目标值，返回 NaN
                    }
                } else {
                    return model.get(GRB.DoubleAttr.ObjVal);
                }
            }
            return Double.NaN;
            
        } catch (GRBException e) {
            return Double.NaN;
        } finally {
            model.dispose();
        }
    }
    
    /**
     * 计算包含点 i 的 TSP 与不包含点 i 的 TSP 的差值（MIP_DIFF）
     * 参考 TSPDualValidator.java 的正确实现
     * 【修改】中心点（depot）总是被包含在点集中，无论其 xi 值如何
     * @param pointListWithI 包含点 i 的点集（已经根据 xi 过滤，但中心点总是被包含）
     * @param basePointSet 基础点集（中心点 + 最近的 n/k 个点，已经根据 xi 过滤）
     * @param pointI 点 i 的索引
     * @param centerIdx 中心点（depot）的索引，总是被包含在点集中
     * @param baseObjVal 基础点集的 TSP 目标值（缓存值，如果点 i 不在基础点集中则使用）
     * @param distMatrix 距离矩阵
     * @param xi 0-1向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     *           如果为 null，则不进行过滤（保持原有行为）
     * @return 差值 = 包含点 i 的 TSP 目标值 - 不包含点 i 的 TSP 目标值
     */
    private static double computeDualValue(List<Integer> pointListWithI,
                                          Set<Integer> basePointSet,
                                          int pointI,
                                          int centerIdx,
                                          double baseObjVal,
                                          double[][] distMatrix,
                                          int[] xi) throws GRBException {
        if (pointListWithI.size() < 2) {
            return 0.0;
        }
        
        // 【修改】如果点 i 是中心点（depot），那么 "with i" 和 "without i" 都包含中心点
        // 根据模型，中心点作为 depot 总是存在，所以差值应该是0
        if (pointI == centerIdx) {
            return 0.0;
        }
        
        // 计算包含点 i 的 TSP（点集已经根据 xi 过滤，但中心点总是被包含）
        double objValWithI = solveTSPMIP(pointListWithI, distMatrix);
        if (Double.isNaN(objValWithI)) {
            return Double.NaN;
        }
        
        // 判断点 i 是否在基础点集中
        boolean iInBaseSet = basePointSet.contains(pointI);
        
        double objValWithoutI;
        if (iInBaseSet) {
            // 如果点 i 在基础点集中，需要计算"不包含点 i"的点集
            // 【修改】中心点（depot）总是被保留，即使点 i 在基础点集中
            Set<Integer> pointSetWithoutI = new HashSet<>(basePointSet);
            pointSetWithoutI.remove(pointI);
            // 确保中心点总是被包含
            if (!pointSetWithoutI.contains(centerIdx) && basePointSet.contains(centerIdx)) {
                pointSetWithoutI.add(centerIdx);
            }
            List<Integer> pointListWithoutI = new ArrayList<>(pointSetWithoutI);
            
            if (pointListWithoutI.size() < 2) {
                objValWithoutI = 0.0;
            } else {
                objValWithoutI = solveTSPMIP(pointListWithoutI, distMatrix);
                if (Double.isNaN(objValWithoutI)) {
                    return Double.NaN;
                }
            }
        } else {
            // 如果点 i 不在基础点集中，使用缓存的基础点集 TSP 值
            if (Double.isNaN(baseObjVal)) {
                return objValWithI;
            }
            objValWithoutI = baseObjVal;
        }
        
        // 返回差值：包含点 i 的成本 - 不包含点 i 的成本
        return objValWithI - objValWithoutI;
    }
    
    /**
     * 重载方法：不提供 xi 向量时，保持原有行为
     */
    private static double computeDualValue(List<Integer> pointListWithI,
                                          Set<Integer> basePointSet,
                                          int pointI,
                                          int centerIdx,
                                          double baseObjVal,
                                          double[][] distMatrix) throws GRBException {
        return computeDualValue(pointListWithI, basePointSet, pointI, centerIdx, baseObjVal, distMatrix, null);
    }

    // ---------------- 主程序入口 ----------------
            // 阈值 > 0.5 确保在LP解中找到主要的“闭环”

    // ---------------- 主程序入口 ----------------
    
    /**
     * 线程安全的日志记录器：同时输出到控制台和文件
     */
    private static class Logger {
        private final PrintWriter fileWriter;
        private final Object lock = new Object();
        
        public Logger(String logFile) throws IOException {
            File file = new File(logFile);
            file.getParentFile().mkdirs(); // 确保目录存在
            this.fileWriter = new PrintWriter(new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(file, true), "UTF-8")));
        }
        
        public void log(String message) {
            synchronized (lock) {
                System.out.print(message);
                fileWriter.print(message);
                fileWriter.flush(); // 立即刷新，确保实时写入
            }
        }
        
        public void logf(String format, Object... args) {
            String message = String.format(format, args);
            log(message);
        }
        
        public void close() {
            synchronized (lock) {
                fileWriter.close();
            }
        }
    }
    
    /**
     * 核心计算逻辑：计算给定 xi 向量的 travelDist 矩阵
     * @param allPoints 所有点的列表
     * @param centers 中心点列表
     * @param distMatrix 距离矩阵
     * @param xi xi 向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     * @param env Gurobi 环境（如果为 null，将创建新的环境）
     * @param logFile 日志文件路径（如果为 null，则不记录到文件）
     * @return travelDist 矩阵，travelDist[i][j] 表示点 i 在区域 j 的对偶值
     */
    public static double[][] computeTravelDist(List<Point> allPoints, List<Point> centers, 
                                               double[][] distMatrix, int[] xi, GRBEnv env, String logFile) {
        int n = allPoints.size();
        int k = centers.size();
        if (k == 0) {
            throw new IllegalArgumentException("区域数 K=0：中心点文件为空或格式有误（第一列须为整数或浮点数，如 70 或 70.0）");
        }
        
        // 如果提供了环境，使用它；否则创建新的环境
        GRBEnv originalEnv = TSPDualExtractor.env;
        boolean needDisposeEnv = false;
        if (env != null) {
            TSPDualExtractor.env = env;
        } else if (TSPDualExtractor.env == null) {
            try {
                TSPDualExtractor.env = new GRBEnv();
                TSPDualExtractor.env.set(GRB.IntParam.OutputFlag, 0);
                TSPDualExtractor.env.set(GRB.IntParam.Threads, 1);
                TSPDualExtractor.env.start();
                needDisposeEnv = true;
            } catch (GRBException e) {
                throw new RuntimeException("无法创建 Gurobi 环境", e);
            }
        }
        
        // 创建日志记录器
        Logger logger = null;
        try {
            if (logFile != null) {
                logger = new Logger(logFile);
            }
        } catch (IOException e) {
            System.err.println("无法创建日志文件: " + logFile + "，将只输出到控制台");
        }
        
        // 创建 final 引用，以便在 lambda 中使用
        final Logger loggerFinal = logger;
        
        try {
            // 结果存储
            double[][] travelDist = new double[n][k];
            
            // 并行求解
            String msgStart = "开始并行求解...\n";
            if (loggerFinal != null) loggerFinal.log(msgStart);
            else System.out.print(msgStart);
            long startTime = System.currentTimeMillis();
            AtomicInteger completed = new AtomicInteger(0);
            // 【优化】总任务数：每个区域j只处理最近的5*n/k个点
            int numPointsPerRegion = 5 * n / k;
            int total = numPointsPerRegion * k;
            
            // 使用 IntStream 并行处理
            final int[] xiFinal = xi; // 为了在 lambda 表达式中使用
            final List<Point> allPointsFinal = allPoints; // 为了在 lambda 表达式中使用
            
            // 【优化】记录每个区域的计算时间，用于分析负载不均衡问题
            Map<Integer, Long> regionTimings = Collections.synchronizedMap(new HashMap<>());
            
            IntStream.range(0, k).parallel().forEach(j -> {
                long regionStartTime = System.currentTimeMillis();
                try {
                    // 步骤1：获取区域 j 的中心点索引
                    int centerId = centers.get(j).id;
                    int centerIdx = findPointIndexById(allPointsFinal, centerId);
                    
                    // 【优化】步骤2：只排序一次，获取排序后的距离列表（可复用）
                    List<DistancePair> sortedDistances = getSortedDistances(centerIdx, distMatrix, n);
                    
                    // 步骤3：从排序后的距离列表中获取基础点集（中心点 + 最近的 n/k 个点）
                    Set<Integer> basePointSet = getBasePointSetFromSorted(centerIdx, sortedDistances, n, k, xiFinal);
                    
                    // 步骤4：计算基础点集的 TSP（只计算一次，对所有点 i 复用）
                    List<Integer> basePointList = new ArrayList<>(basePointSet);
                    double baseObjVal = Double.NaN;
                    if (basePointList.size() >= 2) {
                        baseObjVal = solveTSPMIP(basePointList, distMatrix);
                    } else {
                        baseObjVal = 0.0;
                    }
                    
                    // 【优化】步骤5：从同一个排序后的距离列表中获取最近的5*n/k个点，只对这些点进行循环
                    List<Integer> nearestPoints = getNearestPointsFromSorted(sortedDistances, numPointsPerRegion);
                    
                    // 【优化】先将该区域的所有点初始化
                    Set<Integer> nearestPointsSet = new HashSet<>(nearestPoints);
                    for (int i = 0; i < n; i++) {
                        if (!nearestPointsSet.contains(i)) {
                            if (i == centerIdx) {
                                travelDist[i][j] = 0.0;
                            } else {
                                if (xiFinal != null && xiFinal[i] == 0) {
                                    travelDist[i][j] = Double.NaN;
                                } else {
                                    travelDist[i][j] = OUT_OF_RANGE_MARKER;
                                }
                            }
                        }
                    }
                    
                    if (!nearestPointsSet.contains(centerIdx)) {
                        travelDist[centerIdx][j] = 0.0;
                    }
                    
                    // 对最近的5*n/k个点中的每个点 i，计算包含点 i 的 TSP，并与不包含点 i 的 TSP 做差
                    for (int i : nearestPoints) {
                        int c = completed.incrementAndGet();
                        
                        // 【分析】记录最后200个任务，用于分析哪些区域和点最后完成
                        if (c > total - 200) {
                            String msg1 = String.format("[最后200任务] Region:%d, Point:%d, 进度:%d/%d%n", j, i, c, total);
                            if (loggerFinal != null) loggerFinal.log(msg1);
                            else System.out.print(msg1);
                        }
                        
                        try {
                            if (xiFinal != null && xiFinal[i] == 0 && i != centerIdx) {
                                travelDist[i][j] = Double.NaN;
                                if (c % 100 == 0 || c == total) {
                                    double progress = c * 100.0 / total;
                                    long elapsed = System.currentTimeMillis() - startTime;
                                    long remaining = (c == 0) ? 0 : (long)(elapsed * (double)(total - c) / c);
                                    String msg2 = String.format("进度: %.4f%% (%d/%d) [Region:%d, Point:%d] 已用: %.1fs, 剩余: %.1fs%n", 
                                         progress, c, total, j, i, elapsed/1000.0, remaining/1000.0);
                                    if (loggerFinal != null) loggerFinal.log(msg2);
                                    else System.out.print(msg2);
                                }
                                continue;
                            }
                            
                            List<Integer> pointListWithI = getPointSetWithI(basePointSet, i, centerIdx, xiFinal);
                            double diff = computeDualValue(pointListWithI, basePointSet, i, centerIdx, baseObjVal, distMatrix, xiFinal);
                            travelDist[i][j] = diff;
                        } catch (Exception e) {
                            travelDist[i][j] = Double.NaN;
                        }
                        
                        if (c % 100 == 0 || c == total) {
                            double progress = c * 100.0 / total;
                            long elapsed = System.currentTimeMillis() - startTime;
                            long remaining = (c == 0) ? 0 : (long)(elapsed * (double)(total - c) / c);
                            String msg3 = String.format("进度: %.4f%% (%d/%d) [Region:%d, Point:%d] 已用: %.1fs, 剩余: %.1fs%n", 
                                 progress, c, total, j, i, elapsed/1000.0, remaining/1000.0);
                            if (loggerFinal != null) loggerFinal.log(msg3);
                            else System.out.print(msg3);
                        }
                    }
                    
                } catch (Exception e) {
                    int centerId = centers.get(j).id;
                    int centerIdx = findPointIndexById(allPointsFinal, centerId);
                    if (centerIdx >= 0) {
                        List<Integer> nearestPoints = getNearestPoints(centerIdx, distMatrix, n, numPointsPerRegion);
                        for (int i : nearestPoints) {
                            travelDist[i][j] = Double.NaN;
                            int c = completed.incrementAndGet();
                            if (c % 100 == 0 || c == total) {
                                double progress = c * 100.0 / total;
                                long elapsed = System.currentTimeMillis() - startTime;
                                long remaining = (c == 0) ? 0 : (long)(elapsed * (double)(total - c) / c);
                                String msg4 = String.format("进度: %.4f%% (%d/%d) [Region:%d, Point:%d] (异常) 已用: %.1fs, 剩余: %.1fs%n", 
                                     progress, c, total, j, i, elapsed/1000.0, remaining/1000.0);
                                if (loggerFinal != null) loggerFinal.log(msg4);
                                else System.out.print(msg4);
                            }
                        }
                    }
                } finally {
                    // 【优化】记录区域计算时间
                    long regionTime = System.currentTimeMillis() - regionStartTime;
                    regionTimings.put(j, regionTime);
                    
                    // 如果某个区域计算时间过长，输出警告
                    if (regionTime > 10000) { // 超过10秒
                        String msg5 = String.format("警告: 区域 %d 计算耗时较长: %.2fs%n", j, regionTime / 1000.0);
                        if (loggerFinal != null) loggerFinal.log(msg5);
                        else System.out.print(msg5);
                    }
                }
            });
            
            // 【优化】输出最慢的几个区域，用于分析
            if (regionTimings.size() > 0) {
                List<Map.Entry<Integer, Long>> sortedTimings = new ArrayList<>(regionTimings.entrySet());
                sortedTimings.sort((a, b) -> Long.compare(b.getValue(), a.getValue()));
                String msg6 = "\n最慢的5个区域:\n";
                if (loggerFinal != null) loggerFinal.log(msg6);
                else System.out.print(msg6);
                for (int idx = 0; idx < Math.min(5, sortedTimings.size()); idx++) {
                    Map.Entry<Integer, Long> entry = sortedTimings.get(idx);
                    String msg7 = String.format("  区域 %d: %.2fs%n", entry.getKey(), entry.getValue() / 1000.0);
                    if (loggerFinal != null) loggerFinal.log(msg7);
                    else System.out.print(msg7);
                }
            }
            
            int finalCompleted = completed.get();
            if (finalCompleted != total) {
                double progress = finalCompleted * 100.0 / total;
                long elapsed = System.currentTimeMillis() - startTime;
                String msg8 = String.format("进度: %.4f%% (%d/%d) 已用: %.1fs%n", 
                     progress, finalCompleted, total, elapsed/1000.0);
                if (loggerFinal != null) loggerFinal.log(msg8);
                else System.out.print(msg8);
            } else {
                long elapsed = System.currentTimeMillis() - startTime;
                String msg9 = String.format("进度: 100.0000%% (%d/%d) 已用: %.1fs%n", 
                     finalCompleted, total, elapsed/1000.0);
                if (loggerFinal != null) loggerFinal.log(msg9);
                else System.out.print(msg9);
            }
            
            long totalTime = System.currentTimeMillis() - startTime;
            String msg10 = "求解完成！耗时: " + (totalTime / 1000.0) + "s\n";
            if (loggerFinal != null) loggerFinal.log(msg10);
            else System.out.print(msg10);
            
            return travelDist;
        } finally {
            // 关闭日志记录器
            if (loggerFinal != null) {
                loggerFinal.close();
            }
            // 恢复原始环境
            if (needDisposeEnv && TSPDualExtractor.env != null) {
                try {
                    TSPDualExtractor.env.dispose();
                } catch (GRBException e) {
                    // 忽略错误
                }
            }
            TSPDualExtractor.env = originalEnv;
        }
    }
    
    /**
     * 重载方法：不提供日志文件时，只输出到控制台
     */
    public static double[][] computeTravelDist(List<Point> allPoints, List<Point> centers, 
                                               double[][] distMatrix, int[] xi, GRBEnv env) {
        return computeTravelDist(allPoints, centers, distMatrix, xi, env, null);
    }

    public static void main(String[] args) {
        try {
            System.out.println("=== TSP对偶值提取器 (MIP差值版本) ===");
            
            // 初始化全局环境
            env = new GRBEnv();
            env.set(GRB.IntParam.OutputFlag, 0); // 关闭求解日志
            // 针对多核CPU的优化：外层并行，内层单线程，避免上下文切换开销
            env.set(GRB.IntParam.Threads, 1);
            env.start();

            // 参数设置
            String pointsFile = "data/test/selected_low_ratio_points_top20.csv";
            String centersFile = "data/test/selected_centers_low_ratio_p3.csv";
            String demandMatrixFile = "data/test/demand_matrix_low_ratio.csv";
            String targetDate = "2023-09-03"; // 目标日期：2023年9月3日
            
            // 读取数据
            List<Point> allPoints = readCoordinates(pointsFile);
            List<Point> centers = readCenters(centersFile);
            int n = allPoints.size();
            int k = centers.size();
            
            // 从 demand_matrix.csv 中读取指定日期的列作为 xi 向量
            int[] xi = readXiVectorFromDemandMatrix(demandMatrixFile, targetDate, allPoints);
            
            // 根据中心点个数生成输出文件名
            String outputFile = "output/travel_dist_dual_values_p" + k + ".csv";
            
            System.out.println("N=" + n + ", K=" + k + ". 总任务数: " + (n * k));
            if (xi != null) {
                int accessiblePoints = Arrays.stream(xi).sum();
                System.out.println("使用 xi 向量过滤，可访问的点数: " + accessiblePoints + "/" + n);
            }
            
            // 计算距离矩阵
            double[][] distMatrix = computeDistanceMatrix(allPoints);
            
            // 调用核心计算逻辑
            double[][] travelDist = computeTravelDist(allPoints, centers, distMatrix, xi, env);
            
            // 保存结果
            saveTravelDistToCSV(travelDist, outputFile);
            
            // 清理环境
            env.dispose();
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void saveTravelDistToCSV(double[][] travelDist, String filename) throws IOException {
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "UTF-8"));
        writer.write("\uFEFF"); // BOM
        int n = travelDist.length;
        int k = travelDist[0].length;
        writer.write("PointID");
        for (int j = 0; j < k; j++) writer.write(",Region_" + j);
        writer.write("\n");
        for (int i = 0; i < n; i++) {
            writer.write(String.valueOf(i));
            for (int j = 0; j < k; j++) {
                double value = travelDist[i][j];
                String outputValue;
                if (Double.isNaN(value)) {
                    outputValue = "NaN"; // 计算失败或xi=0等情况
                } else if (value == OUT_OF_RANGE_MARKER) {
                    outputValue = "Null"; // 不在5n/k范围内
                } else {
                    outputValue = String.format("%.10f", value);
                }
                writer.write("," + outputValue);
            }
            writer.write("\n");
        }
        writer.close();
        System.out.println("结果已保存到: " + filename);
    }
    
    /**
     * 格式化基础点集信息（与 TSPDualExtractorValidator 中的格式一致）
     * @param j 区域索引
     * @param basePointSet 基础点集
     * @param centerIdx 中心点索引
     * @param allPoints 所有点的列表
     * @param xi xi 向量
     * @return 格式化的字符串
     */
    private static String formatBasePointSetInfo(int j, Set<Integer> basePointSet, int centerIdx,
                                                 List<Point> allPoints, int[] xi) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("区域 %d 的基础点集:\n", j));
        sb.append(String.format("  点集大小: %d\n", basePointSet.size()));
        sb.append("  点集详情:\n");
        
        // 按索引排序，确保输出顺序一致
        List<Integer> sortedPointSet = new ArrayList<>(basePointSet);
        Collections.sort(sortedPointSet);
        
        for (int pointIdx : sortedPointSet) {
            Point p = allPoints.get(pointIdx);
            boolean isCenter = (pointIdx == centerIdx);
            int xiValue = (xi == null ? 1 : xi[pointIdx]);
            String centerMark = isCenter ? " [中心点/depot]" : "";
            
            sb.append(String.format(
                "    ✓ 点[%d] (经度:%.6f, 纬度:%.6f): ξ[%d]=%d%s\n",
                pointIdx, p.longitude, p.latitude, pointIdx, xiValue, centerMark
            ));
        }
        sb.append("\n");
        
        return sb.toString();
    }
    
    /**
     * 保存基础点集信息到文件
     * @param allBasePointSetInfo 所有区域的基础点集信息列表
     * @param filename 输出文件名
     * @param n 总点数
     * @param k 区域数
     * @param xi xi 向量
     */
    private static void saveBasePointSetInfoToFile(List<String> allBasePointSetInfo, String filename,
                                                   int n, int k, int[] xi) throws IOException {
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "UTF-8"));
        
        // 写入文件头
        writer.write("=== TSP对偶值提取器 - 基础点集信息 ===\n");
        writer.write("生成日期: " + new java.util.Date() + "\n");
        writer.write("\n");
        writer.write("参数设置:\n");
        writer.write("  - 总点数 (N): " + n + "\n");
        writer.write("  - 区域数 (K): " + k + "\n");
        if (xi != null) {
            int accessibleCount = Arrays.stream(xi).sum();
            writer.write("  - 可访问的点数: " + accessibleCount + "/" + n + "\n");
        }
        writer.write("\n");
        writer.write("========================================\n");
        writer.write("所有区域的基础点集（按区域索引排序）:\n");
        writer.write("========================================\n");
        writer.write("\n");
        
        // 按区域索引排序（因为并行处理，顺序可能不一致）
        // 提取区域索引并排序
        List<Pair<Integer, String>> sortedInfo = new ArrayList<>();
        for (String info : allBasePointSetInfo) {
            // 从字符串中提取区域索引
            int regionIdx = extractRegionIndex(info);
            sortedInfo.add(new Pair<>(regionIdx, info));
        }
        sortedInfo.sort(Comparator.comparing(p -> p.first));
        
        // 写入排序后的信息
        for (Pair<Integer, String> pair : sortedInfo) {
            writer.write(pair.second);
        }
        
        writer.write("\n");
        writer.write("========================================\n");
        writer.write("基础点集信息记录完成\n");
        writer.write("========================================\n");
        
        writer.close();
    }
    
    /**
     * 从格式化的字符串中提取区域索引
     */
    private static int extractRegionIndex(String info) {
        // 格式: "区域 32 的基础点集:"
        try {
            int startIdx = info.indexOf("区域 ") + 3;
            int endIdx = info.indexOf(" 的基础点集");
            if (startIdx > 2 && endIdx > startIdx) {
                return Integer.parseInt(info.substring(startIdx, endIdx).trim());
            }
        } catch (Exception e) {
            // 如果解析失败，返回 -1
        }
        return -1;
    }
    
    /**
     * 简单的 Pair 类，用于排序
     */
    private static class Pair<K, V> {
        K first;
        V second;
        Pair(K first, V second) {
            this.first = first;
            this.second = second;
        }
    }
}