import com.gurobi.gurobi.*;
import java.io.*;
import java.util.*;
import java.util.stream.IntStream;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * TSP对偶值提取器验证程序
 * 验证每次加入TSP的点是否满足：
 * 1. x_{ij} = 1（点i被分配给区域j）
 * 2. ξ_i = 1（点i可以被访问）
 */
public class TSPDualExtractorValidator {
    
    private static final double EARTH_RADIUS_METERS = 6371000.0;
    
    /**
     * 点类
     */
    static class Point {
        int id;
        double longitude;
        double latitude;
        
        Point(int id, double longitude, double latitude) {
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
    
    private static List<Point> readCoordinates(String filename) throws IOException {
        List<Point> points = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename), "UTF-8"));
        String line = reader.readLine(); // skip header
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

    private static List<Point> readCenters(String filename) throws IOException {
        return readCoordinates(filename);
    }

    private static int findPointIndexById(List<Point> points, int id) {
        for (int i = 0; i < points.size(); i++) {
            if (points.get(i).id == id) return i;
        }
        return -1;
    }
    
    /**
     * 从 demand_matrix.csv 中读取指定日期的列作为 xi 向量
     */
    private static int[] readXiVectorFromDemandMatrix(String filename, String targetDate, List<Point> allPoints) {
        try {
            int n = allPoints.size();
            int[] xi = new int[n];
            Arrays.fill(xi, 0);
            
            BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename), "UTF-8"));
            
            String headerLine = reader.readLine();
            if (headerLine == null) {
                reader.close();
                System.out.println("警告: demand_matrix.csv 文件为空");
                return null;
            }
            
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
                System.out.println("警告: 在 demand_matrix.csv 中未找到日期 " + targetDate);
                return null;
            }
            
            String line;
            int matchedCount = 0;
            int unmatchedCount = 0;
            double tolerance = 1e-10;  // 【修改】增加容差，从 1e-14 改为 1e-10，提高坐标匹配成功率
            
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                
                String[] parts = line.split(",");
                if (parts.length <= targetDateColIndex) continue;
                
                try {
                    double lon = Double.parseDouble(parts[0].trim());
                    double lat = Double.parseDouble(parts[1].trim());
                    int value = Integer.parseInt(parts[targetDateColIndex].trim());
                    
                    if (value != 0 && value != 1) {
                        value = (value > 0) ? 1 : 0;
                    }
                    
                    boolean found = false;
                    for (int i = 0; i < n; i++) {
                        Point p = allPoints.get(i);
                        if (Math.abs(p.longitude - lon) < tolerance && 
                            Math.abs(p.latitude - lat) < tolerance) {
                            xi[i] = value;
                            matchedCount++;
                            found = true;
                            break;
                        }
                    }
                    
                    // 【新增】如果未找到匹配的点，记录警告（仅对 value = 1 的点）
                    if (!found && value == 1) {
                        unmatchedCount++;
                        if (unmatchedCount <= 10) {  // 只输出前10个未匹配的点，避免输出过多
                            System.out.println("警告: 未找到匹配的点 (经度:" + lon + ", 纬度:" + lat + ", 值:" + value + ")");
                        }
                    }
                } catch (NumberFormatException e) {}
            }
            
            reader.close();
            
            // 【新增】输出匹配统计信息
            int accessibleCount = Arrays.stream(xi).sum();
            System.out.println("xi 向量读取统计:");
            System.out.println("  - 成功匹配的点数: " + matchedCount);
            System.out.println("  - 未匹配的点数: " + unmatchedCount);
            System.out.println("  - 可访问的点数 (ξ_i = 1): " + accessibleCount + "/" + n);
            
            return xi;
        } catch (IOException e) {
            System.out.println("读取 demand_matrix.csv 失败: " + e.getMessage());
            return null;
        }
    }
    
    /**
     * 计算距离矩阵
     */
    private static double[][] computeDistanceMatrix(List<Point> points) {
        int n = points.size();
        double[][] dist = new double[n][n];
        
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
        });
        
        return dist;
    }
    
    /**
     * 计算 x_{ij} 矩阵：对于每个区域j，找到最近的n/k个点，设置 x_{idx,j} = 1
     * 【修改】无论中心点是否可以被访问，都要将其加入点集（中心点充当模型中的depot d）
     * 【重要】在创建基础点集时，只考虑 ξ_i = 1 的点（除了中心点）
     * @param allPoints 所有点的列表
     * @param centers 中心点列表
     * @param dist 距离矩阵
     * @param n 总点数
     * @param k 区域数
     * @param xi 0-1向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     *           如果为 null，则不进行过滤（保持原有行为）
     * @return x_{ij} 矩阵，x[i][j] = 1 表示点 i 属于区域 j
     */
    private static int[][] computeXMatrix(List<Point> allPoints, List<Point> centers, 
                                          double[][] dist, int n, int k, int[] xi) {
        int[][] x = new int[n][k];
        
        for (int j = 0; j < k; j++) {
            int centerId = centers.get(j).id;
            int centerIdx = findPointIndexById(allPoints, centerId);
            if (centerIdx < 0 || centerIdx >= allPoints.size()) continue;
            
            // 【修改】无论中心点是否可以被访问，都要将其加入点集（中心点充当模型中的depot d）
            // 中心点总是包含在内
            x[centerIdx][j] = 1;
            
            // 【重要】收集所有可访问的点（ξ_i = 1），用于计算距离并选择最近的 n/k 个点
            List<AbstractMap.SimpleEntry<Integer, Double>> dists = new ArrayList<>();
            for (int idx = 0; idx < n; idx++) {
                if (idx == centerIdx) {
                    // 跳过中心点（已经添加）
                    continue;
                }
                // 【重要】只考虑 ξ_i = 1 的点，过滤掉 ξ_i = 0 的点
                if (xi != null && xi[idx] == 0) {
                    continue;  // 跳过不可访问的点
                }
                dists.add(new AbstractMap.SimpleEntry<>(idx, dist[centerIdx][idx]));
            }
            
            // 按距离排序，选择最近的 n/k 个可访问的点
            dists.sort(Map.Entry.comparingByValue());
            
            int numBasePoints = n / k;
            for (int idx = 0; idx < numBasePoints && idx < dists.size(); idx++) {
                int pointIdx = dists.get(idx).getKey();
                x[pointIdx][j] = 1;
            }
        }
        
        return x;
    }
    
    /**
     * 获取区域 j 的基础点集（与 TSPDualExtractor 中的逻辑一致）
     * 【修改】无论中心点是否可以被访问，都要将其加入点集（中心点充当模型中的depot d）
     * 【重要】在创建基础点集时，只考虑 ξ_i = 1 的点（除了中心点）
     * @param j 区域索引
     * @param allPoints 所有点的列表
     * @param centers 中心点列表
     * @param distMatrix 距离矩阵
     * @param n 总点数
     * @param k 区域数
     * @param xi 0-1向量，xi[i] = 1 表示点 i 可以被访问，xi[i] = 0 表示点 i 不能被访问
     *           如果为 null，则不进行过滤（保持原有行为）
     * @return 基础点集（中心点 + 最近的 n/k 个可访问的点）
     */
    private static Set<Integer> getBasePointSetForRegion(int j,
                                                          List<Point> allPoints,
                                                          List<Point> centers,
                                                          double[][] distMatrix,
                                                          int n, int k,
                                                          int[] xi) {
        Set<Integer> pointSet = new HashSet<>();
        
        int centerId = centers.get(j).id;
        int centerIdx = findPointIndexById(allPoints, centerId);
        
        if (centerIdx == -1) return pointSet;
        
        // 【修改】无论中心点是否可以被访问，都要将其加入点集（中心点充当模型中的depot d）
        // 中心点总是加入点集
        pointSet.add(centerIdx);
        
        // 【重要】收集所有可访问的点（ξ_i = 1），用于计算距离并选择最近的 n/k 个点
        // 如果中心点不可访问，且没有其他可访问的点，则只返回中心点（点集大小为1，不需要做TSP）
        List<AbstractMap.SimpleEntry<Integer, Double>> distances = new ArrayList<>(n);
        
        for (int idx = 0; idx < n; idx++) {
            if (idx != centerIdx) {
                // 【重要】只考虑 ξ_i = 1 的点，过滤掉 ξ_i = 0 的点
                if (xi != null && xi[idx] == 0) {
                    continue;  // 跳过不可访问的点
                }
                double d = distMatrix[centerIdx][idx];
                distances.add(new AbstractMap.SimpleEntry<>(idx, d));
            }
        }
        
        // 按距离排序，选择最近的 n/k 个可访问的点
        distances.sort(Map.Entry.comparingByValue());
        
        int numBasePoints = n / k;
        for (int idx = 0; idx < numBasePoints && idx < distances.size(); idx++) {
            pointSet.add(distances.get(idx).getKey());
        }
        
        // 注意：如果点集中只有中心点（大小为1），solveTSPMIP 会自动返回 0.0（不需要做TSP）
        
        return pointSet;
    }
    
    /**
     * 获取包含点 i 的点集（与 TSPDualExtractor 中的逻辑一致）
     */
    private static List<Integer> getPointSetWithI(Set<Integer> basePointSet, int i, int[] xi) {
        Set<Integer> pointSet = new HashSet<>();
        
        for (int idx : basePointSet) {
            if (xi == null || xi[idx] == 1) {
                pointSet.add(idx);
            }
        }
        
        if (xi == null || xi[i] == 1) {
            pointSet.add(i);
        }
        
        return new ArrayList<>(pointSet);
    }
    
    /**
     * 验证点集是否满足条件
     * 【修改】允许中心点的 ξ_i = 0（中心点作为 depot 总是被包含）
     * @param pointSet 要验证的点集
     * @param regionIdx 区域索引 j
     * @param centerIdx 中心点（depot）的索引，允许其 ξ_i = 0
     * @param xMatrix x_{ij} 矩阵
     * @param xi xi 向量
     * @param allPoints 所有点的列表
     * @return 验证结果
     */
    private static ValidationResult validatePointSet(Set<Integer> pointSet, int regionIdx, 
                                                     int centerIdx,
                                                     int[][] xMatrix, int[] xi,
                                                     List<Point> allPoints) {
        ValidationResult result = new ValidationResult();
        result.isValid = true;
        result.violations = new ArrayList<>();
        result.details = new ArrayList<>();
        
        result.details.add(String.format("  点集大小: %d", pointSet.size()));
        result.details.add("  验证详情:");
        
        for (int pointIdx : pointSet) {
            // 检查 x_{ij} = 1
            boolean xValid = (xMatrix[pointIdx][regionIdx] == 1);
            
            // 【修改】检查 ξ_i = 1，但中心点（depot）允许 ξ_i = 0
            boolean xiValid;
            if (pointIdx == centerIdx) {
                // 中心点作为 depot，即使 ξ_i = 0 也是有效的
                xiValid = true;
            } else {
                // 非中心点必须满足 ξ_i = 1
                xiValid = (xi == null || xi[pointIdx] == 1);
            }
            
            Point p = allPoints.get(pointIdx);
            String status = (xValid && xiValid) ? "✓" : "✗";
            String detail = String.format(
                "    %s 点[%d] (经度:%.6f, 纬度:%.6f): x[%d][%d]=%d, ξ[%d]=%d%s",
                status, pointIdx, p.longitude, p.latitude, pointIdx, regionIdx, 
                xMatrix[pointIdx][regionIdx], pointIdx, (xi == null ? 1 : xi[pointIdx]),
                (pointIdx == centerIdx) ? " [中心点/depot]" : ""
            );
            result.details.add(detail);
            
            if (!xValid || !xiValid) {
                result.isValid = false;
                String violation = String.format(
                    "点[%d] (经度:%.6f, 纬度:%.6f) 违反条件: x[%d][%d]=%d, ξ[%d]=%d",
                    pointIdx, p.longitude, p.latitude, pointIdx, regionIdx, 
                    xMatrix[pointIdx][regionIdx], pointIdx, (xi == null ? 1 : xi[pointIdx])
                );
                result.violations.add(violation);
            }
        }
        
        result.details.add(String.format("  验证结果: %s", result.isValid ? "通过" : "失败"));
        
        return result;
    }
    
    static class ValidationResult {
        boolean isValid;
        List<String> violations;
        List<String> details; // 详细的验证信息（包括所有点的验证状态）
    }
    
    public static void main(String[] args) {
        try {
            System.out.println("=== TSP对偶值提取器验证程序 ===");
            System.out.println("验证每次加入TSP的点是否满足: x_{ij} = 1 且 ξ_i = 1");
            System.out.println();
            
            // 参数设置（与 TSPDualExtractor 保持一致）
            String pointsFile = "data/unique_coordinates_list.csv";
            String centersFile = "data/selected_centers_p200.csv";
            String demandMatrixFile = "data/demand_matrix.csv";
            String targetDate = "2022-07-17";
            
            // 读取数据
            System.out.println("读取数据...");
            List<Point> allPoints = readCoordinates(pointsFile);
            List<Point> centers = readCenters(centersFile);
            int n = allPoints.size();
            int k = centers.size();
            
            System.out.println("N=" + n + ", K=" + k);
            
            // 读取 xi 向量
            System.out.println("读取 xi 向量...");
            int[] xi = readXiVectorFromDemandMatrix(demandMatrixFile, targetDate, allPoints);
            if (xi != null) {
                int accessibleCount = Arrays.stream(xi).sum();
                System.out.println("可访问的点数: " + accessibleCount + "/" + n);
            } else {
                System.out.println("警告: 未读取到 xi 向量，将只验证 x_{ij} = 1");
            }
            
            // 计算距离矩阵
            System.out.println("计算距离矩阵...");
            double[][] distMatrix = computeDistanceMatrix(allPoints);
            
            // 计算 x_{ij} 矩阵（考虑 ξ_i 过滤，与 TSPDualExtractor 逻辑一致）
            System.out.println("计算 x_{ij} 矩阵（考虑 ξ_i 过滤）...");
            int[][] xMatrix = computeXMatrix(allPoints, centers, distMatrix, n, k, xi);
            
            // 验证
            System.out.println();
            System.out.println("开始验证...");
            System.out.println("验证内容: 对于每个区域 j，验证基础点集中的所有点是否都满足 x_{ij} = 1 且 ξ_i = 1");
            System.out.println("          注意: 中心点（depot）允许 ξ_i = 0，因为中心点作为 depot 总是被包含");
            System.out.println();
            
            AtomicInteger totalChecked = new AtomicInteger(0);
            AtomicInteger totalViolations = new AtomicInteger(0);
            List<String> allResults = Collections.synchronizedList(new ArrayList<>()); // 所有验证结果
            List<String> allViolations = Collections.synchronizedList(new ArrayList<>()); // 违反条件的结果
            
            // 验证每个区域的基础点集
            // 【修改】对于每个区域 j，验证基础点集中的所有点是否都满足 x_{ij} = 1 且 ξ_i = 1（中心点允许 ξ_i = 0）
            IntStream.range(0, k).parallel().forEach(j -> {
                try {
                    // 获取区域 j 的中心点索引
                    int centerId = centers.get(j).id;
                    int centerIdx = findPointIndexById(allPoints, centerId);
                    
                    Set<Integer> basePointSet = getBasePointSetForRegion(j, allPoints, centers, distMatrix, n, k, xi);
                    
                    // 【修改】传入中心点索引，允许中心点的 ξ_i = 0
                    ValidationResult baseResult = validatePointSet(basePointSet, j, centerIdx, xMatrix, xi, allPoints);
                    
                    synchronized (allResults) {
                        totalChecked.incrementAndGet();
                        String header = String.format("区域 %d 的基础点集:", j);
                        allResults.add(header);
                        allResults.addAll(baseResult.details);
                        allResults.add(""); // 空行分隔
                        
                        if (!baseResult.isValid) {
                            totalViolations.incrementAndGet();
                            allViolations.add(header);
                            allViolations.addAll(baseResult.violations);
                        }
                    }
                } catch (Exception e) {
                    System.err.println("验证区域 " + j + " 时出错: " + e.getMessage());
                    e.printStackTrace();
                }
            });
            
            // 输出验证结果
            System.out.println("========================================");
            System.out.println("验证完成！");
            System.out.println("========================================");
            System.out.println("检查的区域数: " + totalChecked.get());
            System.out.println("违反条件的区域数: " + totalViolations.get());
            System.out.println("通过验证的区域数: " + (totalChecked.get() - totalViolations.get()));
            System.out.println();
            
            if (allViolations.isEmpty()) {
                System.out.println("✓ 所有区域的基础点集都满足条件: x_{ij} = 1 且 ξ_i = 1");
            } else {
                System.out.println("✗ 发现违反条件的区域:");
                System.out.println();
                for (String violation : allViolations) {
                    System.out.println(violation);
                }
            }
            
            // 保存验证结果到文件
            String outputFile = "output/validation_result_2.txt";
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputFile), "UTF-8"));
            
            // 写入摘要信息
            writer.write("=== TSP对偶值提取器验证结果 ===\n");
            writer.write("验证日期: " + new java.util.Date() + "\n");
            writer.write("\n");
            writer.write("验证内容: 对于每个区域 j，验证基础点集中的所有点是否都满足 x_{ij} = 1 且 ξ_i = 1\n");
            writer.write("          注意: 中心点（depot）允许 ξ_i = 0，因为中心点作为 depot 总是被包含\n");
            writer.write("\n");
            writer.write("参数设置:\n");
            writer.write("  - 总点数 (N): " + n + "\n");
            writer.write("  - 区域数 (K): " + k + "\n");
            if (xi != null) {
                int accessibleCount = Arrays.stream(xi).sum();
                writer.write("  - 可访问的点数: " + accessibleCount + "/" + n + "\n");
            }
            writer.write("\n");
            writer.write("验证摘要:\n");
            writer.write("  - 检查的区域数: " + totalChecked.get() + "\n");
            writer.write("  - 违反条件的区域数: " + totalViolations.get() + "\n");
            writer.write("  - 通过验证的区域数: " + (totalChecked.get() - totalViolations.get()) + "\n");
            writer.write("\n");
            
            if (allViolations.isEmpty()) {
                writer.write("✓ 所有区域的基础点集都满足条件: x_{ij} = 1 且 ξ_i = 1\n");
            } else {
                writer.write("✗ 发现违反条件的区域:\n");
                writer.write("\n");
                for (String violation : allViolations) {
                    writer.write(violation + "\n");
                }
                writer.write("\n");
            }
            
            writer.write("========================================\n");
            writer.write("完整验证过程（按区域）:\n");
            writer.write("========================================\n");
            writer.write("\n");
            
            // 写入所有验证结果
            for (String result : allResults) {
                writer.write(result + "\n");
            }
            
            writer.write("\n");
            writer.write("========================================\n");
            writer.write("验证完成\n");
            writer.write("========================================\n");
            
            writer.close();
            System.out.println();
            System.out.println("验证结果已保存到: " + outputFile);
            System.out.println("  - 完整验证过程已写入文件");
            System.out.println("  - 文件大小: " + (new File(outputFile).length() / 1024) + " KB");
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
