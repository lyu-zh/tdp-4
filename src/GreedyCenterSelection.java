import java.io.*;
import java.util.*;

/**
 * 基于经纬度坐标的贪心中心点选择算法
 * 使用 Haversine 公式计算地球表面两点间的距离（以米为单位）
 */
public class GreedyCenterSelection {
    
    // 地球半径（米）
    private static final double EARTH_RADIUS_METERS = 6371000.0;
    
    /**
     * 点类，存储经纬度坐标
     */
    static class Point {
        int id;
        double longitude; // 经度
        double latitude;  // 纬度
        
        Point(int id, double longitude, double latitude) {
            this.id = id;
            this.longitude = longitude;
            this.latitude = latitude;
        }
    }
    
    /**
     * 使用 Haversine 公式计算两点间距离（米）
     * @param lon1 点1的经度
     * @param lat1 点1的纬度
     * @param lon2 点2的经度
     * @param lat2 点2的纬度
     * @return 两点间的距离（米）
     */
    private static double haversineDistance(double lon1, double lat1, double lon2, double lat2) {
        // 将角度转换为弧度
        double lat1Rad = Math.toRadians(lat1);
        double lat2Rad = Math.toRadians(lat2);
        double deltaLatRad = Math.toRadians(lat2 - lat1);
        double deltaLonRad = Math.toRadians(lon2 - lon1);
        
        // Haversine 公式
        double a = Math.sin(deltaLatRad / 2) * Math.sin(deltaLatRad / 2) +
                   Math.cos(lat1Rad) * Math.cos(lat2Rad) *
                   Math.sin(deltaLonRad / 2) * Math.sin(deltaLonRad / 2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        
        // 返回距离（米）
        return EARTH_RADIUS_METERS * c;
    }
    
    /**
     * 读取 CSV 文件中的坐标点
     */
    private static List<Point> readCoordinates(String filename) throws IOException {
        List<Point> points = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        
        // 跳过标题行
        String line = reader.readLine();
        
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            
            String[] parts = line.split(",");
            if (parts.length >= 3) {
                try {
                    int id = Integer.parseInt(parts[0].trim());
                    double longitude = Double.parseDouble(parts[1].trim());
                    double latitude = Double.parseDouble(parts[2].trim());
                    points.add(new Point(id, longitude, latitude));
                } catch (NumberFormatException e) {
                    System.err.println("跳过无效行: " + line);
                }
            }
        }
        
        reader.close();
        return points;
    }
    
    /**
     * 贪心选择中心点（与 Algo.java 中的方法相同）
     * 优化版本：按需计算距离，使用缓存避免重复计算
     * @param points 所有点
     * @param k 需要选择的中心点数量
     * @param alpha 阈值参数（与 Algo.java 中相同，初始为0）
     * @param seed 随机种子
     * @return 选中的中心点索引列表
     */
    private static List<Integer> greedySelectCenters(List<Point> points, 
                                                      int k, double alpha, long seed) {
        int n = points.size();
        Random rand = new Random(seed);
        List<Integer> centers = new ArrayList<>();
        boolean[] isCenter = new boolean[n];
        
        // 距离缓存：使用 (i, j) 作为键，其中 i < j
        Map<Long, Double> distanceCache = new HashMap<>();
        
        // 生成缓存键的辅助函数
        java.util.function.BiFunction<Integer, Integer, Long> getCacheKey = (i, j) -> {
            int min = Math.min(i, j);
            int max = Math.max(i, j);
            return ((long) min << 32) | max;
        };
        
        // 获取距离（优先从缓存读取）
        java.util.function.BiFunction<Integer, Integer, Double> getDistance = (i, j) -> {
            if (i == j) return 0.0;
            Long key = getCacheKey.apply(i, j);
            Double cached = distanceCache.get(key);
            if (cached != null) {
                return cached;
            }
            Point p1 = points.get(i);
            Point p2 = points.get(j);
            double dist = haversineDistance(p1.longitude, p1.latitude, 
                                           p2.longitude, p2.latitude);
            distanceCache.put(key, dist);
            return dist;
        };
        
        // 随机选择第一个中心点
        int startId = rand.nextInt(n);
        centers.add(startId);
        isCenter[startId] = true;
        
        // 维护每个点到最近中心点的距离（增量更新优化）
        double[] minDistances = new double[n];
        for (int i = 0; i < n; i++) {
            if (i != startId) {
                minDistances[i] = getDistance.apply(i, startId);
            } else {
                minDistances[i] = 0.0;
            }
        }
        
        System.out.println("开始贪心选择中心点，目标数量: " + k);
        long startTime = System.currentTimeMillis();
        
        // 贪心搜索区域中心点
        while (centers.size() < k) {
            if (centers.size() % 50 == 0) {
                long elapsed = System.currentTimeMillis() - startTime;
                System.out.println("已选择中心点数量: " + centers.size() + "/" + k + 
                                 " (耗时: " + (elapsed / 1000.0) + "秒, 已缓存距离对数: " + 
                                 distanceCache.size() + "个)");
            }
            
            // 优化：使用增量更新，minDistances 已经维护了每个点到最近中心点的距离
            // 不需要重新遍历所有已选中心点，直接使用 minDistances
            
            // 找到最大和最小距离
            double maxMinDist = -1;
            double minMinDist = Double.MAX_VALUE;
            for (int i = 0; i < n; i++) {
                if (!isCenter[i]) {
                    double minDist = minDistances[i];
                    if (minDist > maxMinDist) {
                        maxMinDist = minDist;
                    }
                    if (minDist < minMinDist) {
                        minMinDist = minDist;
                    }
                }
            }
            
            double Thre = maxMinDist - alpha * (maxMinDist - minMinDist);
            
            // 找到所有候选点
            List<Integer> candidates = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                if (!isCenter[i] && minDistances[i] >= Thre) {
                    candidates.add(i);
                }
            }
            
            // 如果候选点为空，选择距离最大的点
            if (candidates.isEmpty()) {
                int bestId = -1;
                double maxDist = -1;
                for (int i = 0; i < n; i++) {
                    if (!isCenter[i] && minDistances[i] > maxDist) {
                        maxDist = minDistances[i];
                        bestId = i;
                    }
                }
                if (bestId != -1) {
                    centers.add(bestId);
                    isCenter[bestId] = true;
                    
                    // 增量更新：更新所有点到新中心点的最小距离
                    minDistances[bestId] = 0.0; // 新中心点到自己的距离为0
                    for (int i = 0; i < n; i++) {
                        if (!isCenter[i]) {
                            double distToNewCenter = getDistance.apply(i, bestId);
                            if (distToNewCenter < minDistances[i]) {
                                minDistances[i] = distToNewCenter;
                            }
                        }
                    }
                } else {
                    break; // 没有更多候选点
                }
            } else {
                // 从候选点中随机选择一个
                int nextId = candidates.get(rand.nextInt(candidates.size()));
                centers.add(nextId);
                isCenter[nextId] = true;
                
                // 增量更新：更新所有点到新中心点的最小距离（关键优化）
                // 只需要与新选的中心点比较，而不是重新遍历所有已选中心点
                minDistances[nextId] = 0.0; // 新中心点到自己的距离为0
                for (int i = 0; i < n; i++) {
                    if (!isCenter[i]) {
                        double distToNewCenter = getDistance.apply(i, nextId);
                        if (distToNewCenter < minDistances[i]) {
                            minDistances[i] = distToNewCenter;
                        }
                    }
                }
            }
        }
        
        System.out.println("中心点选择完成！共选择 " + centers.size() + " 个中心点");
        return centers;
    }
    
    /**
     * 将选中的中心点保存到 CSV 文件
     */
    private static void saveCentersToCSV(List<Point> points, List<Integer> centerIndices, 
                                        String outputFilename) throws IOException {
        // 使用 UTF-8 编码写入，避免中文乱码
        BufferedWriter writer = new BufferedWriter(
            new OutputStreamWriter(new FileOutputStream(outputFilename), "UTF-8"));
        
        // 写入 BOM 标记，确保 Excel 等软件能正确识别 UTF-8 编码
        writer.write("\uFEFF");
        
        // 写入标题行
        writer.write("序号,经度,纬度\n");
        
        // 写入中心点数据
        for (int idx : centerIndices) {
            Point p = points.get(idx);
            writer.write(String.format("%d,%.15f,%.15f\n", p.id, p.longitude, p.latitude));
        }
        
        writer.close();
        System.out.println("中心点已保存到: " + outputFilename);
    }
    
    /**
     * 主函数
     */
    public static void main(String[] args) {
        try {
            // 参数设置
            String inputFile = "data/test/unique_coordinates_list_filtered_new.csv";
            int k = 3; // 需要选择的中心点数量（p值）
            double alpha = 0.0; // 阈值参数（与 Algo.java 中相同）
            long seed = 42; // 随机种子（与 Algo.java 中相同）
            
            // 根据中心点数量 p 动态生成输出文件名
            String outputFile = String.format("data/test/selected_centers_filtered_new_p%d.csv", k);
            
            System.out.println("=== 贪心中心点选择算法 ===");
            System.out.println("输入文件: " + inputFile);
            System.out.println("输出文件: " + outputFile);
            System.out.println("目标中心点数量: " + k);
            System.out.println("随机种子: " + seed);
            System.out.println();
            
            // 读取坐标点
            System.out.println("正在读取坐标点...");
            List<Point> points = readCoordinates(inputFile);
            System.out.println("共读取 " + points.size() + " 个坐标点");
            System.out.println();
            
            // 贪心选择中心点（按需计算距离，无需预先计算整个距离矩阵）
            List<Integer> centerIndices = greedySelectCenters(points, k, alpha, seed);
            System.out.println();
            
            // 保存结果
            saveCentersToCSV(points, centerIndices, outputFile);
            
            System.out.println();
            System.out.println("=== 完成 ===");
            
        } catch (IOException e) {
            System.err.println("错误: " + e.getMessage());
            e.printStackTrace();
        }
    }
}

