import com.gurobi.gurobi.*;
import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.io.StringWriter;

/**
 * TSP对偶值提取器 - 批量处理版本
 * 对 demand_matrix.csv 中的所有日期进行迭代计算
 * 直接调用 TSPDualExtractor 中的算法，不重复实现
 */
public class TSPDualExtractorBatch {
    
    // 特殊值：用于标记不在5n/k范围内的点（在CSV中显示为"Null"）
    private static final double OUT_OF_RANGE_MARKER = Double.NEGATIVE_INFINITY;
    private static final boolean WRITE_TEXT_LOGS = false;
    
    // 日期格式模式：yyyy-MM-dd
    private static final Pattern DATE_PATTERN = Pattern.compile("^\\d{4}-\\d{2}-\\d{2}$");
    
    /**
     * 从 demand_matrix.csv 中读取所有日期列
     * @param filename demand_matrix.csv 文件路径
     * @return 日期列表（格式：yyyy-MM-dd）
     */
    private static List<String> readAllDatesFromDemandMatrix(String filename) throws IOException {
        List<String> dates = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename), "UTF-8"));
        
        // 读取表头
        String headerLine = reader.readLine();
        if (headerLine == null) {
            reader.close();
            return dates;
        }
        
        // 解析表头，提取所有日期列
        String[] headers = headerLine.split(",");
        for (String header : headers) {
            String trimmed = header.trim();
            // 检查是否符合日期格式 yyyy-MM-dd
            if (DATE_PATTERN.matcher(trimmed).matches()) {
                dates.add(trimmed);
            }
        }
        
        reader.close();
        return dates;
    }
    
    /**
     * 保存结果到 CSV 文件
     */
    private static void saveTravelDistToCSV(double[][] travelDist, String filename) throws IOException {
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
                    outputValue = "NaN";
                } else if (value == OUT_OF_RANGE_MARKER) {
                    outputValue = "Null";
                } else {
                    outputValue = String.format("%.10f", value);
                }
                writer.write("," + outputValue);
            }
            writer.write("\n");
        }
        writer.close();
    }
    
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
                new OutputStreamWriter(new FileOutputStream(file, false), "UTF-8")));
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
     * 对单个日期执行完整的计算
     * 【重要】每个日期使用独立的 Gurobi 环境，避免状态累积导致的性能下降
     */
    private static void processSingleDate(String targetDate, 
                                         TSPDualExtractor.Point[] allPointsArray, 
                                         TSPDualExtractor.Point[] centersArray,
                                         double[][] distMatrix, 
                                         String demandMatrixFile, 
                                         int k, 
                                         String outputDir) {
        // 【重要】为每个日期创建独立的 Gurobi 环境，避免状态累积
        GRBEnv dateEnv = null;
        Logger logger = null;
        
        try {
            // 根据日期生成日志文件名
            String dateSuffix = targetDate.replace("-", "");
            File dir = new File(outputDir);
            if (!dir.exists()) {
                dir.mkdirs();
            }
            String logFile = WRITE_TEXT_LOGS
                    ? outputDir + "/log_p" + k + "_" + dateSuffix + ".txt"
                    : "disabled";
            
            // 创建日志记录器
            try {
                if (WRITE_TEXT_LOGS) logger = new Logger(logFile);
            } catch (IOException e) {
                System.err.println("无法创建日志文件: " + logFile + "，将只输出到控制台");
            }
            
            dateEnv = new GRBEnv();
            dateEnv.set(GRB.IntParam.OutputFlag, 0);
            dateEnv.set(GRB.IntParam.Threads, 1);
            dateEnv.start();
            
            String msg = "\n" + "=".repeat(60) + "\n";
            msg += "处理日期: " + targetDate + "\n";
            msg += "=".repeat(60) + "\n";
            if (logger != null) logger.log(msg);
            else System.out.print(msg);
            
            // 将数组转换为列表
            List<TSPDualExtractor.Point> allPoints = Arrays.asList(allPointsArray);
            List<TSPDualExtractor.Point> centers = Arrays.asList(centersArray);
            int n = allPoints.size();
            
            // 从 demand_matrix.csv 中读取指定日期的列作为 xi 向量
            // 注意：readXiVectorFromDemandMatrix 内部使用 System.out.println，无法直接记录到日志
            // 但主要的计算过程会记录到日志文件中
            int[] xi = TSPDualExtractor.readXiVectorFromDemandMatrix(demandMatrixFile, targetDate, allPoints);
            
            if (xi == null) {
                msg = "警告: 无法读取日期 " + targetDate + " 的数据，跳过该日期\n";
                if (logger != null) logger.log(msg);
                else System.out.print(msg);
                if (dateEnv != null) {
                    try {
                        dateEnv.dispose();
                    } catch (GRBException e) {
                        // 忽略错误
                    }
                }
                if (logger != null) logger.close();
                return;
            }
            
            String outputFile = outputDir + "/travel_dist_dual_values_p" + k + "_" + dateSuffix + ".csv";
            
            msg = "N=" + n + ", K=" + k + "\n";
            int accessiblePoints = Arrays.stream(xi).sum();
            msg += "使用 xi 向量过滤，可访问的点数: " + accessiblePoints + "/" + n + "\n";
            msg += "日志文件: " + logFile + "\n";
            if (logger != null) logger.log(msg);
            else System.out.print(msg);
            
            // 调用 TSPDualExtractor 的核心计算逻辑（传入日志文件路径）
            double[][] travelDist = TSPDualExtractor.computeTravelDist(allPoints, centers, distMatrix, xi, dateEnv, null);
            
            // 保存结果
            saveTravelDistToCSV(travelDist, outputFile);
            msg = "结果已保存到: " + outputFile + "\n";
            if (logger != null) logger.log(msg);
            else System.out.print(msg);
            
        } catch (Exception e) {
            String msg = "处理日期 " + targetDate + " 时发生错误: " + e.getMessage() + "\n";
            if (logger != null) {
                logger.log(msg);
                // 将异常堆栈写入日志
                StringWriter sw = new StringWriter();
                PrintWriter pw = new PrintWriter(sw);
                e.printStackTrace(pw);
                logger.log(sw.toString() + "\n");
            } else {
                System.err.print(msg);
                e.printStackTrace();
            }
        } finally {
            // 关闭日志记录器
            if (logger != null) {
                logger.close();
            }
            // 【重要】清理当前日期的环境
            if (dateEnv != null) {
                try {
                    dateEnv.dispose();
                } catch (GRBException e) {
                    System.err.println("清理日期 " + targetDate + " 的 Gurobi 环境时发生错误: " + e.getMessage());
                }
            }
            
            // 强制垃圾回收，清理 Gurobi 模型对象
            System.gc();
        }
    }
    
    // ---------------- 主程序入口 ----------------
    
    public static void main(String[] args) {
        try {
            System.out.println("=== TSP对偶值提取器 - 批量处理版本 ===");
            System.out.println("将对 demand_matrix.csv 中的所有日期进行迭代计算\n");
            
            // 参数设置
            String pointsFile = "data/filtered_top100_active800_lat_le_22_75/unique_coordinates_list.csv";
            String centersFile = "data/filtered_top100_active800_lat_le_22_75/selected_centers_p3.csv";
            String demandMatrixFile = "data/filtered_top100_active800_lat_le_22_75/demand_matrix.csv";
            
            // 创建输出目录
            String outputDir = "output/travel_dist_dual_values_filtered_top100_active800_lat_le_22_75_p3";
            // String pointsFile = "data/test/cluster20_unit_outputs/unique_coordinates_list_cluster20_unit.csv";
            // String centersFile = "data/test/cluster20_unit_outputs/selected_centers_test_p3.csv";
            // String demandMatrixFile = "data/test/cluster20_unit_outputs/demand_matrix_cluster20_dates_filtered.csv";
            
            // // 创建输出目录
            // String outputDir = "output/travel_dist_dual_values_filtered_by_date_low_ratio_cluster20_unit";
            File dir = new File(outputDir);
            if (!dir.exists()) {
                dir.mkdirs();
                System.out.println("已创建输出目录: " + outputDir);
            }
            
            // 读取所有日期
            System.out.println("正在读取 demand_matrix.csv 中的所有日期...");
            List<String> dates = readAllDatesFromDemandMatrix(demandMatrixFile);
            
            if (dates.isEmpty()) {
                System.out.println("错误: 未找到任何日期列，程序退出");
                return;
            }
            
            System.out.println("找到 " + dates.size() + " 个日期:");
            for (int i = 0; i < Math.min(10, dates.size()); i++) {
                System.out.println("  - " + dates.get(i));
            }
            if (dates.size() > 10) {
                System.out.println("  ... 还有 " + (dates.size() - 10) + " 个日期");
            }
            System.out.println();
            
            // 读取数据（只读取一次，所有日期共享）
            System.out.println("正在读取坐标和中心点数据...");
            List<TSPDualExtractor.Point> allPointsList = TSPDualExtractor.readCoordinates(pointsFile);
            List<TSPDualExtractor.Point> centersList = TSPDualExtractor.readCenters(centersFile);
            int n = allPointsList.size();
            int k = centersList.size();
            System.out.println("N=" + n + ", K=" + k + "\n");
            
            // 转换为数组（避免在循环中重复转换）
            TSPDualExtractor.Point[] allPointsArray = allPointsList.toArray(new TSPDualExtractor.Point[0]);
            TSPDualExtractor.Point[] centersArray = centersList.toArray(new TSPDualExtractor.Point[0]);
            
            // 计算距离矩阵（只计算一次，所有日期共享）
            System.out.println("正在计算距离矩阵（所有日期共享）...");
            double[][] distMatrix = TSPDualExtractor.computeDistanceMatrix(allPointsList);
            System.out.println();
            
            // 对每个日期进行处理
            long batchStartTime = System.currentTimeMillis();
            int processedCount = 0;
            int successCount = 0;
            
            for (String date : dates) {
                processedCount++;
                System.out.println("\n[" + processedCount + "/" + dates.size() + "] 开始处理日期: " + date);
                
                try {
                    processSingleDate(date, allPointsArray, centersArray, distMatrix, demandMatrixFile, k, outputDir);
                    successCount++;
                } catch (Exception e) {
                    System.err.println("处理日期 " + date + " 失败: " + e.getMessage());
                    e.printStackTrace();
                }
            }
            
            long batchTotalTime = System.currentTimeMillis() - batchStartTime;
            
            // 总结
            System.out.println("\n" + "=".repeat(60));
            System.out.println("批量处理完成！");
            System.out.println("=".repeat(60));
            System.out.println("总日期数: " + dates.size());
            System.out.println("成功处理: " + successCount);
            System.out.println("失败数量: " + (processedCount - successCount));
            System.out.println("总耗时: " + (batchTotalTime / 1000.0) + "s");
            System.out.println("平均每个日期耗时: " + (batchTotalTime / 1000.0 / processedCount) + "s");
            System.out.println("输出目录: " + outputDir);
            
        } catch (Exception e) {
            System.err.println("批量处理过程中发生错误: " + e.getMessage());
            e.printStackTrace();
        }
    }
}
