import java.io.File;
import java.io.IOException;

/**
 * Gabriel Graph 批量生成器
 * 生成多个 Gabriel Graph 实例并保存到 Instances_new 目录
 */
public class GabrielGraphBatchGenerator {
    
    public static void main(String[] args) {
        // 创建输出目录
        File outputDir = new File("Instances_new");
        if (!outputDir.exists()) {
            outputDir.mkdirs();
            System.out.println("创建目录: " + outputDir.getAbsolutePath());
        }
        
        // 定义要生成的实例参数
        // 规则：|I| = 100: p ∈ {6, 8}
        //       |I| = 500: p ∈ {10, 20}
        //       |I| = 1000: p ∈ {20, 40}
        int[] nodeCounts = {100, 500, 1000}; // 节点数量
        int[][] districtCounts = {
            {6, 8},      // 对应 100 节点
            {10, 20},    // 对应 500 节点
            {20, 40}     // 对应 1000 节点
        };
        
        System.out.println("开始批量生成 Gabriel Graph 实例...");
        System.out.println("输出目录: " + outputDir.getAbsolutePath());
        
        int totalGenerated = 0;
        int instanceCounter = 1;
        
        for (int idx = 0; idx < nodeCounts.length; idx++) {
            int numNodes = nodeCounts[idx];
            int[] pValues = districtCounts[idx];
            
            System.out.println("\n=== 生成 " + numNodes + " 个节点的实例 ===");
            
            for (int p : pValues) {
                // 使用固定的种子模式，确保可重复性
                long seed = (long) numNodes * 10000 + p * 1000 + instanceCounter * 100;
                
                try {
                    GabrielGraphGenerator generator = new GabrielGraphGenerator(numNodes, seed);
                    
                    // 在[0, 1000] × [0, 1000]的正方形区域内生成随机点（均匀分布）
                    generator.generateRandomPoints(0, 1000, 0, 1000);
                    
                    // 构造Gabriel图
                    generator.constructGabrielGraph();
                    
                    // 生成文件名，格式类似：GG100-6-1.dat, GG100-8-1.dat, ...
                    String filename = String.format("Instances_new/GG%d-%d-%d.dat", numNodes, p, instanceCounter);
                    
                    // 写入文件（新格式：包含点数、区域数量、坐标、边）
                    generator.writeInstanceFile(filename, p);
                    
                    totalGenerated++;
                    instanceCounter++;
                    
                    System.out.println("  已生成: " + filename + " (节点数=" + numNodes + ", 区域数=" + p + ")");
                } catch (IOException e) {
                    System.err.println("生成实例 " + numNodes + "-" + p + "-" + instanceCounter + " 时出错: " + e.getMessage());
                    e.printStackTrace();
                }
            }
            
            System.out.println("完成 " + numNodes + " 个节点的实例生成");
        }
        
        System.out.println("\n=== 批量生成完成 ===");
        System.out.println("总共生成 " + totalGenerated + " 个实例");
        System.out.println("文件保存在: " + outputDir.getAbsolutePath());
    }
}

