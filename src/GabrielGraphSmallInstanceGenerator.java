import java.io.File;
import java.io.IOException;

/**
 * Gabriel Graph 小规模实例生成器
 * 生成150点的Gabriel Graph实例（6个区域和8个区域）
 */
public class GabrielGraphSmallInstanceGenerator {
    
    public static void main(String[] args) {
        // 创建输出目录
        File outputDir = new File("Instances_test");
        if (!outputDir.exists()) {
            outputDir.mkdirs();
            System.out.println("创建目录: " + outputDir.getAbsolutePath());
        }
        
        // 定义要生成的实例参数
        // 150点：中心点数 6 和 8
        int[] nodeCounts = {100}; // 节点数量
        int[][] districtCounts = {
            {6, 6, 8, 8}      // 对应 150 节点elsarticle-template-harv
        };
        
        System.out.println("开始批量生成 Gabriel Graph 小规模实例...");
        System.out.println("输出目录: " + outputDir.getAbsolutePath());
        
        int totalGenerated = 0;
        int instanceCounter = 1;
        
        for (int idx = 0; idx < nodeCounts.length; idx++) {
            int numNodes = nodeCounts[idx];
            int[] pValues = districtCounts[idx];
            
            System.out.println("\n=== 生成 " + numNodes + " 个节点的实例 ===");
            
            for (int p : pValues) {
                // 使用固定的种子模式，确保可重复性
                long seed = (long) numNodes * 9999 + p * 1000 + instanceCounter * 100;
                
                try {
                    GabrielGraphGenerator generator = new GabrielGraphGenerator(numNodes, seed);
                    
                    // 在[0, 1000] × [0, 1000]的正方形区域内生成随机点（均匀分布）
                    generator.generateRandomPoints(0, 1000, 0, 1000);
                    
                    // 构造Gabriel图
                    generator.constructGabrielGraph();
                    
                    // 生成文件名，格式类似：GG20-3-1.dat, GG20-5-1.dat, ...
                    String filename = String.format("Instances_test/GG%d-%d-%d.dat", numNodes, p, instanceCounter);
                    
                    // 写入文件（新格式：包含点数、区域数量、坐标、边）
                    generator.writeInstanceFile(filename, p);
                    
                    // 打印统计信息
                    generator.printStatistics();
                    
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

