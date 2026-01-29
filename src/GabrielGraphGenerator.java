import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Gabriel Graph Generator
 * 
 * 根据论文 "A NEW STATISTICAL APPROACH TO GEOGRAPHIC VARIATION ANALYSIS"
 * 生成Gabriel图。
 * 
 * Gabriel图的定义：
 * 两个地点A和B是相邻的（contiguous），当且仅当所有其他地点都在A-B圆之外，
 * 其中A-B圆是以A和B为直径的圆。
 * 
 * 等价地，如果存在另一个地点C，使得在三角形ABC中，C处的角度≥90度，
 * 则A和B不相邻。
 * 
 * 点坐标生成方式：
 * - 在 [0, 1000] × [0, 1000] 的正方形区域内
 * - 按照均匀分布随机生成
 * - 支持的节点数量：|I| ∈ {100, 500, 1000}
 */
public class GabrielGraphGenerator {
    
    private int numNodes;
    private double[] x;
    private double[] y;
    private boolean[][] adjacencyMatrix;
    private Random random;
    
    /**
     * 构造函数
     * @param numNodes 节点数量
     * @param seed 随机种子（用于可重复性）
     */
    public GabrielGraphGenerator(int numNodes, long seed) {
        this.numNodes = numNodes;
        this.x = new double[numNodes];
        this.y = new double[numNodes];
        this.adjacencyMatrix = new boolean[numNodes][numNodes];
        this.random = new Random(seed);
    }
    
    /**
     * 生成随机点坐标
     * @param minX 最小X坐标
     * @param maxX 最大X坐标
     * @param minY 最小Y坐标
     * @param maxY 最大Y坐标
     */
    public void generateRandomPoints(double minX, double maxX, double minY, double maxY) {
        for (int i = 0; i < numNodes; i++) {
            x[i] = minX + random.nextDouble() * (maxX - minX);
            y[i] = minY + random.nextDouble() * (maxY - minY);
        }
    }
    
    /**
     * 构造Gabriel图
     * 对于每对点(i,j)，检查是否所有其他点都在以i和j为直径的圆外
     */
    public void constructGabrielGraph() {
        System.out.println("正在构造Gabriel图...");
        
        // 初始化邻接矩阵
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                adjacencyMatrix[i][j] = false;
            }
        }
        
        // 对于每对点(i,j)
        for (int i = 0; i < numNodes; i++) {
            for (int j = i + 1; j < numNodes; j++) {
                // 计算以i和j为直径的圆的中心和半径
                double centerX = (x[i] + x[j]) / 2.0;
                double centerY = (y[i] + y[j]) / 2.0;
                double radius = Math.sqrt(Math.pow(x[i] - x[j], 2) + Math.pow(y[i] - y[j], 2)) / 2.0;
                
                // 检查是否有其他点在这个圆内
                boolean isGabrielEdge = true;
                for (int k = 0; k < numNodes; k++) {
                    if (k != i && k != j) {
                        // 计算点k到圆心的距离
                        double distanceToCenter = Math.sqrt(
                            Math.pow(x[k] - centerX, 2) + Math.pow(y[k] - centerY, 2)
                        );
                        
                        // 如果点k在圆内（使用小的epsilon来处理浮点数精度问题）
                        // 注意：如果点在圆周上，我们认为是Gabriel边
                        if (distanceToCenter < radius - 1e-10) {
                            isGabrielEdge = false;
                            break;
                        }
                    }
                }
                
                // 如果没有其他点在圆内，则(i,j)是Gabriel边
                if (isGabrielEdge) {
                    adjacencyMatrix[i][j] = true;
                    adjacencyMatrix[j][i] = true;
                }
            }
        }
        
        System.out.println("Gabriel图构造完成");
    }
    
    /**
     * 获取边的列表
     * @return 边的列表，每条边用两个节点ID表示
     */
    public List<int[]> getEdges() {
        List<int[]> edges = new ArrayList<>();
        for (int i = 0; i < numNodes; i++) {
            for (int j = i + 1; j < numNodes; j++) {
                if (adjacencyMatrix[i][j]) {
                    edges.add(new int[]{i, j});
                }
            }
        }
        return edges;
    }
    
    /**
     * 生成活跃度指标
     * @return 活跃度指标数组，每个节点有3个指标
     */
    public double[][] generateActiveness() {
        double[][] activeness = new double[numNodes][3];
        double baseDemand = 500.0;
        double randomFactor = 0.2;
        
        for (int i = 0; i < numNodes; i++) {
            // 第一个指标：基础需求 + 随机变化
            activeness[i][0] = baseDemand + (random.nextDouble() * 2 * randomFactor - randomFactor) * baseDemand;
            // 第二个指标：1500 + 随机值(0-300)
            activeness[i][1] = 1500 + random.nextDouble() * 300;
            // 第三个指标：1500 + 随机值(0-300)
            activeness[i][2] = 1500 + random.nextDouble() * 300;
        }
        
        return activeness;
    }
    
    /**
     * 将实例数据写入文件（新格式：包含点数、区域数量、点坐标、边信息）
     * @param filename 输出文件名
     * @param k 区域数量
     * @throws IOException 文件写入异常
     */
    public void writeInstanceFile(String filename, int k) throws IOException {
        System.out.println("正在写入文件: " + filename);
        
        FileWriter writer = new FileWriter(filename);
        
        // 写入节点数量
        writer.write(numNodes + "\n");
        
        // 写入区域数量
        writer.write(k + "\n");
        
        // 写入节点数据 (id, x, y) - 只包含序号和坐标
        for (int i = 0; i < numNodes; i++) {
            writer.write(String.format("%d %.6f %.6f\n", i, x[i], y[i]));
        }
        
        // 获取边列表
        List<int[]> edges = getEdges();
        int numEdges = edges.size();
        
        // 写入边的数量
        writer.write(numEdges + "\n");
        
        // 写入边数据
        for (int[] edge : edges) {
            writer.write(edge[0] + " " + edge[1] + "\n");
        }
        
        writer.close();
        System.out.println("文件写入完成");
    }
    
    /**
     * 将实例数据写入文件（旧格式：包含活跃度等额外信息）
     * @param filename 输出文件名
     * @param k 区域数量
     * @throws IOException 文件写入异常
     */
    public void writeInstanceFileOld(String filename, int k) throws IOException {
        System.out.println("正在写入文件: " + filename);
        
        FileWriter writer = new FileWriter(filename);
        
        // 写入节点数量
        writer.write(numNodes + "\n");
        
        // 生成活跃度指标
        double[][] activeness = generateActiveness();
        
        // 写入节点数据 (id, x, y, activeness1, activeness2, activeness3)
        for (int i = 0; i < numNodes; i++) {
            writer.write(String.format("%d %.6f %.6f %.0f %.0f %.0f\n",
                i, x[i], y[i], activeness[i][0], activeness[i][1], activeness[i][2]));
        }
        
        // 获取边列表
        List<int[]> edges = getEdges();
        int numEdges = edges.size();
        
        // 写入边的数量
        writer.write(numEdges + "\n");
        
        // 写入边数据
        for (int[] edge : edges) {
            writer.write(edge[0] + " " + edge[1] + "\n");
        }
        
        // 写入k和其他参数
        writer.write(String.format("%d %.6f %.6f %.6f\n", k, 6.0, 7.0, 0.05));
        
        // 写入随机种子
        writer.write("semilla: " + random.nextInt(Integer.MAX_VALUE) + "\n");
        
        writer.close();
        System.out.println("文件写入完成");
    }
    
    /**
     * 打印图的统计信息
     */
    public void printStatistics() {
        List<int[]> edges = getEdges();
        int numEdges = edges.size();
        
        System.out.println("=== Gabriel图统计信息 ===");
        System.out.println("节点数量: " + numNodes);
        System.out.println("边的数量: " + numEdges);
        System.out.println("平均度数: " + (2.0 * numEdges / numNodes));
        
        // 计算连通性
        boolean[] visited = new boolean[numNodes];
        int components = 0;
        for (int i = 0; i < numNodes; i++) {
            if (!visited[i]) {
                dfs(i, visited);
                components++;
            }
        }
        System.out.println("连通分量数量: " + components);
    }
    
    /**
     * 深度优先搜索，用于计算连通分量
     */
    private void dfs(int node, boolean[] visited) {
        visited[node] = true;
        for (int i = 0; i < numNodes; i++) {
            if (adjacencyMatrix[node][i] && !visited[i]) {
                dfs(i, visited);
            }
        }
    }
    
    /**
     * 主函数 - 用于测试
     */
    public static void main(String[] args) {
        // 默认节点数量：支持 100, 500, 1000
        int numNodes = 100;
        long seed = System.currentTimeMillis();
        int k = 5; // 区域数量
        
        if (args.length >= 1) {
            numNodes = Integer.parseInt(args[0]);
        }
        if (args.length >= 2) {
            seed = Long.parseLong(args[1]);
        }
        if (args.length >= 3) {
            k = Integer.parseInt(args[2]);
        }
        
        GabrielGraphGenerator generator = new GabrielGraphGenerator(numNodes, seed);
        
        // 在[0, 1000] × [0, 1000]的正方形区域内生成随机点（均匀分布）
        generator.generateRandomPoints(0, 1000, 0, 1000);
        
        // 构造Gabriel图
        generator.constructGabrielGraph();
        
        // 打印统计信息
        generator.printStatistics();
        
        // 写入文件（使用新格式：包含点数、区域数量、坐标、边）
        try {
            String filename = String.format("Instances_new/gabriel_%d_%d.dat", numNodes, seed % 10000);
            generator.writeInstanceFile(filename, k);
            System.out.println("实例生成完成！文件保存为: " + filename);
        } catch (IOException e) {
            System.err.println("写入文件时出错: " + e.getMessage());
            e.printStackTrace();
        }
    }
}

