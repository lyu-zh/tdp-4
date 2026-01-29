import java.util.Random;

public class ChanceConstrainedMain {
    public static void main(String[] args) throws Exception {
        // 指定输入文件
        String instanceFile = "./Instances/DU200-05-10.dat";
        String outputFileName = "DU200-05-10_cc";
        double gamma = 0.05; // 机会约束风险参数
        int numScenarios = 10000; // 场景数量
        long seed = 12345678; // 添加一个固定种子

        // 需求期望值和相对标准差
        double E = 50.0; // 期望值
        double[] RSDValues = {0.125, 0.25, 0.5}; // 相对标准差数组

        // 选择特定的RSD值
        double RSD = RSDValues[0]; // 可以通过索引选择不同的RSD值

        // 加载实例
        Instance instance = new Instance(instanceFile);

        // 生成随机场景 - 使用均匀分布
        double[][] scenarios = generateScenarios(instance.getN(), numScenarios, E, RSD, seed);

        // 创建并运行算法 - 传入固定种子
        ChanceConstrainedAlgo algo = new ChanceConstrainedAlgo(instance, scenarios, gamma, seed,0.1);
        algo.run(outputFileName,false);
        System.out.println("基于机会约束的配送区域划分问题求解完成。");
        String outputImagePath = "./output/" + outputFileName + "_visualization.png";

    }

    // 生成随机场景 - 使用均匀分布
    private static double[][] generateScenarios(int n, int numScenarios, double E, double RSD, long seed) {
        double[][] scenarios = new double[numScenarios][n];
        Random rand = new Random(seed); // 固定种子以保证可重复性

        // 计算均匀分布的左右端点
        double lowerBound = E * (1 - Math.sqrt(3) * RSD);
        double upperBound = E * (1 + Math.sqrt(3) * RSD);

        // 生成场景
        for (int s = 0; s < numScenarios; s++) {
            for (int i = 0; i < n; i++) {
                // 均匀分布生成需求
                double demand = lowerBound + rand.nextDouble() * (upperBound - lowerBound);

                // 确保需求为正值
                scenarios[s][i] = Math.max(1, demand);
            }
        }

        return scenarios;
    }
}