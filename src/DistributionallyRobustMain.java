import java.util.Random;

public class DistributionallyRobustMain {
    public static void main(String[] args) throws Exception {
        // 指定输入文件和输出文件名
        String instanceFile = "./Instances/2DU60-05-1.dat";  // 使用与ChanceConstrainedMain类似的实例
        String outputFileName = "2DU60-05-1_drcc";

        // 设置分布鲁棒优化参数
        double gamma = 0.1;  // 风险参数（违反概率）
        int numScenarios = 100;  // 场景数量
        double r = 0.1;
        long seed = 12345678;  // 随机种子

        // 需求期望值和相对标准差
        double E = 50.0; // 期望值
        double[] RSDValues = {0.125, 0.25, 0.5}; // 相对标准差数组

        // 选择特定的RSD值
        double RSD = RSDValues[0]; // 可以通过索引选择不同的RSD值

        // 分布鲁棒优化其他参数
        boolean useD1 = true;  // 是否使用D_1模糊集（false则使用D_2模糊集）
        double delta1 = 2;  // D_2模糊集参数
        double delta2 = 4;   // D_2模糊集参数（需保证delta2 > max{delta1, 1}）
        boolean useJointChance = false;  // 是否使用联合机会约束（Bonferroni近似）

        // 加载实例
        Instance instance = new Instance(instanceFile);

        // 生成随机场景 - 使用均匀分布
        double[][] scenarios = generateScenarios(instance.getN(), numScenarios, E, RSD, seed);


        // 默认使用近似方法（MISOCP），可以通过参数修改
        boolean useExactMethod = false; // 设置为true使用精确方法（割约束）
        
        DistributionallyRobustAlgo algo = new DistributionallyRobustAlgo(
                instance, scenarios, gamma, seed, useD1, delta1, delta2, useJointChance, r, useExactMethod);

        // 运行算法并生成结果
        algo.run(outputFileName);

        System.out.println("基于分布鲁棒机会约束的配送区域划分问题求解完成。");

        // 可视化结果
        String outputImagePath = "./output/" + outputFileName + "_visualization.png";
        DistrictVisualizer visualizer = new DistrictVisualizer(instance, algo.getZones(), algo.getCenters());
        visualizer.saveVisualization(outputImagePath);

        System.out.println("可视化图像已保存至: " + outputImagePath);
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

    // 比较不同模型的性能
//    private static void runComparativeAnalysis(Instance instance, double[][] scenarios,
//                                               double gamma, long seed, String baseFileName) throws Exception {
//        System.out.println("\n开始进行比较分析...");
//
//        // 1. 运行确定性模型
//        Algo deterministicAlgo = new Algo(instance);
//        deterministicAlgo.run(baseFileName + "_det");
//
//        // 2. 运行场景近似模型
//        ChanceConstrainedAlgo ccaAlgo = new ChanceConstrainedAlgo(instance, scenarios, gamma, seed,0.1);
//        ccaAlgo.run(baseFileName + "_cc");
//
//        // 3. 运行D_2模糊集分布鲁棒模型
//        DistributionallyRobustAlgo drAlgo = new DistributionallyRobustAlgo(
//                instance, scenarios, gamma, seed, false, 0.05, 0.5, false);
//        drAlgo.run(baseFileName + "_dr_d2");
//
//        // 输出各模型的目标函数值和求解时间等信息
//        System.out.println("比较分析完成。详细结果请参见输出文件。");
//    }
}