import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Locale;
import java.util.Random;

public class CCTest1 {
    public static void main(String[] args) throws Exception {
        // Experiment configuration
        double E = 50.0; // Expected value
        double[] RSDValues = {0.125, 0.25, 0.5}; // Relative standard deviation array
        double[] rValues = {0.5, 0.4, 0.3, 0.2, 0.1}; // Tolerance parameter values
        double[] gammaValues = {0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01}; // Chance constraint risk parameter
        // 原始 d_i 模式：每个元素为随机场景条数。d_ij+CSV 时场景数以训练目录下 CSV 个数为准，此处数组在 d_ij 下仅首项会触发一次真实实验（避免与 numScenarios 重复跑同一数据）。
        int[] scenarioNumValues = {1000};
        boolean[] useScenarioGeneration = {true}; // Whether to use scenario generation
        /** 为 true 时使用 d_{ijs} 工作量与 T_s 相对平衡模型（与 DR 中 assignment-dependent 口径一致）；false 保持原 d_i 模型。 */
        boolean useDijWorkloadModel = false;
        long seed = 12345678; // Random seed

        long testSeed = seed + 1000;
        int numTestScenarios = 5000;

        // Get start time for filename
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmmss");
        String timestamp = dateFormat.format(new Date());

        // d_ij 模式：结果写入 assignment_dependent 子目录（与仓库下 output/csv/chance_constrained_results/assignment_dependent 对应）
        String outputCSVPath = useDijWorkloadModel
                ? "./output/csv/chance_constrained_results/assignment_dependent/chance_constrained_results_" + timestamp + ".csv"
                : "./output/csv/chance_constrained_results/chance_constrained_results_" + timestamp + ".csv";
        File outputCsvFile = new File(outputCSVPath);
        File outputParent = outputCsvFile.getParentFile();
        if (outputParent != null && !outputParent.exists() && !outputParent.mkdirs()) {
            System.err.println("无法创建输出目录: " + outputParent.getAbsolutePath());
            return;
        }

        System.out.println("CSV文件将保存到: " + outputCsvFile.getAbsolutePath());
        System.out.println("实验开始时间: " + new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));

        // Prepare CSV file
        // Modified code for CCTest1 class - CSV output section

// Prepare CSV file
        try (BufferedWriter csvWriter = new BufferedWriter(new FileWriter(outputCsvFile))) {
            // Write CSV header with modified columns - replacing Instance with NumUnits and NumRegions
            csvWriter.write("InstanceName,NumUnits,NumRegions,RSD,r,gamma,Scenarios,UseScenarioGeneration,UseDijWorkload,Runtime(s),Objective,OutOfSamplePerformance,Status,OptimalAssignment,FinalCenters");
            csvWriter.newLine();
            csvWriter.flush(); // Flush header immediately to ensure it's saved even if program terminates early

            // 实例来源：d_ij 与 DRTest1 / 合成 ad_* CSV 一致用 cluster20；否则仍扫 Instances_new 的 GG 前缀
            ArrayList<File> instanceFileList = new ArrayList<>();
            if (useDijWorkloadModel) {
                File cluster20Dat = new File("./data/test/cluster20_unit_outputs/unique_coordinates_list_cluster20_unit.dat");
                if (!cluster20Dat.exists()) {
                    System.err.println("未找到实例（与 DRTest1 assignment-dependent 分支一致）: " + cluster20Dat.getPath());
                    return;
                }
                instanceFileList.add(cluster20Dat);
                System.out.println("d_ij 模式：使用实例 " + cluster20Dat.getAbsolutePath());
            } else {
                File instDir = new File("./Instances_new");
                File[] instanceFiles = instDir.listFiles((d, name) -> name.endsWith(".dat"));
                if (instanceFiles == null || instanceFiles.length == 0) {
                    System.out.println("No .dat files found in ./Instances_new directory.");
                    return;
                }
                String[] targetPrefixes = {"GG20-2", "GG20-4"};
                System.out.println("只处理文件名以以下前缀开头的实例文件:");
                for (String prefix : targetPrefixes) {
                    System.out.println("  - " + prefix);
                }
                for (File f : instanceFiles) {
                    String fileName = f.getName();
                    boolean ok = false;
                    for (String prefix : targetPrefixes) {
                        if (fileName.startsWith(prefix)) {
                            ok = true;
                            break;
                        }
                    }
                    if (ok) {
                        instanceFileList.add(f);
                    } else {
                        System.out.println("跳过实例 " + fileName + " (不匹配任何目标前缀)");
                    }
                }
                if (instanceFileList.isEmpty()) {
                    System.out.println("Instances_new 中未找到匹配 GG20-/GG50- 的实例。");
                    return;
                }
                System.out.println("将处理 " + instanceFileList.size() + " 个 Instances_new 实例。");
            }

            int cnt = 0;
            for (File instanceFile : instanceFileList) {
                String instanceName = instanceFile.getName();

                System.out.println("处理实例 " + instanceName);

                Instance instance = new Instance(instanceFile.getPath());
                int numUnits = instance.getN();
                int numRegions = instance.k;

                // Iterate through all parameter combinations
                for (double RSD : RSDValues) {
                    for (double r : rValues) {
                        // 与 DRTest1 类似的 gamma 剪枝：固定 (RSD,r) 下，若某 gammaIndex 已非最优/失败，则更大 index 的 gamma 不再真实求解，只写占位行
                        int minGammaIndexToSkip = Integer.MAX_VALUE;
                        for (int gammaIndex = 0; gammaIndex < gammaValues.length; gammaIndex++) {
                            double gamma = gammaValues[gammaIndex];
                            for (int numScenarios : scenarioNumValues) {
                                for (boolean useScenario : useScenarioGeneration) {
                                    if (gammaIndex >= minGammaIndexToSkip) {
                                        writeSkippedChanceConstrainedRow(csvWriter, instanceName, numUnits, numRegions,
                                                RSD, r, gamma, numScenarios, useScenario, useDijWorkloadModel);
                                        continue;
                                    }

                                    int trainCsvCountForDij = -1;
                                    if (useDijWorkloadModel) {
                                        syncAssignmentDependentDirsForRsd(RSD);
                                        trainCsvCountForDij = ChanceConstrainedAlgo.countAssignmentDependentTrainingCsvFiles();
                                        if (trainCsvCountForDij <= 0) {
                                            System.err.println("d_ij: 训练目录无 CSV，跳过: "
                                                    + ChanceConstrainedAlgo.getAssignmentDependentTrainCsvDirectory());
                                            writeSkippedChanceConstrainedRow(csvWriter, instanceName, numUnits, numRegions,
                                                    RSD, r, gamma, numScenarios, useScenario, useDijWorkloadModel);
                                            continue;
                                        }
                                        // 场景数完全由目录决定；scenarioNumValues 若多项则只对首项跑一次，避免重复相同训练集
                                        if (scenarioNumValues.length > 0 && numScenarios != scenarioNumValues[0]) {
                                            continue;
                                        }
                                    }

                                    System.out.println("current inst number:" + cnt + "    >>>>>>>>>>>Running experiment:");
                                    cnt++;
                                    System.out.println("Instance: " + instanceName);
                                    System.out.println("RSD: " + RSD);
                                    System.out.println("r: " + r);
                                    System.out.println("gamma: " + gamma);
                                    if (useDijWorkloadModel) {
                                        System.out.println("Scenarios: " + trainCsvCountForDij + " (训练目录 CSV 数，忽略 scenarioNumValues 中的 " + numScenarios + ")");
                                    } else {
                                        System.out.println("Scenarios: " + numScenarios);
                                    }
                                    System.out.println("Use Scenario Generation: " + useScenario);
                                    System.out.println("Use D_ij workload model: " + useDijWorkloadModel);
                                    if (useDijWorkloadModel) {
                                        System.out.println("d_ij 训练目录: "
                                                + ChanceConstrainedAlgo.getAssignmentDependentTrainCsvDirectory());
                                        System.out.println("d_ij 样本外目录: "
                                                + ChanceConstrainedAlgo.getAssignmentDependentOosCsvDirectory());
                                    }

                                    long startTime = System.currentTimeMillis();

                                    ChanceConstrainedAlgo algo = null;
                                    String errorMessage = null;
                                    if (useDijWorkloadModel) {
                                        try {
                                            algo = ChanceConstrainedAlgo.fromAssignmentDependentTrainingCsv(
                                                    instance, gamma, seed, r);
                                        } catch (IOException e) {
                                            System.err.println("d_ij 加载训练 CSV 失败: " + e.getMessage());
                                            errorMessage = e.getMessage();
                                        }
                                    } else {
                                        double[][] scenarios = generateScenarios(
                                                numUnits, numScenarios, E, RSD, seed);
                                        algo = new ChanceConstrainedAlgo(instance, scenarios, gamma, seed, r);
                                    }

                                    double objectiveValue = 0;
                                    try {
                                        if (errorMessage == null && algo != null) {
                                            objectiveValue = algo.run("", useScenario);
                                        } else {
                                            objectiveValue = -1;
                                        }
                                    } catch (Exception e) {
                                        System.err.println("Error running experiment: " + e.getMessage());
                                        errorMessage = e.getMessage();
                                        objectiveValue = -1;
                                    }
                                    if (objectiveValue == -1 && errorMessage == null) {
                                        System.out.println("Error: Objective value is -1 for instance: " + instanceName);
                                    }

                                    long endTime = System.currentTimeMillis();
                                    double runtime = (endTime - startTime) / 1000.0;

                                    boolean solveSuccess = errorMessage == null
                                            && objectiveValue != -1
                                            && objectiveValue != Double.MAX_VALUE;

                                    double outOfSamplePerformance = -1.0;
                                    if (solveSuccess && algo != null) {
                                        try {
                                            if (useDijWorkloadModel) {
                                                outOfSamplePerformance =
                                                        ChanceConstrainedAlgo.testOutOfSamplePerformanceAssignmentDependentFromConfiguredCsv(
                                                                instance, algo, r);
                                            } else {
                                                outOfSamplePerformance = testOutOfSamplePerformance(
                                                        instance, algo, E, RSD, testSeed, r, numTestScenarios);
                                            }
                                        } catch (Exception e) {
                                            System.err.println("Error testing out-of-sample performance: " + e.getMessage());
                                            outOfSamplePerformance = -1.0;
                                        }
                                    }

                                    int statusCode = algo != null ? algo.getLastModelStatus() : -1;

                                    int scenariosReported = (useDijWorkloadModel && algo != null)
                                            ? algo.getNumScenarios()
                                            : numScenarios;

                                    csvWriter.write(String.format(Locale.US,
                                            "%s,%d,%d,%.3f,%.1f,%.4f,%d,%s,%s,%.3f,%.4f,%.4f,%d,%s,%s",
                                            instanceName, numUnits, numRegions, RSD, r, gamma, scenariosReported,
                                            useScenario ? "true" : "false",
                                            useDijWorkloadModel ? "true" : "false",
                                            runtime,
                                            (objectiveValue == -1 || objectiveValue == Double.MAX_VALUE) ? -1.0 : objectiveValue,
                                            outOfSamplePerformance,
                                            statusCode,
                                            formatOptimalAssignmentForCsv(algo, solveSuccess, objectiveValue, numRegions),
                                            formatFinalCentersForCsv(algo, solveSuccess, objectiveValue)
                                    ));
                                    csvWriter.newLine();
                                    csvWriter.flush();

                                    if (errorMessage != null) {
                                        System.out.println("Warning: Experiment failed with error: " + errorMessage);
                                    } else if (objectiveValue == -1 || objectiveValue == Double.MAX_VALUE) {
                                        System.out.println("Warning: Objective value is invalid for instance: " + instanceName);
                                    }

                                    if (isNonOptimalResult(errorMessage, solveSuccess, objectiveValue)) {
                                        int nextGammaIndex = gammaIndex + 1;
                                        if (nextGammaIndex < minGammaIndexToSkip) {
                                            minGammaIndexToSkip = nextGammaIndex;
                                            System.out.println("Prune enabled: skip gamma indices >= " + nextGammaIndex
                                                    + " for current (RSD, r)=(" + RSD + ", " + r + ").");
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            System.err.println("Error writing results to CSV: " + e.getMessage());
        }

        System.out.println("Experiment completed, results saved to " + outputCsvFile.getAbsolutePath());
    }

    /** 与 DRTest1.isNonOptimalResult 一致：超时、异常或目标无效则视为可触发剪枝的失败。 */
    private static boolean isNonOptimalResult(String errorMessage, boolean solveSuccess, double objectiveValue) {
        return errorMessage != null
                || !solveSuccess
                || objectiveValue == -1
                || objectiveValue == Double.MAX_VALUE;
    }

    /**
     * 剪枝跳过：仍写入一行 CSV，实验参数为计划值，结果列与 DRTest1 占位约定一致（数值 -1，划分列为 "-1"）。
     */
    private static void writeSkippedChanceConstrainedRow(
            BufferedWriter csvWriter,
            String instanceName,
            int numUnits,
            int numRegions,
            double rsd,
            double r,
            double gamma,
            int numScenarios,
            boolean useScenario,
            boolean useDijWorkload) throws IOException {
        csvWriter.write(String.format(Locale.US,
                "%s,%d,%d,%.3f,%.1f,%.4f,%d,%s,%s,%.3f,%.4f,%.4f,%d,%s,%s",
                instanceName, numUnits, numRegions, rsd, r, gamma, numScenarios,
                useScenario ? "true" : "false",
                useDijWorkload ? "true" : "false",
                -1.0, -1.0, -1.0, -1,
                "-1", "-1"));
        csvWriter.newLine();
        csvWriter.flush();
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

    /**
     * 按 RSD 设置 {@link ChanceConstrainedAlgo} 的 assignment-dependent 训练/样本外目录
     *（路径约定与 DRTest1 中 {@code syncAssignmentDependentDirsForRsd} 一致，仅作用于 CC 侧接口）。
     */
    private static void syncAssignmentDependentDirsForRsd(double rsd) {
        if (Math.abs(rsd - 0.125) < 1e-9) {
            ChanceConstrainedAlgo.configureAssignmentDependentCsvDirectories(
                    "data/ad_t125", "data/ad_e125");
        } else if (Math.abs(rsd - 0.25) < 1e-9) {
            ChanceConstrainedAlgo.configureAssignmentDependentCsvDirectories(
                    "data/ad_t25", "data/ad_e25");
        } else {
            ChanceConstrainedAlgo.configureAssignmentDependentCsvDirectories(
                    "data/travel_dist_dual_values_filtered_by_date_cluster20_unit_synth142", null);
        }
    }

    /**
     * 将最优区域划分序列化为 CSV 单列：{@code z0:u1 u2|z1:u3 ...}。
     * 与 DRTest1 中逻辑一致；仅在求解成功且目标值有效时返回非空字符串。
     */
    private static String formatOptimalAssignmentForCsv(
            ChanceConstrainedAlgo algo, boolean solveSuccess, double objectiveValue, int numRegions) {
        if (!solveSuccess || objectiveValue == -1 || objectiveValue == Double.MAX_VALUE || algo == null) {
            return "";
        }
        ArrayList<Integer>[] zones = algo.getZones();
        if (zones == null) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int z = 0; z < numRegions; z++) {
            if (z > 0) {
                sb.append('|');
            }
            sb.append(z).append(':');
            if (z < zones.length && zones[z] != null && !zones[z].isEmpty()) {
                for (int i = 0; i < zones[z].size(); i++) {
                    if (i > 0) {
                        sb.append(' ');
                    }
                    sb.append(zones[z].get(i));
                }
            }
        }
        String raw = sb.toString();
        if (raw.indexOf(',') >= 0 || raw.indexOf('"') >= 0 || raw.indexOf('\n') >= 0 || raw.indexOf('\r') >= 0) {
            return "\"" + raw.replace("\"", "\"\"") + "\"";
        }
        return raw;
    }

    /**
     * 记录最后一次的中心点组合（按当前 centers 顺序）。与 DRTest1 中逻辑一致。
     */
    private static String formatFinalCentersForCsv(
            ChanceConstrainedAlgo algo, boolean solveSuccess, double objectiveValue) {
        if (!solveSuccess || objectiveValue == -1 || objectiveValue == Double.MAX_VALUE || algo == null) {
            return "";
        }
        ArrayList<Area> centers = algo.getCenters();
        if (centers == null || centers.isEmpty()) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < centers.size(); i++) {
            if (i > 0) {
                sb.append(' ');
            }
            sb.append(centers.get(i).getId());
        }
        String raw = sb.toString();
        if (raw.indexOf(',') >= 0 || raw.indexOf('"') >= 0 || raw.indexOf('\n') >= 0 || raw.indexOf('\r') >= 0) {
            return "\"" + raw.replace("\"", "\"\"") + "\"";
        }
        return raw;
    }

    /**
     * Tests the out-of-sample performance of a solution by checking constraint satisfaction
     * across newly generated scenarios not used during optimization.
     *
     * @param instance         The problem instance
     * @param algo             The algorithm with a computed solution
     * @param E                Expected demand value
     * @param RSD              Relative standard deviation for demand generation
     * @param testSeed         Random seed for test scenario generation
     * @param r                Capacity tolerance parameter
     * @param numTestScenarios Number of test scenarios to generate
     * @return The percentage of test scenarios where constraints are satisfied
     *         （相对平衡原始需求：各区需求在 [(1-r)/k,(1+r)/k]×场景总需求 内）
     */
    private static double testOutOfSamplePerformance(
            Instance instance,
            ChanceConstrainedAlgo algo,
            double E,
            double RSD,
            long testSeed,
            double r,
            int numTestScenarios) {

        // Generate test scenarios with a different seed than training
        double[][] testScenarios = generateScenarios(
                instance.getN(), numTestScenarios, E, RSD, testSeed);

        // Get the solution
        ArrayList<Integer>[] zones = algo.getZones();

        // Count satisfied scenarios
        int satisfiedScenarios = 0;

        // For each test scenario
        for (int s = 0; s < numTestScenarios; s++) {
            // Calculate this scenario's total demand
            double scenarioTotalDemand = 0;
            for (int i = 0; i < instance.getN(); i++) {
                scenarioTotalDemand += testScenarios[s][i];
            }

            // 相对平衡（原始需求）：与 ChanceConstrainedAlgo / DRTest1 一致，每区需求 ∈ [(1-r)/k,(1+r)/k]·D_s
            int k = instance.k;
            double lowerCap = (1.0 - r) * (scenarioTotalDemand / k);
            double upperCap = (1.0 + r) * (scenarioTotalDemand / k);

            boolean scenarioSatisfied = true;

            // Check each district
            for (int j = 0; j < zones.length; j++) {
                if (zones[j] == null || zones[j].isEmpty()) {
                    continue;
                }

                // Calculate total demand for this district in this scenario
                double districtDemand = 0;
                for (int areaId : zones[j]) {
                    districtDemand += testScenarios[s][areaId];
                }

                if (districtDemand < lowerCap || districtDemand > upperCap) {
                    scenarioSatisfied = false;
                    break;
                }
            }

            if (scenarioSatisfied) {
                satisfiedScenarios++;
            }
        }

        // Return percentage of satisfied scenarios
        return (double) satisfiedScenarios / numTestScenarios;
    }

}