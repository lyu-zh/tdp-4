import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

public class DRTest1 {
    private static class DistrictDiagnostics {
        double minShare;
        double q05Share;
        double medianShare;
        double q95Share;
        double maxShare;
        int lowerViolations;
        int upperViolations;
        double minMargin;
        double medianMargin;
    }

    private static class AssignmentDependentDiagnostics {
        int totalScenarios;
        int satisfiedScenarios;
        double lowerBoundShare;
        double upperBoundShare;
        DistrictDiagnostics[] districtStats;
    }

    public static void main(String[] args) throws Exception {
        // Experiment configuration
        double E = 50.0; // Expected value
        double[] RSDValues = {0.125}; // Relative standard deviation array
        double[] rValues = {0.5, 0.4, 0.3}; // Tolerance parameter values
        double[] gammaValues = {0.4, 0.3, 0.2, 0.1}; // Chance constraint risk parameter
        int[] scenarioNumValues = {1000}; // Number of scenarios
        boolean[] useD1Values = {true}; // Whether to use D_1 or D_2 ambiguity set
        double[] delta1Values = {2}; // D_2 ambiguity set parameter delta1
        double[] delta2Values = {4}; // D_2 ambiguity set parameter delta2
        boolean[] useJointChanceValues = {false}; // Whether to use joint chance constraint
        boolean[] useExactMethodValues = {false}; // Whether to use exact method or approximation
        boolean[] useImprovedModelValues = {false}; // Whether to use improved model (shortest path distance + workload constraint)
        int[] avgDistMethodValues = {3}; // avgDist计算方式：1=方式1（所有j到i的距离平均），2=方式2（j的k近邻平均距离的平均）
        boolean useRelativeBalanceValue = true; // Whether to use relative balance constraint
        boolean useSupportingHyperplaneCutsValue = true; // Whether to use supporting-hyperplane cuts for relative balance (heuristic)
        boolean useAssignmentDependentValue = true; // Whether to use assignment-dependent workload model
        long seed = 12345678; // Random seed

        long testSeed = seed + 1000;
        int numTestScenarios = 5000; // For standard model, number of test scenarios to generate
        int numTrainingScenarios = 100; // For assignment-dependent model, 已弃用（现在使用全部数据作为训练和测试）

        // Get start time for filename
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmmss");
        String timestamp = dateFormat.format(new Date());

        // Output CSV file - save to output/csv folder with timestamp
        // If useImprovedModelValues contains true, save to improved subfolder
        // If useAssignmentDependentValue is true, save to assignment_dependent subfolder
        String outputCSVPath;
        if (useAssignmentDependentValue) {
            outputCSVPath = "./output/csv/assignment_dependent/distributionally_robust_results_" + timestamp + ".csv";
        } else if (useImprovedModelValues.length > 0 && useImprovedModelValues[0]) {
            outputCSVPath = "./output/csv/improved/distributionally_robust_results_" + timestamp + ".csv";
        } else {
            outputCSVPath = "./output/csv/distributionally_robust_results_" + timestamp + ".csv";
        }
        String diagnosticsCSVPath = useAssignmentDependentValue
                ? outputCSVPath.replace(".csv", "_diagnostics.csv")
                : null;
        
        // Create output directory if it doesn't exist
        File outputDir = new File(outputCSVPath).getParentFile();
        if (outputDir != null && !outputDir.exists()) {
            outputDir.mkdirs();
        }

        System.out.println("CSV文件将保存到: " + outputCSVPath);
        System.out.println("实验开始时间: " + new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));

        BufferedWriter diagnosticsWriter = null;
        try (BufferedWriter csvWriter = new BufferedWriter(new FileWriter(outputCSVPath))) {
            if (diagnosticsCSVPath != null) {
                diagnosticsWriter = new BufferedWriter(new FileWriter(diagnosticsCSVPath));
                diagnosticsWriter.write("InstanceName,RSD,r,gamma,UseD1,Delta1,Delta2,UseJointChance,UseExactMethod,Objective,OutOfSamplePerformance,TotalScenarios,SatisfiedScenarios,LowerBoundShare,UpperBoundShare,DistrictIndex,ZoneMembers,MinShare,Q05Share,MedianShare,Q95Share,MaxShare,LowerViolations,UpperViolations,MinMargin,MedianMargin");
                diagnosticsWriter.newLine();
                diagnosticsWriter.flush();
                System.out.println("诊断CSV文件将保存到: " + diagnosticsCSVPath);
            }
            // Write CSV header
            csvWriter.write("InstanceName,NumUnits,NumRegions,RSD,r,gamma,Scenarios,UseD1,Delta1,Delta2,UseJointChance,UseExactMethod,UseImprovedModel,UseAssignmentDependent,AvgDistMethod,AvgDistMean,Runtime(s),Objective,OutOfSamplePerformance,SolveSuccess,StatusCode,CutIterations,FailureStage");
            csvWriter.newLine();
            csvWriter.flush(); // Flush header immediately to ensure it's saved even if program terminates early

            // Get instance files based on model type
            ArrayList<File> instanceFiles = new ArrayList<>();
            if (useAssignmentDependentValue) {
                // Assignment-dependent model: use cluster20 unit instance data
                File instanceFile = new File("./data/test/cluster20_unit_outputs/unique_coordinates_list_cluster20_unit.dat");
                if (instanceFile.exists()) {
                    instanceFiles.add(instanceFile);
                    System.out.println("找到assignment-dependent模型数据文件: " + instanceFile.getName());
                } else {
                    System.out.println("未找到assignment-dependent模型数据文件: ./data/test/cluster20_unit_outputs/unique_coordinates_list_cluster20_unit.dat");
                    return;
                }
            } else {
                // Standard model: use Instances_new directory
                File dir = new File("./Instances_new");
                File[] allFiles = dir.listFiles((d, name) -> name.endsWith(".dat"));

                // If no instance files are found, provide a hint
                if (allFiles == null || allFiles.length == 0) {
                    System.out.println("No .dat files found in ./Instances_new directory.");
                    return;
                }

                // Filter to only test 20-point and 50-point instances (GG20- and GG50-)
                for (File file : allFiles) {
                    String fileName = file.getName();
                    // 匹配GG20-和GG50-开头的文件
                    if (fileName.startsWith("GG50-")) {
                        instanceFiles.add(file);
                        System.out.println("找到数据文件: " + fileName);
                    }
                }

                if (instanceFiles.isEmpty()) {
                    System.out.println("未找到20点的数据文件（GG20-开头）");
                    return;
                }
            }

            System.out.println("将测试 " + instanceFiles.size() + " 个数据文件");
            int cnt = 0;
            // Iterate through all instance files
            for (File instanceFile : instanceFiles) {
                String instanceName = instanceFile.getName();

                // Iterate through all parameter combinations
                for (double RSD : RSDValues) {
                    for (double r : rValues) {
                        for (double gamma : gammaValues) {
                            for (int numScenarios : scenarioNumValues) {
                                for (boolean useD1 : useD1Values) {
                                    // For D_1, delta values are not used, but we still need to iterate
                                    // For D_2, iterate through delta values
                                    if (useD1) {
                                        // Use D_1: delta values are not used, set to 0
                                        for (boolean useJointChance : useJointChanceValues) {
                                            for (boolean useExactMethod : useExactMethodValues) {
                                                for (boolean useImprovedModel : useImprovedModelValues) {
                                                        for (int avgDistMethod : avgDistMethodValues) {
                                                        // 只有在使用改进模型时才测试avgDistMethod
                                                        // if (!useImprovedModel && avgDistMethod != 1) continue;
                                                        // Print current experiment information
                                                        System.out.println("current inst number:" + cnt + "    >>>>>>>>>>>Running experiment:");
                                                        cnt++;
                                                        System.out.println("Instance: " + instanceName);
                                                        System.out.println("RSD: " + RSD);
                                                        System.out.println("r: " + r);
                                                        System.out.println("gamma: " + gamma);
                                                        System.out.println("Scenarios: " + numScenarios);
                                                        System.out.println("Use D_1: " + useD1);
                                                        System.out.println("Use Joint Chance: " + useJointChance);
                                                        System.out.println("Use Exact Method: " + useExactMethod);
                                                        System.out.println("Use Improved Model: " + useImprovedModel);
                                                        System.out.println("AvgDist Method: " + avgDistMethod);
                                                        System.out.println("Use Assignment Dependent: " + useAssignmentDependentValue);

                                                        // Generate random scenarios (only for standard model)
                                                        Instance instance = new Instance(instanceFile.getPath());

                                                        // Get the number of basic units and regions from the instance
                                                        int numUnits = instance.getN();
                                                        int numRegions = instance.k;

                                                        double[][] scenarios = null;
                                                        if (!useAssignmentDependentValue) {
                                                            scenarios = generateScenarios(
                                                                    numUnits, numScenarios, E, RSD, seed
                                                            );
                                                        }

                                                        // Record start time
                                                        long startTime = System.currentTimeMillis();

                                                        // Create algorithm instance
                                                        DistributionallyRobustAlgo algo;
                                                        if (useAssignmentDependentValue) {
                                                            // Assignment-dependent model: use new constructor
                                                            algo = new DistributionallyRobustAlgo(
                                                                    instance, gamma, seed, useD1, 0.0, 0.0, useJointChance, r, useExactMethod, useRelativeBalanceValue, useSupportingHyperplaneCutsValue
                                                            );
                                                        } else {
                                                            // Standard model: use original constructor
                                                            algo = new DistributionallyRobustAlgo(
                                                                    instance, scenarios, gamma, seed, useD1, 0.0, 0.0, useJointChance, r, useExactMethod, useImprovedModel, avgDistMethod, useRelativeBalanceValue, useSupportingHyperplaneCutsValue
                                                            );
                                                        }

                                                        // Run algorithm and get objective value
                                                        double objectiveValue = 0;
                                                        boolean solveSuccess = false;
                                                        int statusCode = -1;
                                                        int cutIterations = -1;
                                                        int failureStage = -1;
                                                        String errorMessage = null;
                                                        try {
                                                            algo.run(instanceName);
                                                            objectiveValue = algo.getBestObjective();
                                                            solveSuccess = algo.isSolveSuccess();
                                                            statusCode = algo.getStatusCode();
                                                            cutIterations = algo.getCutIterations();
                                                            failureStage = algo.getFailureStage();
                                                        } catch (Exception e) {
                                                            System.err.println("Error running experiment: " + e.getMessage());
                                                            errorMessage = e.getMessage();
                                                            solveSuccess = false;
                                                            objectiveValue = -1;
                                                            failureStage = -1; // 异常情况下设为-1表示未知
                                                        }

                                                        // Calculate runtime
                                                        long endTime = System.currentTimeMillis();
                                                        double runtime = (endTime - startTime) / 1000.0;

                                                        // Test out-of-sample performance (only if algorithm ran successfully)
                                                        double outOfSamplePerformance = -1.0;
                                                        if (errorMessage == null && objectiveValue != -1 && objectiveValue != Double.MAX_VALUE) {
                                                            try {
                                                                if (useAssignmentDependentValue) {
                                                                    outOfSamplePerformance = testOutOfSamplePerformanceAssignmentDependent(
                                                                            instance, algo, r, numTrainingScenarios);
                                                                } else {
                                                                    outOfSamplePerformance = testOutOfSamplePerformance(
                                                                            instance, algo, E, RSD, testSeed, r, numTestScenarios);
                                                                }
                                                            } catch (Exception e) {
                                                                System.err.println("Error testing out-of-sample performance: " + e.getMessage());
                                                                outOfSamplePerformance = -1.0;
                                                            }
                                                        }

                                                        if (useAssignmentDependentValue && diagnosticsWriter != null
                                                                && errorMessage == null && objectiveValue != -1
                                                                && objectiveValue != Double.MAX_VALUE) {
                                                            try {
                                                                AssignmentDependentDiagnostics diagnostics =
                                                                        collectAssignmentDependentDiagnostics(instance, algo, r);
                                                                writeAssignmentDependentDiagnostics(
                                                                        diagnosticsWriter, diagnostics,
                                                                        instanceName, RSD, r, gamma, useD1,
                                                                        0.0, 0.0, useJointChance, useExactMethod,
                                                                        objectiveValue, outOfSamplePerformance, algo);
                                                            } catch (Exception e) {
                                                                System.err.println("Error writing assignment-dependent diagnostics: " + e.getMessage());
                                                            }
                                                        }

                                                        // Get avgDistMean
                                                        double avgDistMean = algo.getAvgDistMean();
                                                        
                                                        // Write to CSV file (record all experiments, including failures)
                                                        csvWriter.write(String.format(
                                                                "%s,%d,%d,%.3f,%.1f,%.1f,%d,%s,%.2f,%.2f,%s,%s,%s,%s,%d,%.4f,%.3f,%.4f,%.4f,%s,%d,%d,%d",
                                                                instanceName, numUnits, numRegions, RSD, r, gamma, numScenarios,
                                                                useD1 ? "true" : "false",
                                                                0.0, 0.0, // Delta values not used for D_1
                                                                useJointChance ? "true" : "false",
                                                                useExactMethod ? "true" : "false",
                                                                useImprovedModel ? "true" : "false",
                                                                useAssignmentDependentValue ? "true" : "false",
                                                                avgDistMethod,
                                                                (avgDistMean < 0) ? -1.0 : avgDistMean,
                                                                runtime,
                                                                (objectiveValue == -1 || objectiveValue == Double.MAX_VALUE) ? -1.0 : objectiveValue,
                                                                outOfSamplePerformance,
                                                                solveSuccess ? "true" : "false",
                                                                statusCode,
                                                                cutIterations,
                                                                failureStage
                                                        ));
                                                        csvWriter.newLine();

                                                        // Flush to ensure real-time writing
                                                        csvWriter.flush();

                                                        // Print warning if experiment failed
                                                        if (errorMessage != null) {
                                                            System.out.println("Warning: Experiment failed with error: " + errorMessage);
                                                        } else if (objectiveValue == -1 || objectiveValue == Double.MAX_VALUE) {
                                                            System.out.println("Warning: Objective value is invalid for instance: " + instanceName);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    } else {
                                        // Use D_2: iterate through delta values
                                        for (double delta1 : delta1Values) {
                                            for (double delta2 : delta2Values) {
                                                if (delta2 <= delta1) continue; // delta2 must be greater than delta1
                                                for (boolean useJointChance : useJointChanceValues) {
                                                    for (boolean useExactMethod : useExactMethodValues) {
                                                        for (boolean useImprovedModel : useImprovedModelValues) {
                                                                for (int avgDistMethod : avgDistMethodValues) {
                                                                // 只有在使用改进模型时才测试avgDistMethod
                                                                // if (!useImprovedModel && avgDistMethod != 1) continue;
                                                                // Print current experiment information
                                                                System.out.println("current inst number:" + cnt + "    >>>>>>>>>>>Running experiment:");
                                                                cnt++;
                                                                System.out.println("Instance: " + instanceName);
                                                                System.out.println("RSD: " + RSD);
                                                                System.out.println("r: " + r);
                                                                System.out.println("gamma: " + gamma);
                                                                System.out.println("Scenarios: " + numScenarios);
                                                                System.out.println("Use D_1: " + useD1);
                                                                System.out.println("Delta1: " + delta1);
                                                                System.out.println("Delta2: " + delta2);
                                                                System.out.println("Use Joint Chance: " + useJointChance);
                                                                System.out.println("Use Exact Method: " + useExactMethod);
                                                                System.out.println("Use Improved Model: " + useImprovedModel);
                                                                System.out.println("AvgDist Method: " + avgDistMethod);
                                                                System.out.println("Use Assignment Dependent: " + useAssignmentDependentValue);

                                                                // Generate random scenarios (only for standard model)
                                                                Instance instance = new Instance(instanceFile.getPath());

                                                                // Get the number of basic units and regions from the instance
                                                                int numUnits = instance.getN();
                                                                int numRegions = instance.k;

                                                                double[][] scenarios = null;
                                                                if (!useAssignmentDependentValue) {
                                                                    scenarios = generateScenarios(
                                                                            numUnits, numScenarios, E, RSD, seed
                                                                    );
                                                                }

                                                                // Record start time
                                                                long startTime = System.currentTimeMillis();

                                                                // Create algorithm instance
                                                                DistributionallyRobustAlgo algo;
                                                                if (useAssignmentDependentValue) {
                                                                    // Assignment-dependent model: use new constructor
                                                                    algo = new DistributionallyRobustAlgo(
                                                                            instance, gamma, seed, useD1, delta1, delta2, useJointChance, r, useExactMethod, useRelativeBalanceValue, useSupportingHyperplaneCutsValue
                                                                    );
                                                                } else {
                                                                    // Standard model: use original constructor
                                                                    algo = new DistributionallyRobustAlgo(
                                                                            instance, scenarios, gamma, seed, useD1, delta1, delta2, useJointChance, r, useExactMethod, useImprovedModel, avgDistMethod, useRelativeBalanceValue, useSupportingHyperplaneCutsValue
                                                                    );
                                                                }

                                                                // Run algorithm and get objective value
                                                                double objectiveValue = 0;
                                                                boolean solveSuccess = false;
                                                                int statusCode = -1;
                                                                int cutIterations = -1;
                                                                int failureStage = -1;
                                                                String errorMessage = null;
                                                                try {
                                                                    algo.run(instanceName);
                                                                    objectiveValue = algo.getBestObjective();
                                                                    solveSuccess = algo.isSolveSuccess();
                                                                    statusCode = algo.getStatusCode();
                                                                    cutIterations = algo.getCutIterations();
                                                                    failureStage = algo.getFailureStage();
                                                                } catch (Exception e) {
                                                                    System.err.println("Error running experiment: " + e.getMessage());
                                                                    errorMessage = e.getMessage();
                                                                    solveSuccess = false;
                                                                    objectiveValue = -1;
                                                                    failureStage = -1; // 异常情况下设为-1表示未知
                                                                }

                                                                // Calculate runtime
                                                                long endTime = System.currentTimeMillis();
                                                                double runtime = (endTime - startTime) / 1000.0;

                                                                // Test out-of-sample performance (only if algorithm ran successfully)
                                                                double outOfSamplePerformance = -1.0;
                                                                if (errorMessage == null && objectiveValue != -1 && objectiveValue != Double.MAX_VALUE) {
                                                                    try {
                                                                        if (useAssignmentDependentValue) {
                                                                            outOfSamplePerformance = testOutOfSamplePerformanceAssignmentDependent(
                                                                                    instance, algo, r, numTrainingScenarios);
                                                                        } else {
                                                                            outOfSamplePerformance = testOutOfSamplePerformance(
                                                                                    instance, algo, E, RSD, testSeed, r, numTestScenarios);
                                                                        }
                                                                    } catch (Exception e) {
                                                                        System.err.println("Error testing out-of-sample performance: " + e.getMessage());
                                                                        outOfSamplePerformance = -1.0;
                                                                    }
                                                                }

                                                                if (useAssignmentDependentValue && diagnosticsWriter != null
                                                                        && errorMessage == null && objectiveValue != -1
                                                                        && objectiveValue != Double.MAX_VALUE) {
                                                                    try {
                                                                        AssignmentDependentDiagnostics diagnostics =
                                                                                collectAssignmentDependentDiagnostics(instance, algo, r);
                                                                        writeAssignmentDependentDiagnostics(
                                                                                diagnosticsWriter, diagnostics,
                                                                                instanceName, RSD, r, gamma, useD1,
                                                                                delta1, delta2, useJointChance, useExactMethod,
                                                                                objectiveValue, outOfSamplePerformance, algo);
                                                                    } catch (Exception e) {
                                                                        System.err.println("Error writing assignment-dependent diagnostics: " + e.getMessage());
                                                                    }
                                                                }

                                                                // Get avgDistMean
                                                                double avgDistMean = algo.getAvgDistMean();
                                                                
                                                                // Write to CSV file (record all experiments, including failures)
                                                                csvWriter.write(String.format(
                                                                        "%s,%d,%d,%.3f,%.1f,%.1f,%d,%s,%.2f,%.2f,%s,%s,%s,%s,%d,%.4f,%.3f,%.4f,%.4f,%s,%d,%d,%d",
                                                                        instanceName, numUnits, numRegions, RSD, r, gamma, numScenarios,
                                                                        useD1 ? "true" : "false",
                                                                        delta1, delta2,
                                                                        useJointChance ? "true" : "false",
                                                                        useExactMethod ? "true" : "false",
                                                                        useImprovedModel ? "true" : "false",
                                                                        useAssignmentDependentValue ? "true" : "false",
                                                                        avgDistMethod,
                                                                        (avgDistMean < 0) ? -1.0 : avgDistMean,
                                                                        runtime,
                                                                        (objectiveValue == -1 || objectiveValue == Double.MAX_VALUE) ? -1.0 : objectiveValue,
                                                                        outOfSamplePerformance,
                                                                        solveSuccess ? "true" : "false",
                                                                        statusCode,
                                                                        cutIterations,
                                                                        failureStage
                                                                ));
                                                                csvWriter.newLine();

                                                                // Flush to ensure real-time writing
                                                                csvWriter.flush();

                                                                // Print warning if experiment failed
                                                                if (errorMessage != null) {
                                                                    System.out.println("Warning: Experiment failed with error: " + errorMessage);
                                                                } else if (objectiveValue == -1 || objectiveValue == Double.MAX_VALUE) {
                                                                    System.out.println("Warning: Objective value is invalid for instance: " + instanceName);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
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
        } finally {
            if (diagnosticsWriter != null) {
                try {
                    diagnosticsWriter.close();
                } catch (IOException ignored) {
                }
            }
        }

        System.out.println("Experiment completed, results saved to " + outputCSVPath);
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
     * Tests the out-of-sample performance of a solution by checking constraint satisfaction
     * across newly generated scenarios not used during optimization.
     *
     * @param instance         The problem instance
     * @param algo             The algorithm with a computed solution
     * @param E                Expected demand value
     * @param RSD              Relative standard deviation for demand generation
     * @param testSeed          Random seed for test scenario generation
     * @param r                 Tolerance parameter
     * @param numTestScenarios  Number of test scenarios to generate
     * @return The percentage of test scenarios where constraints are satisfied
     */
    private static double testOutOfSamplePerformance(
            Instance instance,
            DistributionallyRobustAlgo algo,
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
        ArrayList<Area> centers = algo.getCenters();

        // Check if using improved model
        boolean useImprovedModel = algo.getUseImprovedModel();
        double[][] shortestPathDist = algo.getShortestPathDist();
        double[] meanVector = algo.getMeanVector();
        boolean useRelativeBalance = algo.getUseRelativeBalance();

        // 获取主算法中调整后的最终上下界（仅用于绝对平衡性约束）
        double demandLowerBound = algo.getDemandLowerBound();
        double demandUpperBound = algo.getDemandUpperBound();

        int p = zones.length; // 区域数量
        double coeffLower = (1.0 - r) / p; // (1-r)/p
        double coeffUpper = (1.0 + r) / p; // (1+r)/p

        // Count satisfied scenarios
        int satisfiedScenarios = 0;

        // For each test scenario
        for (int s = 0; s < numTestScenarios; s++) {
            boolean scenarioSatisfied = true;

            if (useRelativeBalance) {
                // 相对平衡性约束：对于每个场景，计算场景特定的上下界
                if (useImprovedModel && shortestPathDist != null && meanVector != null) {
                    // 改进模型：使用工作量约束
                    // 计算该场景的总工作量：sum_k sum_j d_j^s * r_jk
                    double totalWorkload = 0.0;
                    for (int k = 0; k < zones.length; k++) {
                        if (zones[k] == null || zones[k].isEmpty()) {
                            continue;
                        }
                        int centerId_k = centers.get(k).getId();
                        for (int areaId : zones[k]) {
                            double r_jk = shortestPathDist[areaId][centerId_k];
                            totalWorkload += testScenarios[s][areaId] * r_jk;
                        }
                    }

                    // Check each district
                    for (int j = 0; j < zones.length; j++) {
                        if (zones[j] == null || zones[j].isEmpty()) {
                            continue;
                        }

                        int centerId = centers.get(j).getId();
                        // 计算区域 j 的工作量：∑_{i ∈ zone_j} d_i^s * r_ij
                        double districtWorkload = 0.0;
                        for (int areaId : zones[j]) {
                            double r_ij = shortestPathDist[areaId][centerId];
                            districtWorkload += testScenarios[s][areaId] * r_ij;
                        }

                        // 相对平衡性约束：(1-r)/p * totalWorkload <= districtWorkload <= (1+r)/p * totalWorkload
                        double lowerBound = coeffLower * totalWorkload;
                        double upperBound = coeffUpper * totalWorkload;
                        if (districtWorkload < lowerBound || districtWorkload > upperBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                } else {
                    // 原始模型：使用需求约束
                    // 计算该场景的总需求：sum_k sum_j d_j^s
                    double totalDemand = 0.0;
                    for (int i = 0; i < instance.getN(); i++) {
                        totalDemand += testScenarios[s][i];
                    }

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

                        // 相对平衡性约束：(1-r)/p * totalDemand <= districtDemand <= (1+r)/p * totalDemand
                        double lowerBound = coeffLower * totalDemand;
                        double upperBound = coeffUpper * totalDemand;
                        if (districtDemand < lowerBound || districtDemand > upperBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                }
            } else {
                // 绝对平衡性约束：使用固定的上下界
                if (useImprovedModel && shortestPathDist != null && meanVector != null) {
                    // 改进模型：使用工作量约束
                    // 使用主算法中调整后的最终上下界（而不是为每个场景重新计算）

                    // Check each district
                    for (int j = 0; j < zones.length; j++) {
                        if (zones[j] == null || zones[j].isEmpty()) {
                            continue;
                        }

                        int centerId = centers.get(j).getId();
                        // 计算区域 j 的工作量：∑_{i ∈ zone_j} d_i^s * r_ij
                        double districtWorkload = 0.0;
                        for (int areaId : zones[j]) {
                            double r_ij = shortestPathDist[areaId][centerId];
                            districtWorkload += testScenarios[s][areaId] * r_ij;
                        }

                        // 使用主算法中调整后的最终上下界检查约束
                        if (districtWorkload > demandUpperBound || districtWorkload < demandLowerBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                } else {
                    // 原始模型：使用需求约束
                    // 使用主算法中调整后的最终上下界（而不是为每个场景重新计算）

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

                        // 使用主算法中调整后的最终上下界检查约束
                        if (districtDemand > demandUpperBound || districtDemand < demandLowerBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                }
            }

            if (scenarioSatisfied) {
                satisfiedScenarios++;
            }
        }

        // Calculate and return percentage of satisfied scenarios
        double outOfSamplePerformance = (double) satisfiedScenarios / numTestScenarios;

        // Output to console
        System.out.println(String.format("样本外性能: %.4f (%d/%d 场景满足约束)",
                outOfSamplePerformance, satisfiedScenarios, numTestScenarios));

        return outOfSamplePerformance;
    }

    /**
     * Tests the out-of-sample performance for assignment-dependent model
     * Uses all CSV files as test data (与训练数据相同，使用全部数据)
     *
     * @param instance         The problem instance
     * @param algo             The algorithm with a computed solution
     * @param r                 Tolerance parameter
     * @param numTrainingScenarios  已弃用参数（保留以保持接口兼容性）
     * @return The percentage of test scenarios where constraints are satisfied
     */
    private static double testOutOfSamplePerformanceAssignmentDependent(
            Instance instance,
            DistributionallyRobustAlgo algo,
            double r,
            int numTrainingScenarios) {

        // Get all CSV files from the directory (使用与训练数据相同的目录)
        String dataDir = "data/travel_dist_dual_values_filtered_by_date_cluster20_unit_filled";
        File dir = new File(dataDir);
        File[] allFiles = dir.listFiles((d, name) -> name.endsWith(".csv") && name.startsWith("travel_dist_dual_values_p3"));
        
        if (allFiles == null || allFiles.length == 0) {
            System.err.println("错误: 在目录 " + dataDir + " 中未找到CSV文件");
            return -1.0;
        }
        
        // Sort files to ensure consistent ordering
        java.util.Arrays.sort(allFiles, (f1, f2) -> f1.getName().compareTo(f2.getName()));
        
        // 使用全部文件作为测试数据
        int numTestFiles = allFiles.length;
        System.out.println("使用全部 " + numTestFiles + " 个CSV文件作为测试数据");

        // Get the solution
        ArrayList<Integer>[] zones = algo.getZones();
        ArrayList<Area> centers = algo.getCenters();
        boolean useRelativeBalance = algo.getUseRelativeBalance();

        int p = zones.length; // 区域数量
        double coeffLower = (1.0 - r) / p; // (1-r)/p
        double coeffUpper = (1.0 + r) / p; // (1+r)/p

        // Count satisfied scenarios
        int satisfiedScenarios = 0;
        int totalTestScenarios = 0;

        // For each test CSV file (使用全部文件)
        for (int fileIdx = 0; fileIdx < allFiles.length; fileIdx++) {
            File testFile = allFiles[fileIdx];
            
            try {
                // Load test scenario data from CSV
                double[][] testScenarioData = loadAssignmentDependentDataFromCSV(testFile, instance.getN());
                if (testScenarioData == null || testScenarioData.length == 0) {
                    System.err.println("警告: 无法加载测试文件 " + testFile.getName() + " 或数据为空");
                    continue;
                }
                
                totalTestScenarios++;
                boolean scenarioSatisfied = true;

                if (useRelativeBalance) {
                    // 相对平衡性约束：对于assignment-dependent模型，工作量是 d_ij
                    // 计算该场景的总工作量：sum_k sum_{i ∈ zone_k} d_ik
                    // 注意：在assignment-dependent模型中，区域索引j对应CSV文件中的Region_j列
                    double totalWorkload = 0.0;
                    for (int k = 0; k < zones.length; k++) {
                        if (zones[k] == null || zones[k].isEmpty()) {
                            continue;
                        }
                        // 区域索引k对应CSV文件中的Region_k列
                        for (int areaId : zones[k]) {
                            // 对于assignment-dependent模型，d_ik 直接从场景数据中获取
                            // testScenarioData[areaId][k] 就是 d_{areaId, k}（区域k）
                            if (areaId >= 0 && areaId < testScenarioData.length && 
                                testScenarioData[areaId] != null && 
                                k < testScenarioData[areaId].length) {
                                totalWorkload += testScenarioData[areaId][k];
                            }
                        }
                    }

                    // Check each district
                    for (int j = 0; j < zones.length; j++) {
                        if (zones[j] == null || zones[j].isEmpty()) {
                            continue;
                        }

                        // 计算区域 j 的工作量：∑_{i ∈ zone_j} d_ij
                        // 区域索引j对应CSV文件中的Region_j列
                        double districtWorkload = 0.0;
                        for (int areaId : zones[j]) {
                            // testScenarioData[areaId][j] 就是 d_{areaId, j}（区域j）
                            if (areaId >= 0 && areaId < testScenarioData.length && 
                                testScenarioData[areaId] != null && 
                                j < testScenarioData[areaId].length) {
                                districtWorkload += testScenarioData[areaId][j];
                            }
                        }

                        // 相对平衡性约束：(1-r)/p * totalWorkload <= districtWorkload <= (1+r)/p * totalWorkload
                        double lowerBound = coeffLower * totalWorkload;
                        double upperBound = coeffUpper * totalWorkload;
                        if (districtWorkload < lowerBound || districtWorkload > upperBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                } else {
                    // 绝对平衡性约束：使用固定的上下界
                    double demandLowerBound = algo.getDemandLowerBound();
                    double demandUpperBound = algo.getDemandUpperBound();

                    // Check each district
                    for (int j = 0; j < zones.length; j++) {
                        if (zones[j] == null || zones[j].isEmpty()) {
                            continue;
                        }

                        // 计算区域 j 的工作量：∑_{i ∈ zone_j} d_ij
                        // 区域索引j对应CSV文件中的Region_j列
                        double districtWorkload = 0.0;
                        for (int areaId : zones[j]) {
                            if (areaId >= 0 && areaId < testScenarioData.length && 
                                testScenarioData[areaId] != null && 
                                j < testScenarioData[areaId].length) {
                                districtWorkload += testScenarioData[areaId][j];
                            }
                        }

                        // 使用主算法中调整后的最终上下界检查约束
                        if (districtWorkload > demandUpperBound || districtWorkload < demandLowerBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                }

                if (scenarioSatisfied) {
                    satisfiedScenarios++;
                }
            } catch (Exception e) {
                System.err.println("错误: 处理测试文件 " + testFile.getName() + " 时出错: " + e.getMessage());
                e.printStackTrace();
            }
        }

        if (totalTestScenarios == 0) {
            System.err.println("错误: 没有有效的测试场景");
            return -1.0;
        }

        // Calculate and return percentage of satisfied scenarios
        double outOfSamplePerformance = (double) satisfiedScenarios / totalTestScenarios;

        // Output to console
        System.out.println(String.format("样本外性能: %.4f (%d/%d 场景满足约束)",
                outOfSamplePerformance, satisfiedScenarios, totalTestScenarios));

        return outOfSamplePerformance;
    }

    @SuppressWarnings("unchecked")
    private static AssignmentDependentDiagnostics collectAssignmentDependentDiagnostics(
            Instance instance,
            DistributionallyRobustAlgo algo,
            double r) {
        String dataDir = "data/travel_dist_dual_values_filtered_by_date_cluster20_unit_filled";
        File dir = new File(dataDir);
        File[] allFiles = dir.listFiles((d, name) -> name.endsWith(".csv") && name.startsWith("travel_dist_dual_values_p3"));

        if (allFiles == null || allFiles.length == 0) {
            throw new RuntimeException("在目录 " + dataDir + " 中未找到CSV文件");
        }

        Arrays.sort(allFiles, (f1, f2) -> f1.getName().compareTo(f2.getName()));

        ArrayList<Integer>[] zones = algo.getZones();
        int p = zones.length;
        double coeffLower = (1.0 - r) / p;
        double coeffUpper = (1.0 + r) / p;

        ArrayList<Double>[] districtShares = new ArrayList[p];
        ArrayList<Double>[] districtMargins = new ArrayList[p];
        int[] lowerViolations = new int[p];
        int[] upperViolations = new int[p];
        for (int j = 0; j < p; j++) {
            districtShares[j] = new ArrayList<>();
            districtMargins[j] = new ArrayList<>();
        }

        int satisfiedScenarios = 0;
        int totalScenarios = 0;

        for (File testFile : allFiles) {
            double[][] testScenarioData = loadAssignmentDependentDataFromCSV(testFile, instance.getN());
            if (testScenarioData == null || testScenarioData.length == 0) {
                continue;
            }

            totalScenarios++;
            boolean scenarioSatisfied = true;

            double totalWorkload = 0.0;
            double[] districtWorkloads = new double[p];

            for (int j = 0; j < p; j++) {
                if (zones[j] == null || zones[j].isEmpty()) {
                    continue;
                }

                double districtWorkload = 0.0;
                for (int areaId : zones[j]) {
                    if (areaId >= 0 && areaId < testScenarioData.length
                            && testScenarioData[areaId] != null
                            && j < testScenarioData[areaId].length) {
                        districtWorkload += testScenarioData[areaId][j];
                    }
                }
                districtWorkloads[j] = districtWorkload;
                totalWorkload += districtWorkload;
            }

            if (totalWorkload <= 1e-12) {
                for (int j = 0; j < p; j++) {
                    districtShares[j].add(0.0);
                    districtMargins[j].add(-coeffLower);
                    lowerViolations[j]++;
                }
                continue;
            }

            for (int j = 0; j < p; j++) {
                double share = districtWorkloads[j] / totalWorkload;
                double margin = Math.min(share - coeffLower, coeffUpper - share);
                districtShares[j].add(share);
                districtMargins[j].add(margin);

                if (share < coeffLower) {
                    lowerViolations[j]++;
                    scenarioSatisfied = false;
                } else if (share > coeffUpper) {
                    upperViolations[j]++;
                    scenarioSatisfied = false;
                }
            }

            if (scenarioSatisfied) {
                satisfiedScenarios++;
            }
        }

        AssignmentDependentDiagnostics diagnostics = new AssignmentDependentDiagnostics();
        diagnostics.totalScenarios = totalScenarios;
        diagnostics.satisfiedScenarios = satisfiedScenarios;
        diagnostics.lowerBoundShare = coeffLower;
        diagnostics.upperBoundShare = coeffUpper;
        diagnostics.districtStats = new DistrictDiagnostics[p];

        for (int j = 0; j < p; j++) {
            diagnostics.districtStats[j] = buildDistrictDiagnostics(
                    districtShares[j], districtMargins[j], lowerViolations[j], upperViolations[j]);
        }

        return diagnostics;
    }

    private static DistrictDiagnostics buildDistrictDiagnostics(
            ArrayList<Double> shares,
            ArrayList<Double> margins,
            int lowerViolations,
            int upperViolations) {
        DistrictDiagnostics stats = new DistrictDiagnostics();
        stats.lowerViolations = lowerViolations;
        stats.upperViolations = upperViolations;

        if (shares.isEmpty()) {
            stats.minShare = Double.NaN;
            stats.q05Share = Double.NaN;
            stats.medianShare = Double.NaN;
            stats.q95Share = Double.NaN;
            stats.maxShare = Double.NaN;
            stats.minMargin = Double.NaN;
            stats.medianMargin = Double.NaN;
            return stats;
        }

        double[] shareArr = toSortedArray(shares);
        double[] marginArr = toSortedArray(margins);

        stats.minShare = shareArr[0];
        stats.q05Share = quantileFromSorted(shareArr, 0.05);
        stats.medianShare = quantileFromSorted(shareArr, 0.50);
        stats.q95Share = quantileFromSorted(shareArr, 0.95);
        stats.maxShare = shareArr[shareArr.length - 1];
        stats.minMargin = marginArr[0];
        stats.medianMargin = quantileFromSorted(marginArr, 0.50);
        return stats;
    }

    private static double[] toSortedArray(ArrayList<Double> values) {
        double[] arr = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            arr[i] = values.get(i);
        }
        Arrays.sort(arr);
        return arr;
    }

    private static double quantileFromSorted(double[] sortedValues, double q) {
        if (sortedValues.length == 0) {
            return Double.NaN;
        }
        if (sortedValues.length == 1) {
            return sortedValues[0];
        }
        int idx = (int) Math.floor(q * (sortedValues.length - 1));
        if (idx < 0) {
            idx = 0;
        }
        if (idx >= sortedValues.length) {
            idx = sortedValues.length - 1;
        }
        return sortedValues[idx];
    }

    private static void writeAssignmentDependentDiagnostics(
            BufferedWriter diagnosticsWriter,
            AssignmentDependentDiagnostics diagnostics,
            String instanceName,
            double rsd,
            double r,
            double gamma,
            boolean useD1,
            double delta1,
            double delta2,
            boolean useJointChance,
            boolean useExactMethod,
            double objectiveValue,
            double outOfSamplePerformance,
            DistributionallyRobustAlgo algo) throws IOException {
        ArrayList<Integer>[] zones = algo.getZones();
        for (int j = 0; j < diagnostics.districtStats.length; j++) {
            DistrictDiagnostics stats = diagnostics.districtStats[j];
            diagnosticsWriter.write(String.format(
                    "%s,%.3f,%.1f,%.1f,%s,%.2f,%.2f,%s,%s,%.4f,%.4f,%d,%d,%.4f,%.4f,%d,\"%s\",%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%.4f,%.4f",
                    instanceName,
                    rsd,
                    r,
                    gamma,
                    useD1 ? "true" : "false",
                    delta1,
                    delta2,
                    useJointChance ? "true" : "false",
                    useExactMethod ? "true" : "false",
                    objectiveValue,
                    outOfSamplePerformance,
                    diagnostics.totalScenarios,
                    diagnostics.satisfiedScenarios,
                    diagnostics.lowerBoundShare,
                    diagnostics.upperBoundShare,
                    j,
                    formatZoneMembers(zones[j]),
                    stats.minShare,
                    stats.q05Share,
                    stats.medianShare,
                    stats.q95Share,
                    stats.maxShare,
                    stats.lowerViolations,
                    stats.upperViolations,
                    stats.minMargin,
                    stats.medianMargin
            ));
            diagnosticsWriter.newLine();
        }
        diagnosticsWriter.flush();
    }

    private static String formatZoneMembers(ArrayList<Integer> zone) {
        if (zone == null || zone.isEmpty()) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < zone.size(); i++) {
            if (i > 0) {
                sb.append(' ');
            }
            sb.append(zone.get(i));
        }
        return sb.toString();
    }

    /**
     * Load assignment-dependent data from a single CSV file
     * @param csvFile The CSV file to load
     * @param n Number of basic units
     * @return A 2D array [i][j] representing d_ij values, or null if error
     */
    private static double[][] loadAssignmentDependentDataFromCSV(File csvFile, int n) {
        // 检查n是否有效
        if (n <= 0) {
            System.err.println("警告: n = " + n + " 无效，无法加载数据");
            return null;
        }
        
        double[][] data = new double[n][n];
        
        // Initialize all values to 0
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                data[i][j] = 0.0;
            }
        }
        
        try {
            java.io.BufferedReader reader = new java.io.BufferedReader(
                new java.io.InputStreamReader(new java.io.FileInputStream(csvFile), "UTF-8"));
            
            String line = reader.readLine(); // 读取表头
            if (line == null) {
                System.err.println("警告: 文件 " + csvFile.getName() + " 为空");
                reader.close();
                return null;
            }
            
            // 解析表头，获取区域数量
            String[] headers = line.split(",");
            int numRegions = headers.length - 1; // 减去PointID列
            
            // 读取数据
            while ((line = reader.readLine()) != null) {
                String[] values = line.split(",");
                if (values.length < 2) continue;
                
                try {
                    int pointId = Integer.parseInt(values[0].trim());
                    if (pointId < 0 || pointId >= n) continue;
                    
                    // 读取每个区域的值
                    for (int j = 0; j < Math.min(numRegions, n); j++) {
                        String valueStr = values[j + 1].trim();
                        if (!valueStr.equals("Null") && !valueStr.isEmpty()) {
                            try {
                                double value = Double.parseDouble(valueStr);
                                data[pointId][j] = value;
                            } catch (NumberFormatException e) {
                                // 忽略无效数值
                            }
                        }
                    }
                } catch (NumberFormatException e) {
                    // 忽略无效行
                }
            }
            
            reader.close();
        } catch (Exception e) {
            System.err.println("错误: 读取文件 " + csvFile.getName() + " 时出错: " + e.getMessage());
            e.printStackTrace();
            return null;
        }
        
        return data;
    }
}