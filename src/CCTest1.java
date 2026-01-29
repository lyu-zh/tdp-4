import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;

public class CCTest1 {
    public static void main(String[] args) throws Exception {
        // Experiment configuration
        double E = 50.0; // Expected value
        double[] RSDValues = {0.125, 0.25}; // Relative standard deviation array
        double[] rValues = {0.5, 0.4, 0.3}; // Tolerance parameter values
        double[] gammaValues = {0.4, 0.3, 0.2, 0.1}; // Chance constraint risk parameter
        int[] scenarioNumValues = {1000}; // Number of scenarios
        boolean[] useScenarioGeneration = {true}; // Whether to use scenario generation
        long seed = 12345678; // Random seed

        long testSeed = seed + 1000;
        int numTestScenarios = 5000;

        // Get start time for filename
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmmss");
        String timestamp = dateFormat.format(new Date());
        
        // Output CSV file - save to output/csv folder with timestamp
        String outputCSVPath = "./output/csv/chance_constrained_results/chance_constrained_results_" + timestamp + ".csv";
        
        System.out.println("CSV文件将保存到: " + outputCSVPath);
        System.out.println("实验开始时间: " + new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));

        // Prepare CSV file
        // Modified code for CCTest1 class - CSV output section

// Prepare CSV file
        try (BufferedWriter csvWriter = new BufferedWriter(new FileWriter(outputCSVPath))) {
            // Write CSV header with modified columns - replacing Instance with NumUnits and NumRegions
            csvWriter.write("InstanceName,NumUnits,NumRegions,RSD,r,gamma,Scenarios,UseScenarioGeneration,Runtime(s),Objective,OutOfSamplePerformance,Status");
            csvWriter.newLine();
            csvWriter.flush(); // Flush header immediately to ensure it's saved even if program terminates early

            // Get all .dat files in Instances_new directory
            File dir = new File("./Instances_new");
            File[] instanceFiles = dir.listFiles((d, name) -> name.endsWith(".dat"));

            // If no instance files are found, provide a hint
            if (instanceFiles == null || instanceFiles.length == 0) {
                System.out.println("No .dat files found in ./Instances_new directory.");
                return;
            }
            int cnt = 0;
            // 定义需要处理的文件名前缀列表
            String[] targetPrefixes = {"GG20", "GG50", "GG80"};
            System.out.println("只处理文件名以以下前缀开头的实例文件:");
            for (String prefix : targetPrefixes) {
                System.out.println("  - " + prefix);
            }
            
            // Iterate through all instance files
            for (File instanceFile : instanceFiles) {
                String instanceName = instanceFile.getName();
                String fileName = instanceName; // 文件名（不含路径）

                // Filter: only process instances with matching prefix
                boolean shouldProcess = false;
                for (String prefix : targetPrefixes) {
                    if (fileName.startsWith(prefix)) {
                        shouldProcess = true;
                        break;
                    }
                }
                
                if (!shouldProcess) {
                    System.out.println("跳过实例 " + instanceName + " (不匹配任何目标前缀)");
                    continue;
                }
                
                System.out.println("处理实例 " + instanceName);

                // Load instance
                Instance instance = new Instance(instanceFile.getPath());
                int numUnits = instance.getN();

                // Iterate through all parameter combinations
                for (double RSD : RSDValues) {
                    for (double r : rValues) {
                        for (double gamma : gammaValues) {
                            for (int numScenarios : scenarioNumValues) {
                                for (boolean useScenario : useScenarioGeneration) {
                                    // Print current experiment information
                                    System.out.println("current inst number:" + cnt + "    >>>>>>>>>>>Running experiment:");
                                    cnt++;
                                    System.out.println("Instance: " + instanceName);
                                    System.out.println("RSD: " + RSD);
                                    System.out.println("r: " + r);
                                    System.out.println("gamma: " + gamma);
                                    System.out.println("Scenarios: " + numScenarios);
                                    System.out.println("Use Scenario Generation: " + useScenario);

                                    // Get the number of basic units and regions from the instance
                                    int numRegions = instance.k;

                                    double[][] scenarios = generateScenarios(
                                            numUnits, numScenarios, E, RSD, seed
                                    );

                                    // Record start time
                                    long startTime = System.currentTimeMillis();

                                    // Create algorithm instance
                                    ChanceConstrainedAlgo algo = new ChanceConstrainedAlgo(
                                            instance, scenarios, gamma, seed, r
                                    );

                                    // Run algorithm and get objective value
                                    double objectiveValue = 0;
                                    try {
                                        // Pass the useScenario parameter to control scenario generation
                                        objectiveValue = algo.run("", useScenario);
                                    } catch (Exception e) {
                                        System.err.println("Error running experiment: " + e.getMessage());
                                        continue;
                                    }
                                    if (objectiveValue == -1) {
                                        System.out.println("Error: Objective value is -1, skipping this instance:" + instanceName);
                                    }
                                    // Calculate runtime
                                    long endTime = System.currentTimeMillis();
                                    double runtime = (endTime - startTime) / 1000.0;

                                    // Test out-of-sample performance
                                    double outOfSamplePerformance = testOutOfSamplePerformance(
                                            instance, algo, E, RSD, testSeed, r, numTestScenarios);

                                    // Get the last model status code
                                    int statusCode = algo.getLastModelStatus();

                                    // Write to CSV file with the modified fields (numUnits and numRegions instead of instanceName)
                                    csvWriter.write(String.format(
                                            "%s,%d,%d,%.3f,%.1f,%.1f,%d,%s,%.3f,%.4f,%.4f,%d",
                                            instanceName, numUnits, numRegions, RSD, r, gamma, numScenarios,
                                            useScenario ? "true" : "false",
                                            runtime, objectiveValue, outOfSamplePerformance, statusCode
                                    ));
                                    csvWriter.newLine();

                                    // Flush to ensure real-time writing
                                    csvWriter.flush();
                                }
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            System.err.println("Error writing results to CSV: " + e.getMessage());
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
     * @param testSeed         Random seed for test scenario generation
     * @param r                Capacity tolerance parameter
     * @param numTestScenarios Number of test scenarios to generate
     * @return The percentage of test scenarios where constraints are satisfied
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
        ArrayList<Area> centers = algo.getCenters();

        // Count satisfied scenarios
        int satisfiedScenarios = 0;

        // For each test scenario
        for (int s = 0; s < numTestScenarios; s++) {
            // Calculate this scenario's total demand
            double scenarioTotalDemand = 0;
            for (int i = 0; i < instance.getN(); i++) {
                scenarioTotalDemand += testScenarios[s][i];
            }

            // Calculate scenario-specific capacity limit
            double scenarioCapacityLimit = (1 + r) * (scenarioTotalDemand / instance.k);

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

                // Check if capacity constraint is violated
                if (districtDemand > scenarioCapacityLimit) {
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