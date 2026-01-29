import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class CCTest2 {
    public static void main(String[] args) throws Exception {
        // Fixed parameters for the test
        String instanceFile = "./Instances/2DU120-05-11.dat";  // Specify the instance file
        double E = 50.0;                                    // Expected value
        double RSD = 0.5;                                  // Relative standard deviation
        double r = 0.2;                                     // Tolerance parameter
        double gamma = 0.2;                                 // Chance constraint risk parameter
        int numScenarios = 500;                            // Number of scenarios
        boolean useScenarioGeneration = true;               // Whether to use scenario generation
        long seed = 12345678;                               // Random seed

        long testSeed = seed + 1000;
        int numTestScenarios = 1000;

        // Output CSV file
        String outputCSVPath = "./output/chance_constrained_single_test.csv";

        try (BufferedWriter csvWriter = new BufferedWriter(new FileWriter(outputCSVPath))) {
            // Write CSV header
            csvWriter.write("Instance,RSD,r,gamma,Scenarios,UseScenarioGeneration,Runtime(s),Objective,OutOfSamplePerformance");
            csvWriter.newLine();

            System.out.println("Running experiment:");
            System.out.println("Instance: " + instanceFile);
            System.out.println("RSD: " + RSD);
            System.out.println("r: " + r);
            System.out.println("gamma: " + gamma);
            System.out.println("Scenarios: " + numScenarios);
            System.out.println("Use Scenario Generation: " + useScenarioGeneration);

            // Load the instance
            Instance instance = new Instance(instanceFile);

            // Generate random scenarios
            double[][] scenarios = generateScenarios(
                    instance.getN(), numScenarios, E, RSD, seed
            );

            // Record start time
            long startTime = System.currentTimeMillis();

            // Create algorithm instance
            ChanceConstrainedAlgo algo = new ChanceConstrainedAlgo(
                    instance, scenarios, gamma, seed, r
            );

            // Run algorithm
            double objectiveValue = 0;
            try {
                objectiveValue = algo.run("", useScenarioGeneration);
            } catch (Exception e) {
                System.err.println("Error running experiment: " + e.getMessage());
                e.printStackTrace();
                return;
            }

            // Calculate runtime
            long endTime = System.currentTimeMillis();
            double runtime = (endTime - startTime) / 1000.0;

            // Test out-of-sample performance
            double outOfSamplePerformance = testOutOfSamplePerformance(
                    instance, algo, E, RSD, testSeed, r, numTestScenarios);

            // Write to CSV file
            csvWriter.write(String.format(
                    "%s,%.3f,%.1f,%.1f,%d,%s,%.3f,%.4f,%.4f",
                    instanceFile, RSD, r, gamma, numScenarios,
                    useScenarioGeneration ? "true" : "false",
                    runtime, objectiveValue, outOfSamplePerformance
            ));
            csvWriter.newLine();

            // Print results to console
            System.out.println("Experiment completed");
            System.out.println("Runtime: " + runtime + " seconds");
            System.out.println("Objective value: " + objectiveValue);
            System.out.println("Out-of-sample performance: " + outOfSamplePerformance);

        } catch (IOException e) {
            System.err.println("Error writing results to CSV: " + e.getMessage());
        }

        System.out.println("Experiment completed, results saved to " + outputCSVPath);
    }

    // Generate random scenarios - using uniform distribution
    private static double[][] generateScenarios(int n, int numScenarios, double E, double RSD, long seed) {
        double[][] scenarios = new double[numScenarios][n];
        Random rand = new Random(seed); // Fixed seed for reproducibility

        // Calculate uniform distribution endpoints
        double lowerBound = E * (1 - Math.sqrt(3) * RSD);
        double upperBound = E * (1 + Math.sqrt(3) * RSD);

        // Generate scenarios
        for (int s = 0; s < numScenarios; s++) {
            for (int i = 0; i < n; i++) {
                // Generate demand from uniform distribution
                double demand = lowerBound + rand.nextDouble() * (upperBound - lowerBound);

                // Ensure demand is positive
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