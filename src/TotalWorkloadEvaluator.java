import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;

public class TotalWorkloadEvaluator {
    private static final List<String> DEFAULT_INPUT_DIRS = Arrays.asList(
            "./output/csv",
            "./output/csv/chance_constrained_results"
    );
    private static final String DEFAULT_OUTPUT_DIR = "./output/csv/total workload";
    private static final String DEFAULT_INSTANCE_DIR = "./Instances_new";

    private static final double E = 50.0;
    private static final long TRAIN_SEED = 12345678L;
    private static final long TEST_SEED = TRAIN_SEED + 1000L;
    private static final int NUM_TEST_SCENARIOS = 1000;

    private static final List<String> DROPPED_OUTPUT_COLUMNS = Arrays.asList(
            "UseJointChance",
            "UseExactMethod",
            "UseImprovedModel",
            "UseAssignmentDependent",
            "AvgDistMethod",
            "AvgDistMean",
            "Runtime(s)",
            "StatusCode",
            "CutIterations",
            "FailureStage",
            "Status",
            "OptimalAssignment",
            "FinalCenters"
    );

    private static final DecimalFormat RSD_FORMAT = new DecimalFormat("0.###");
    private static final DecimalFormat WORKLOAD_FORMAT = new DecimalFormat("0.0000");

    private static class CsvTable {
        List<String> headers = new ArrayList<>();
        List<LinkedHashMap<String, String>> rows = new ArrayList<>();
    }

    private static class WorkloadResult {
        String testRsdText = "-1";
        String meanTotalWorkloadText = "-1";
        String scenarioCountText = "0";
        String status = "unknown";
    }

    public static void main(String[] args) throws Exception {
        List<File> inputFiles = resolveInputFiles(args);
        if (inputFiles.isEmpty()) {
            System.out.println("No matching CSV files found to evaluate.");
            return;
        }

        File outputDir = new File(DEFAULT_OUTPUT_DIR);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Failed to create output directory: " + outputDir.getPath());
        }

        Map<String, Instance> instanceCache = new LinkedHashMap<>();
        for (File inputFile : inputFiles) {
            processFile(inputFile, outputDir, instanceCache);
        }

        System.out.println("Total workload evaluation completed.");
    }

    private static List<File> resolveInputFiles(String[] args) throws IOException {
        LinkedHashSet<String> seenPaths = new LinkedHashSet<>();
        List<File> files = new ArrayList<>();

        if (args != null && args.length > 0) {
            for (String arg : args) {
                addInputPath(new File(arg), files, seenPaths);
            }
        } else {
            for (String dirPath : DEFAULT_INPUT_DIRS) {
                addInputPath(new File(dirPath), files, seenPaths);
            }
        }
        return files;
    }

    private static void addInputPath(File path, List<File> files, LinkedHashSet<String> seenPaths)
            throws IOException {
        if (path.isDirectory()) {
            File[] matched = path.listFiles((dir, name) -> isTargetResultCsv(name));
            if (matched != null) {
                Arrays.sort(matched, (a, b) -> a.getName().compareToIgnoreCase(b.getName()));
                for (File file : matched) {
                    addFileIfNew(file, files, seenPaths);
                }
            }
        } else if (path.isFile()) {
            addFileIfNew(path, files, seenPaths);
        }
    }

    private static void addFileIfNew(File file, List<File> files, LinkedHashSet<String> seenPaths)
            throws IOException {
        String canonicalPath = file.getCanonicalPath();
        if (seenPaths.add(canonicalPath)) {
            files.add(file);
        }
    }

    private static boolean isTargetResultCsv(String name) {
        String lower = name.toLowerCase(Locale.ROOT);
        return lower.endsWith("1000test.csv")
                && !lower.startsWith("cc_vs_")
                && !lower.endsWith("_harder_oos.csv")
                && !lower.endsWith("_total_workload.csv");
    }

    private static void processFile(File inputFile, File outputDir, Map<String, Instance> instanceCache)
            throws Exception {
        System.out.println("Processing file: " + inputFile.getPath());
        CsvTable table = readCsv(inputFile);

        if (!table.headers.contains("OptimalAssignment") || !table.headers.contains("FinalCenters")) {
            System.out.println("Skip file without required solution columns: " + inputFile.getName());
            return;
        }

        List<String> outputHeaders = buildOutputHeaders(table.headers);
        List<LinkedHashMap<String, String>> outputRows = new ArrayList<>();
        int successCount = 0;

        for (LinkedHashMap<String, String> row : table.rows) {
            LinkedHashMap<String, String> outputRow = new LinkedHashMap<>();
            for (String header : outputHeaders) {
                if (!isGeneratedOutputColumn(header)) {
                    outputRow.put(header, row.getOrDefault(header, ""));
                }
            }

            WorkloadResult result = evaluateRow(row, instanceCache);
            outputRow.put("WorkloadTestRSD", result.testRsdText);
            outputRow.put("MeanTotalWorkload", result.meanTotalWorkloadText);
            outputRow.put("WorkloadScenarioCount", result.scenarioCountText);
            outputRow.put("TotalWorkloadStatus", result.status);

            if ("success".equals(result.status)) {
                successCount++;
            }
            outputRows.add(outputRow);
        }

        File outputFile = new File(outputDir, inputFile.getName().replace(".csv", "_total_workload.csv"));
        writeCsv(outputFile, outputHeaders, outputRows);
        System.out.println("Wrote " + outputRows.size() + " rows, successful evaluations: " + successCount
                + ", output: " + outputFile.getPath());
    }

    private static List<String> buildOutputHeaders(List<String> inputHeaders) {
        List<String> outputHeaders = new ArrayList<>();
        for (String header : inputHeaders) {
            if (!DROPPED_OUTPUT_COLUMNS.contains(header)) {
                outputHeaders.add(header);
            }
        }
        appendIfMissing(outputHeaders, "WorkloadTestRSD");
        appendIfMissing(outputHeaders, "MeanTotalWorkload");
        appendIfMissing(outputHeaders, "WorkloadScenarioCount");
        appendIfMissing(outputHeaders, "TotalWorkloadStatus");
        return outputHeaders;
    }

    private static boolean isGeneratedOutputColumn(String header) {
        return "WorkloadTestRSD".equals(header)
                || "MeanTotalWorkload".equals(header)
                || "WorkloadScenarioCount".equals(header)
                || "TotalWorkloadStatus".equals(header);
    }

    private static WorkloadResult evaluateRow(Map<String, String> row, Map<String, Instance> instanceCache) {
        WorkloadResult result = new WorkloadResult();
        try {
            if (parseBoolean(row.get("UseAssignmentDependent"))) {
                result.status = "assignment_dependent_unsupported";
                return result;
            }

            double rsd = parseDouble(row.get("RSD"), Double.NaN);
            if (Double.isNaN(rsd)) {
                result.status = "invalid_rsd";
                return result;
            }
            result.testRsdText = RSD_FORMAT.format(rsd);

            String optimalAssignment = safeTrim(row.get("OptimalAssignment"));
            String finalCenters = safeTrim(row.get("FinalCenters"));
            if (isMissingSolution(optimalAssignment) || isMissingSolution(finalCenters)) {
                result.status = "solve_failed";
                return result;
            }

            String instanceName = safeTrim(row.get("InstanceName"));
            Instance instance = loadInstance(instanceName, instanceCache);
            if (instance == null) {
                result.status = "instance_not_found";
                return result;
            }

            int numRegions = parseInt(row.get("NumRegions"), instance.k);
            ArrayList<Integer>[] zones = parseOptimalAssignment(optimalAssignment, numRegions);
            ArrayList<Integer> centers = parseFinalCenters(finalCenters);
            if (zones == null || centers == null || centers.size() != numRegions) {
                result.status = "invalid_solution_format";
                return result;
            }

            boolean useImprovedModel = parseBoolean(row.get("UseImprovedModel"));
            double[][] distanceMatrix = useImprovedModel
                    ? computeShortestPathDistances(instance)
                    : instance.dist;

            double[][] testScenarios = generateScenarios(
                    instance.getN(), NUM_TEST_SCENARIOS, E, rsd, TEST_SEED);
            double meanTotalWorkload = computeMeanTotalWorkload(
                    zones, centers, testScenarios, distanceMatrix);

            if (meanTotalWorkload < 0.0) {
                result.status = "invalid_distance_or_solution";
                return result;
            }

            result.meanTotalWorkloadText = WORKLOAD_FORMAT.format(meanTotalWorkload);
            result.scenarioCountText = Integer.toString(NUM_TEST_SCENARIOS);
            result.status = "success";
            return result;
        } catch (Exception e) {
            result.status = "error:" + e.getClass().getSimpleName();
            return result;
        }
    }

    private static double computeMeanTotalWorkload(
            ArrayList<Integer>[] zones,
            ArrayList<Integer> centers,
            double[][] testScenarios,
            double[][] distanceMatrix) {
        double totalAcrossScenarios = 0.0;

        for (double[] scenario : testScenarios) {
            double scenarioTotalWorkload = 0.0;

            for (int district = 0; district < zones.length; district++) {
                if (zones[district] == null || zones[district].isEmpty()) {
                    continue;
                }

                int centerId = centers.get(district);
                if (centerId < 0 || centerId >= distanceMatrix.length) {
                    return -1.0;
                }

                for (int areaId : zones[district]) {
                    if (areaId < 0 || areaId >= scenario.length || areaId >= distanceMatrix.length) {
                        return -1.0;
                    }
                    double distance = distanceMatrix[areaId][centerId];
                    if (distance == Double.MAX_VALUE || Double.isInfinite(distance) || Double.isNaN(distance)) {
                        return -1.0;
                    }
                    scenarioTotalWorkload += scenario[areaId] * distance;
                }
            }

            totalAcrossScenarios += scenarioTotalWorkload;
        }

        return totalAcrossScenarios / testScenarios.length;
    }

    private static boolean isMissingSolution(String text) {
        return text.isEmpty() || "-1".equals(text);
    }

    private static String safeTrim(String value) {
        return value == null ? "" : value.trim();
    }

    private static void appendIfMissing(List<String> headers, String header) {
        if (!headers.contains(header)) {
            headers.add(header);
        }
    }

    private static Instance loadInstance(String instanceName, Map<String, Instance> instanceCache)
            throws Exception {
        if (instanceCache.containsKey(instanceName)) {
            return instanceCache.get(instanceName);
        }

        File instanceFile = new File(DEFAULT_INSTANCE_DIR, instanceName);
        if (!instanceFile.exists()) {
            return null;
        }

        Instance instance = new Instance(instanceFile.getPath());
        instanceCache.put(instanceName, instance);
        return instance;
    }

    private static ArrayList<Integer>[] parseOptimalAssignment(String text, int numRegions) {
        @SuppressWarnings("unchecked")
        ArrayList<Integer>[] zones = new ArrayList[numRegions];
        for (int i = 0; i < numRegions; i++) {
            zones[i] = new ArrayList<>();
        }

        String[] parts = text.split("\\|");
        for (String part : parts) {
            String trimmed = part.trim();
            if (trimmed.isEmpty()) {
                continue;
            }
            int colonIndex = trimmed.indexOf(':');
            if (colonIndex < 0) {
                return null;
            }

            int zoneIndex = Integer.parseInt(trimmed.substring(0, colonIndex).trim());
            if (zoneIndex < 0 || zoneIndex >= numRegions) {
                return null;
            }

            String membersText = trimmed.substring(colonIndex + 1).trim();
            if (membersText.isEmpty()) {
                continue;
            }

            String[] members = membersText.split("\\s+");
            for (String member : members) {
                if (!member.isEmpty()) {
                    zones[zoneIndex].add(Integer.parseInt(member));
                }
            }
        }
        return zones;
    }

    private static ArrayList<Integer> parseFinalCenters(String text) {
        ArrayList<Integer> centers = new ArrayList<>();
        String[] parts = text.trim().split("\\s+");
        for (String part : parts) {
            if (!part.isEmpty()) {
                centers.add(Integer.parseInt(part));
            }
        }
        return centers;
    }

    private static double[][] generateScenarios(int n, int numScenarios, double expectedValue, double rsd, long seed) {
        double[][] scenarios = new double[numScenarios][n];
        Random rand = new Random(seed);

        double lowerBound = expectedValue * (1 - Math.sqrt(3) * rsd);
        double upperBound = expectedValue * (1 + Math.sqrt(3) * rsd);

        for (int s = 0; s < numScenarios; s++) {
            for (int i = 0; i < n; i++) {
                double demand = lowerBound + rand.nextDouble() * (upperBound - lowerBound);
                scenarios[s][i] = Math.max(1, demand);
            }
        }
        return scenarios;
    }

    private static double[][] computeShortestPathDistances(Instance instance) {
        int n = instance.getN();
        double[][] shortestPathDist = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    shortestPathDist[i][j] = 0.0;
                } else if (instance.getEdges()[i][j] == 1) {
                    shortestPathDist[i][j] = instance.dist[i][j];
                } else {
                    shortestPathDist[i][j] = Double.MAX_VALUE;
                }
            }
        }

        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                if (shortestPathDist[i][k] == Double.MAX_VALUE) {
                    continue;
                }
                for (int j = 0; j < n; j++) {
                    if (shortestPathDist[k][j] == Double.MAX_VALUE) {
                        continue;
                    }
                    double newDist = shortestPathDist[i][k] + shortestPathDist[k][j];
                    if (newDist < shortestPathDist[i][j]) {
                        shortestPathDist[i][j] = newDist;
                    }
                }
            }
        }
        return shortestPathDist;
    }

    private static CsvTable readCsv(File file) throws IOException {
        CsvTable table = new CsvTable();
        try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8))) {
            String line = reader.readLine();
            if (line == null) {
                return table;
            }
            table.headers = parseCsvLine(stripBom(line));

            while ((line = reader.readLine()) != null) {
                List<String> values = parseCsvLine(line);
                LinkedHashMap<String, String> row = new LinkedHashMap<>();
                for (int i = 0; i < table.headers.size(); i++) {
                    String value = i < values.size() ? values.get(i) : "";
                    row.put(table.headers.get(i), value);
                }
                table.rows.add(row);
            }
        }
        return table;
    }

    private static void writeCsv(File file, List<String> headers, List<LinkedHashMap<String, String>> rows)
            throws IOException {
        try (BufferedWriter writer = new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(file), StandardCharsets.UTF_8))) {
            writer.write('\uFEFF');
            writer.write(toCsvLine(headers));
            writer.newLine();
            for (LinkedHashMap<String, String> row : rows) {
                List<String> values = new ArrayList<>();
                for (String header : headers) {
                    values.add(row.getOrDefault(header, ""));
                }
                writer.write(toCsvLine(values));
                writer.newLine();
            }
        }
    }

    private static String stripBom(String text) {
        if (!text.isEmpty() && text.charAt(0) == '\uFEFF') {
            return text.substring(1);
        }
        return text;
    }

    private static List<String> parseCsvLine(String line) {
        List<String> values = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean inQuotes = false;

        for (int i = 0; i < line.length(); i++) {
            char ch = line.charAt(i);
            if (ch == '"') {
                if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                    current.append('"');
                    i++;
                } else {
                    inQuotes = !inQuotes;
                }
            } else if (ch == ',' && !inQuotes) {
                values.add(current.toString());
                current.setLength(0);
            } else {
                current.append(ch);
            }
        }
        values.add(current.toString());
        return values;
    }

    private static String toCsvLine(List<String> values) {
        List<String> escaped = new ArrayList<>(values.size());
        for (String value : values) {
            String safe = value == null ? "" : value;
            if (safe.contains(",") || safe.contains("\"") || safe.contains("\n") || safe.contains("\r")) {
                escaped.add("\"" + safe.replace("\"", "\"\"") + "\"");
            } else {
                escaped.add(safe);
            }
        }
        return String.join(",", escaped);
    }

    private static double parseDouble(String text, double defaultValue) {
        try {
            return Double.parseDouble(safeTrim(text));
        } catch (Exception e) {
            return defaultValue;
        }
    }

    private static int parseInt(String text, int defaultValue) {
        try {
            String trimmed = safeTrim(text);
            if (trimmed.contains(".")) {
                return (int) Math.round(Double.parseDouble(trimmed));
            }
            return Integer.parseInt(trimmed);
        } catch (Exception e) {
            return defaultValue;
        }
    }

    private static boolean parseBoolean(String text) {
        return "true".equalsIgnoreCase(safeTrim(text));
    }
}
