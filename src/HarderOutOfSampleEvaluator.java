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
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;

public class HarderOutOfSampleEvaluator {
    private static final String DEFAULT_INPUT_DIR = "./output/csv/chance_constrained_results";
    private static final String DEFAULT_OUTPUT_DIR = "./output/csv/harder out of sample";
    private static final String DEFAULT_INSTANCE_DIR = "./Instances_new";

    private static final double E = 50.0;
    private static final long TRAIN_SEED = 12345678L;
    private static final long TEST_SEED = TRAIN_SEED + 1000L;
    private static final int NUM_TEST_SCENARIOS = 5000;
    private static final boolean USE_RELATIVE_BALANCE = true;
    private static final double EPS = 1e-9;
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
            "OptimalAssignment",
            "FinalCenters"
    );

    private static final DecimalFormat RSD_FORMAT = new DecimalFormat("0.###");
    private static final DecimalFormat OOS_FORMAT = new DecimalFormat("0.0000");

    private static class CsvTable {
        List<String> headers = new ArrayList<>();
        List<LinkedHashMap<String, String>> rows = new ArrayList<>();
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

        System.out.println("Harder out-of-sample evaluation completed.");
    }

    private static List<File> resolveInputFiles(String[] args) {
        List<File> files = new ArrayList<>();
        if (args != null && args.length > 0) {
            for (String arg : args) {
                File file = new File(arg);
                if (file.isDirectory()) {
                    File[] matched = file.listFiles((dir, name) -> isTargetResultCsv(name));
                    if (matched != null) {
                        Arrays.sort(matched, (a, b) -> a.getName().compareToIgnoreCase(b.getName()));
                        files.addAll(Arrays.asList(matched));
                    }
                } else if (file.isFile()) {
                    files.add(file);
                }
            }
        } else {
            File inputDir = new File(DEFAULT_INPUT_DIR);
            File[] matched = inputDir.listFiles((dir, name) -> isTargetResultCsv(name));
            if (matched != null) {
                Arrays.sort(matched, (a, b) -> a.getName().compareToIgnoreCase(b.getName()));
                files.addAll(Arrays.asList(matched));
            }
        }
        return files;
    }

    private static boolean isTargetResultCsv(String name) {
        String lower = name.toLowerCase(Locale.ROOT);
        return lower.endsWith("1000train.csv") && !lower.startsWith("cc_vs_");
    }

    private static void processFile(File inputFile, File outputDir, Map<String, Instance> instanceCache) throws Exception {
        System.out.println("Processing file: " + inputFile.getName());
        CsvTable table = readCsv(inputFile);

        if (!table.headers.contains("OptimalAssignment") || !table.headers.contains("FinalCenters")) {
            System.out.println("Skip file without required columns: " + inputFile.getName());
            return;
        }

        List<String> outputHeaders = new ArrayList<>();
        for (String header : table.headers) {
            if (!DROPPED_OUTPUT_COLUMNS.contains(header)) {
                outputHeaders.add(header);
            }
        }
        appendIfMissing(outputHeaders, "HarderTestRSD");

        List<LinkedHashMap<String, String>> outputRows = new ArrayList<>();
        int successCount = 0;

        for (LinkedHashMap<String, String> row : table.rows) {
            LinkedHashMap<String, String> outputRow = new LinkedHashMap<>();
            for (String header : outputHeaders) {
                if (!header.startsWith("Harder")) {
                    outputRow.put(header, row.getOrDefault(header, ""));
                }
            }

            EvaluationResult result = evaluateRow(row, instanceCache);
            outputRow.put("HarderTestRSD", result.harderTestRsdText);
            outputRow.put("OutOfSamplePerformance", result.harderOosText);

            if ("success".equals(result.status)) {
                successCount++;
            }
            outputRows.add(outputRow);
        }

        File outputFile = new File(outputDir, inputFile.getName().replace(".csv", "_harder_oos.csv"));
        writeCsv(outputFile, outputHeaders, outputRows);
        System.out.println("Wrote " + outputRows.size() + " rows, successful evaluations: " + successCount
                + ", output: " + outputFile.getPath());
    }

    private static class EvaluationResult {
        String harderTestRsdText = "-1";
        String harderOosText = "-1";
        String status = "unknown";
    }

    private static EvaluationResult evaluateRow(Map<String, String> row, Map<String, Instance> instanceCache) {
        EvaluationResult result = new EvaluationResult();
        try {
            double originalRsd = parseDouble(row.get("RSD"), Double.NaN);
            double harderRsd = mapToHarderRsd(originalRsd);
            if (Double.isNaN(harderRsd)) {
                result.status = "unsupported_rsd";
                return result;
            }
            result.harderTestRsdText = RSD_FORMAT.format(harderRsd);

            if (parseBoolean(row.get("UseAssignmentDependent"))) {
                result.status = "assignment_dependent_unsupported";
                return result;
            }

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
            double r = parseDouble(row.get("r"), Double.NaN);
            int numTrainingScenarios = parseInt(row.get("Scenarios"), 1000);
            int avgDistMethod = parseInt(row.get("AvgDistMethod"), 3);

            double[][] harderTestScenarios = generateScenarios(
                    instance.getN(), NUM_TEST_SCENARIOS, E, harderRsd, TEST_SEED);

            double[][] shortestPathDist = null;
            Bounds bounds = null;
            if (useImprovedModel || !USE_RELATIVE_BALANCE) {
                double[][] trainingScenarios = generateScenarios(
                        instance.getN(), numTrainingScenarios, E, originalRsd, TRAIN_SEED);
                double[] meanVector = computeMeanVector(trainingScenarios);

                if (useImprovedModel) {
                    shortestPathDist = computeShortestPathDistances(instance);
                    bounds = computeImprovedBounds(instance, shortestPathDist, meanVector, r, avgDistMethod);
                } else {
                    bounds = computeOriginalBounds(instance, meanVector, r);
                }
            }

            double harderOos = computeOutOfSamplePerformance(
                    instance, zones, centers, harderTestScenarios, r,
                    useImprovedModel, USE_RELATIVE_BALANCE, shortestPathDist, bounds);
            result.harderOosText = OOS_FORMAT.format(harderOos);
            result.status = "success";
            return result;
        } catch (Exception e) {
            result.status = "error:" + e.getClass().getSimpleName();
            return result;
        }
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

    private static Instance loadInstance(String instanceName, Map<String, Instance> instanceCache) throws Exception {
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

    private static double mapToHarderRsd(double rsd) {
        if (approxEquals(rsd, 0.125)) {
            return 0.25;
        }
        if (approxEquals(rsd, 0.25)) {
            return 0.5;
        }
        if (approxEquals(rsd, 0.5)) {
            return 0.75;
        }
        return Double.NaN;
    }

    private static boolean approxEquals(double a, double b) {
        return Math.abs(a - b) < EPS;
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

    private static double[] computeMeanVector(double[][] scenarios) {
        int numScenarios = scenarios.length;
        int n = scenarios[0].length;
        double[] meanVector = new double[n];

        for (double[] scenario : scenarios) {
            for (int i = 0; i < n; i++) {
                meanVector[i] += scenario[i];
            }
        }

        for (int i = 0; i < n; i++) {
            meanVector[i] /= numScenarios;
        }
        return meanVector;
    }

    private static class Bounds {
        double lowerBound;
        double upperBound;
    }

    private static Bounds computeOriginalBounds(Instance instance, double[] meanVector, double r) {
        Bounds bounds = new Bounds();
        double totalMeanDemand = 0.0;
        for (double value : meanVector) {
            totalMeanDemand += value;
        }
        double mu = totalMeanDemand / instance.k;
        bounds.upperBound = (1 + r) * mu;
        bounds.lowerBound = (1 - r) * mu;
        return bounds;
    }

    private static Bounds computeImprovedBounds(
            Instance instance, double[][] shortestPathDist, double[] meanVector, double r, int avgDistMethod) {
        int n = instance.getN();
        int p = instance.k;
        int k = (int) Math.ceil((double) n / p);
        k = Math.min(k, n);

        double[] avgDist = computeAverageDistances(shortestPathDist, n, k, avgDistMethod);
        double totalExpectedWorkload = 0.0;
        for (int i = 0; i < n; i++) {
            totalExpectedWorkload += meanVector[i] * avgDist[i];
        }

        double avgWorkload = totalExpectedWorkload / p;
        Bounds bounds = new Bounds();
        bounds.upperBound = (1 + r) * avgWorkload;
        bounds.lowerBound = (1 - r) * avgWorkload;
        return bounds;
    }

    private static double[] computeAverageDistances(double[][] shortestPathDist, int n, int k, int method) {
        double[] avgDist = new double[n];

        if (method == 1) {
            for (int i = 0; i < n; i++) {
                ArrayList<Integer> validJs = findValidJs(shortestPathDist, n, k, i);
                if (!validJs.isEmpty()) {
                    double sum = 0.0;
                    for (int j : validJs) {
                        sum += shortestPathDist[j][i];
                    }
                    avgDist[i] = sum / validJs.size();
                } else {
                    avgDist[i] = averageOfKNearest(shortestPathDist[i], k);
                }
            }
        } else if (method == 2) {
            for (int i = 0; i < n; i++) {
                ArrayList<Integer> validJs = findValidJs(shortestPathDist, n, k, i);
                if (!validJs.isEmpty()) {
                    double sumOfAverages = 0.0;
                    for (int j : validJs) {
                        sumOfAverages += averageOfKNearest(shortestPathDist[j], k);
                    }
                    avgDist[i] = sumOfAverages / validJs.size();
                } else {
                    avgDist[i] = averageOfKNearest(shortestPathDist[i], k);
                }
            }
        } else {
            for (int i = 0; i < n; i++) {
                avgDist[i] = averageOfKNearest(shortestPathDist[i], k);
            }
        }

        return avgDist;
    }

    private static ArrayList<Integer> findValidJs(double[][] shortestPathDist, int n, int k, int targetIndex) {
        ArrayList<Integer> validJs = new ArrayList<>();
        for (int j = 0; j < n; j++) {
            double[] distances = new double[n];
            int[] indices = new int[n];
            for (int idx = 0; idx < n; idx++) {
                distances[idx] = shortestPathDist[j][idx];
                indices[idx] = idx;
            }
            sortByDistanceWithIndex(distances, indices);

            boolean found = false;
            for (int idx = 0; idx < k; idx++) {
                if (indices[idx] == targetIndex) {
                    found = true;
                    break;
                }
            }
            if (found) {
                validJs.add(j);
            }
        }
        return validJs;
    }

    private static void sortByDistanceWithIndex(double[] distances, int[] indices) {
        for (int a = 0; a < distances.length - 1; a++) {
            for (int b = a + 1; b < distances.length; b++) {
                if (distances[a] > distances[b]) {
                    double tempDist = distances[a];
                    distances[a] = distances[b];
                    distances[b] = tempDist;

                    int tempIdx = indices[a];
                    indices[a] = indices[b];
                    indices[b] = tempIdx;
                }
            }
        }
    }

    private static double averageOfKNearest(double[] rowDistances, int k) {
        double[] copy = Arrays.copyOf(rowDistances, rowDistances.length);
        Arrays.sort(copy);
        double sum = 0.0;
        for (int idx = 0; idx < k; idx++) {
            sum += copy[idx];
        }
        return sum / k;
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

    private static double computeOutOfSamplePerformance(
            Instance instance,
            ArrayList<Integer>[] zones,
            ArrayList<Integer> centers,
            double[][] testScenarios,
            double r,
            boolean useImprovedModel,
            boolean useRelativeBalance,
            double[][] shortestPathDist,
            Bounds bounds) {

        int p = zones.length;
        double coeffLower = (1.0 - r) / p;
        double coeffUpper = (1.0 + r) / p;
        int satisfiedScenarios = 0;

        for (int s = 0; s < testScenarios.length; s++) {
            boolean scenarioSatisfied = true;

            if (useRelativeBalance) {
                if (useImprovedModel && shortestPathDist != null) {
                    double totalWorkload = 0.0;
                    for (int k = 0; k < zones.length; k++) {
                        if (zones[k] == null || zones[k].isEmpty()) {
                            continue;
                        }
                        int centerId = centers.get(k);
                        for (int areaId : zones[k]) {
                            totalWorkload += testScenarios[s][areaId] * shortestPathDist[areaId][centerId];
                        }
                    }

                    for (int j = 0; j < zones.length; j++) {
                        if (zones[j] == null || zones[j].isEmpty()) {
                            continue;
                        }
                        int centerId = centers.get(j);
                        double districtWorkload = 0.0;
                        for (int areaId : zones[j]) {
                            districtWorkload += testScenarios[s][areaId] * shortestPathDist[areaId][centerId];
                        }
                        double lowerBound = coeffLower * totalWorkload;
                        double upperBound = coeffUpper * totalWorkload;
                        if (districtWorkload < lowerBound || districtWorkload > upperBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                } else {
                    double totalDemand = 0.0;
                    for (int i = 0; i < instance.getN(); i++) {
                        totalDemand += testScenarios[s][i];
                    }

                    for (int j = 0; j < zones.length; j++) {
                        if (zones[j] == null || zones[j].isEmpty()) {
                            continue;
                        }
                        double districtDemand = 0.0;
                        for (int areaId : zones[j]) {
                            districtDemand += testScenarios[s][areaId];
                        }
                        double lowerBound = coeffLower * totalDemand;
                        double upperBound = coeffUpper * totalDemand;
                        if (districtDemand < lowerBound || districtDemand > upperBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                }
            } else {
                if (bounds == null) {
                    return -1.0;
                }
                if (useImprovedModel && shortestPathDist != null) {
                    for (int j = 0; j < zones.length; j++) {
                        if (zones[j] == null || zones[j].isEmpty()) {
                            continue;
                        }
                        int centerId = centers.get(j);
                        double districtWorkload = 0.0;
                        for (int areaId : zones[j]) {
                            districtWorkload += testScenarios[s][areaId] * shortestPathDist[areaId][centerId];
                        }
                        if (districtWorkload < bounds.lowerBound || districtWorkload > bounds.upperBound) {
                            scenarioSatisfied = false;
                            break;
                        }
                    }
                } else {
                    for (int j = 0; j < zones.length; j++) {
                        if (zones[j] == null || zones[j].isEmpty()) {
                            continue;
                        }
                        double districtDemand = 0.0;
                        for (int areaId : zones[j]) {
                            districtDemand += testScenarios[s][areaId];
                        }
                        if (districtDemand < bounds.lowerBound || districtDemand > bounds.upperBound) {
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

        return (double) satisfiedScenarios / testScenarios.length;
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
