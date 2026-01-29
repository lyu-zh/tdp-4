import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.*;

public class DistrictPolygonVisualizer {
    private Instance instance;
    private ArrayList<Integer>[] zones;
    private ArrayList<Area> centers;
    private int width = 1200;
    private int height = 1000;
    private double padding = 50;
    private double minX, maxX, minY, maxY;
    private double scaleX, scaleY;
    private Map<Integer, Polygon2D> areaPolygons;

    public DistrictPolygonVisualizer(Instance instance, ArrayList<Integer>[] zones, ArrayList<Area> centers) {
        this.instance = instance;
        this.zones = zones;
        this.centers = centers;
        this.areaPolygons = new HashMap<>();

        // Calculate bounds for scaling
        Area[] areas = instance.getAreas();
        minX = Double.MAX_VALUE;
        maxX = Double.MIN_VALUE;
        minY = Double.MAX_VALUE;
        maxY = Double.MIN_VALUE;

        for (Area area : areas) {
            minX = Math.min(minX, area.getX());
            maxX = Math.max(maxX, area.getX());
            minY = Math.min(minY, area.getY());
            maxY = Math.max(maxY, area.getY());
        }

        // Calculate scaling factors
        scaleX = (width - 2 * padding) / (maxX - minX);
        scaleY = (height - 2 * padding) / (maxY - minY);

        // Generate polygons for all areas
        generateAreaPolygons();
    }

    private void generateAreaPolygons() {
        Area[] areas = instance.getAreas();

        // Step 1: Create a graph representation
        Map<Integer, Set<Integer>> adjacencyGraph = new HashMap<>();
        for (int i = 0; i < areas.length; i++) {
            Area area = areas[i];
            adjacencyGraph.put(area.getId(), new HashSet<>(area.getNeighbors()));
        }

        // Step 2: Use a force-directed layout to position the areas
        Map<Integer, Point2D.Double> areaPositions = forceDirectedLayout(areas, adjacencyGraph);

        // Step 3: Generate Voronoi-like polygons for each area
        for (Area area : areas) {
            int areaId = area.getId();
            Point2D.Double pos = areaPositions.get(areaId);

            // Get neighbors
            Set<Integer> neighbors = adjacencyGraph.get(areaId);
            List<Point2D.Double> neighborPositions = new ArrayList<>();

            for (int neighborId : neighbors) {
                neighborPositions.add(areaPositions.get(neighborId));
            }

            // Generate polygon for this area
            Polygon2D polygon = generateAreaPolygon(pos, neighborPositions, areas.length);
            areaPolygons.put(areaId, polygon);
        }
    }

    private Map<Integer, Point2D.Double> forceDirectedLayout(Area[] areas, Map<Integer, Set<Integer>> adjacencyGraph) {
        Map<Integer, Point2D.Double> positions = new HashMap<>();

        // Initialize positions using actual coordinates but scaled
        for (Area area : areas) {
            double x = padding + (area.getX() - minX) * scaleX;
            double y = height - padding - (area.getY() - minY) * scaleY;
            positions.put(area.getId(), new Point2D.Double(x, y));
        }

        // Perform force-directed layout to position nodes
        int iterations = 100;
        double k = 30.0; // Optimal distance
        double temperature = 0.1 * Math.min(width, height);

        for (int i = 0; i < iterations; i++) {
            // Calculate repulsive forces (nodes repel each other)
            Map<Integer, Point2D.Double> forces = new HashMap<>();
            for (Area area : areas) {
                forces.put(area.getId(), new Point2D.Double(0, 0));
            }

            for (int a = 0; a < areas.length; a++) {
                int idA = areas[a].getId();
                Point2D.Double posA = positions.get(idA);

                for (int b = a + 1; b < areas.length; b++) {
                    int idB = areas[b].getId();
                    Point2D.Double posB = positions.get(idB);

                    // Calculate displacement
                    double dx = posB.x - posA.x;
                    double dy = posB.y - posA.y;
                    double distance = Math.max(0.1, Math.sqrt(dx * dx + dy * dy));

                    // Repulsive force is inversely proportional to distance
                    double force = k * k / distance;
                    double fx = force * dx / distance;
                    double fy = force * dy / distance;

                    Point2D.Double forceA = forces.get(idA);
                    Point2D.Double forceB = forces.get(idB);

                    forceA.x -= fx;
                    forceA.y -= fy;
                    forceB.x += fx;
                    forceB.y += fy;
                }
            }

            // Calculate attractive forces (edges pull connected nodes together)
            for (int a = 0; a < areas.length; a++) {
                int idA = areas[a].getId();
                Point2D.Double posA = positions.get(idA);

                for (int neighborId : adjacencyGraph.get(idA)) {
                    Point2D.Double posB = positions.get(neighborId);

                    // Calculate displacement
                    double dx = posB.x - posA.x;
                    double dy = posB.y - posA.y;
                    double distance = Math.max(0.1, Math.sqrt(dx * dx + dy * dy));

                    // Attractive force is proportional to distance
                    double force = distance * distance / k;
                    double fx = force * dx / distance;
                    double fy = force * dy / distance;

                    Point2D.Double forceA = forces.get(idA);
                    forceA.x += fx;
                    forceA.y += fy;
                }
            }

            // Apply forces with temperature
            for (Area area : areas) {
                int id = area.getId();
                Point2D.Double pos = positions.get(id);
                Point2D.Double force = forces.get(id);

                double magnitude = Math.sqrt(force.x * force.x + force.y * force.y);
                if (magnitude > 0) {
                    double limitedMagnitude = Math.min(magnitude, temperature);
                    pos.x += force.x / magnitude * limitedMagnitude;
                    pos.y += force.y / magnitude * limitedMagnitude;

                    // Keep within bounds
                    pos.x = Math.max(padding, Math.min(width - padding, pos.x));
                    pos.y = Math.max(padding, Math.min(height - padding, pos.y));
                }
            }

            // Cool down
            temperature *= 0.95;
        }

        return positions;
    }

    private Polygon2D generateAreaPolygon(Point2D.Double center, List<Point2D.Double> neighborPositions, int totalAreas) {
        // If no neighbors, create a circle
        if (neighborPositions.isEmpty()) {
            double radius = 20.0;
            return createRegularPolygon(center, radius, 8);
        }

        // Create polygon points using Voronoi-like cell construction
        List<Point2D.Double> polygonPoints = new ArrayList<>();

        // Add points between this area and each neighbor
        for (int i = 0; i < neighborPositions.size(); i++) {
            Point2D.Double neighborPos = neighborPositions.get(i);
            Point2D.Double nextNeighborPos = neighborPositions.get((i + 1) % neighborPositions.size());

            // Midpoint between center and neighbor
            Point2D.Double midpoint = new Point2D.Double(
                    (center.x + neighborPos.x) / 2,
                    (center.y + neighborPos.y) / 2
            );

            // Add a point rotated around the midpoint
            double angle = Math.atan2(neighborPos.y - center.y, neighborPos.x - center.x);
            angle += Math.PI / 2; // Rotate 90 degrees

            double distance = 0.3 * center.distance(neighborPos);
            Point2D.Double rotatedPoint = new Point2D.Double(
                    midpoint.x + distance * Math.cos(angle),
                    midpoint.y + distance * Math.sin(angle)
            );

            polygonPoints.add(rotatedPoint);

            // If there are at least 3 neighbors, add an extra point to round the polygon
            if (neighborPositions.size() >= 3) {
                // Calculate a point between this neighbor and the next
                Point2D.Double between = new Point2D.Double(
                        (neighborPos.x + nextNeighborPos.x + center.x) / 3,
                        (neighborPos.y + nextNeighborPos.y + center.y) / 3
                );
                polygonPoints.add(between);
            }
        }

        // Sort the points in counterclockwise order around the center
        final Point2D.Double finalCenter = center;
        polygonPoints.sort((p1, p2) -> {
            double angle1 = Math.atan2(p1.y - finalCenter.y, p1.x - finalCenter.x);
            double angle2 = Math.atan2(p2.y - finalCenter.y, p2.x - finalCenter.x);
            return Double.compare(angle1, angle2);
        });

        return new Polygon2D(polygonPoints);
    }

    private Polygon2D createRegularPolygon(Point2D.Double center, double radius, int numSides) {
        List<Point2D.Double> points = new ArrayList<>();
        for (int i = 0; i < numSides; i++) {
            double angle = 2 * Math.PI * i / numSides;
            points.add(new Point2D.Double(
                    center.x + radius * Math.cos(angle),
                    center.y + radius * Math.sin(angle)
            ));
        }
        return new Polygon2D(points);
    }

    public void saveVisualization(String outputPath) {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = image.createGraphics();

        // Set rendering hints for better quality
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);

        // Fill background
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, width, height);

        // Generate colors for each district
        Color[] districtColors = generateDistrictColors(zones.length);

        // Draw polygons with district colors
        drawDistrictPolygons(g2d, districtColors);

        // Draw borders between areas
        drawAreaBorders(g2d);

        // Mark centers with distinctive marking
        drawCenters(g2d);

        // Add legend
        drawLegend(g2d, districtColors);

        // Add title
        g2d.setColor(Color.BLACK);
        g2d.setFont(new Font("Arial", Font.BOLD, 16));
        g2d.drawString("District Visualization - " + zones.length + " districts", width / 2 - 150, 30);

        g2d.dispose();

        try {
            File outputFile = new File(outputPath);
            ImageIO.write(image, "png", outputFile);
            System.out.println("Visualization saved to: " + outputPath);
        } catch (IOException e) {
            System.err.println("Error saving visualization: " + e.getMessage());
        }
    }

    private Color[] generateDistrictColors(int numDistricts) {
        Color[] colors = new Color[numDistricts];
        Random random = new Random(42); // Fixed seed for consistent colors

        for (int i = 0; i < numDistricts; i++) {
            // Generate distinct colors with good contrast
            float hue = (float) i / numDistricts;
            float saturation = 0.7f + random.nextFloat() * 0.3f; // 0.7-1.0
            float brightness = 0.8f + random.nextFloat() * 0.2f; // 0.8-1.0

            colors[i] = Color.getHSBColor(hue, saturation, brightness);
        }

        return colors;
    }

    private void drawDistrictPolygons(Graphics2D g2d, Color[] districtColors) {
        Area[] areas = instance.getAreas();

        // Draw each area with its district color
        for (int areaId = 0; areaId < areas.length; areaId++) {
            // Find which district this area belongs to
            int districtIndex = -1;
            for (int j = 0; j < zones.length; j++) {
                if (zones[j].contains(areaId)) {
                    districtIndex = j;
                    break;
                }
            }

            if (districtIndex >= 0 && areaPolygons.containsKey(areaId)) {
                Polygon2D polygon = areaPolygons.get(areaId);

                // Fill polygon with district color
                g2d.setColor(districtColors[districtIndex]);
                g2d.fill(polygon.toAWTPath());

                // Draw area ID label
                Point2D.Double center = polygon.getCenter();
                g2d.setColor(Color.BLACK);
                g2d.setFont(new Font("Arial", Font.PLAIN, 10));
                g2d.drawString(String.valueOf(areaId), (float) center.x - 4, (float) center.y + 4);
            }
        }
    }

    private void drawAreaBorders(Graphics2D g2d) {
        g2d.setColor(Color.BLACK);
        g2d.setStroke(new BasicStroke(0.5f));

        for (Polygon2D polygon : areaPolygons.values()) {
            g2d.draw(polygon.toAWTPath());
        }
    }

    private void drawCenters(Graphics2D g2d) {
        for (Area center : centers) {
            if (areaPolygons.containsKey(center.getId())) {
                Polygon2D polygon = areaPolygons.get(center.getId());
                Point2D.Double centerPoint = polygon.getCenter();

                double centerRadius = 12.0;

                // Draw center with a star symbol
                g2d.setColor(Color.BLACK);
                g2d.setStroke(new BasicStroke(2.0f));

                g2d.draw(new Ellipse2D.Double(
                        centerPoint.x - centerRadius,
                        centerPoint.y - centerRadius,
                        2 * centerRadius,
                        2 * centerRadius
                ));

                g2d.drawLine(
                        (int) (centerPoint.x - centerRadius),
                        (int) centerPoint.y,
                        (int) (centerPoint.x + centerRadius),
                        (int) centerPoint.y
                );

                g2d.drawLine(
                        (int) centerPoint.x,
                        (int) (centerPoint.y - centerRadius),
                        (int) centerPoint.x,
                        (int) (centerPoint.y + centerRadius)
                );

                // Add center ID
                g2d.setFont(new Font("Arial", Font.BOLD, 12));
                g2d.drawString("C" + center.getId(),
                        (float) (centerPoint.x + centerRadius + 2),
                        (float) (centerPoint.y - 2));
            }
        }
    }

    private void drawLegend(Graphics2D g2d, Color[] districtColors) {
        int legendX = width - 180;
        int legendY = 70;
        int itemHeight = 20;
        int colorBoxWidth = 15;

        g2d.setFont(new Font("Arial", Font.BOLD, 12));
        g2d.setColor(Color.BLACK);
        g2d.drawString("Districts:", legendX, legendY - 10);

        for (int i = 0; i < zones.length; i++) {
            // Draw color box
            g2d.setColor(districtColors[i]);
            g2d.fillRect(legendX, legendY + i * itemHeight, colorBoxWidth, colorBoxWidth);

            // Draw outline
            g2d.setColor(Color.BLACK);
            g2d.drawRect(legendX, legendY + i * itemHeight, colorBoxWidth, colorBoxWidth);

            // Draw district label
            g2d.drawString("District " + (i + 1) + " (C" + centers.get(i).getId() + ")",
                    legendX + colorBoxWidth + 5, legendY + i * itemHeight + 12);
        }
    }

    // Helper class to represent a polygon
    private static class Polygon2D {
        private List<Point2D.Double> points;

        public Polygon2D(List<Point2D.Double> points) {
            this.points = points;
        }

        public Path2D.Double toAWTPath() {
            Path2D.Double path = new Path2D.Double();

            if (!points.isEmpty()) {
                path.moveTo(points.get(0).x, points.get(0).y);
                for (int i = 1; i < points.size(); i++) {
                    path.lineTo(points.get(i).x, points.get(i).y);
                }
                path.closePath();
            }

            return path;
        }

        public Point2D.Double getCenter() {
            double sumX = 0, sumY = 0;
            for (Point2D.Double point : points) {
                sumX += point.x;
                sumY += point.y;
            }
            return new Point2D.Double(sumX / points.size(), sumY / points.size());
        }
    }

    public static void visualizeDistrictSolution(String instanceFile, String solutionFile, String outputImagePath) {
        try {
            // Load instance
            Instance instance = new Instance(instanceFile);

            // Parse solution file to extract zones and centers
            ArrayList<Integer>[] zones = parseZonesFromSolutionFile(solutionFile, instance.k);
            ArrayList<Area> centers = parseCentersFromSolutionFile(solutionFile, instance);

            // Create visualizer and generate image
            DistrictPolygonVisualizer visualizer = new DistrictPolygonVisualizer(instance, zones, centers);
            visualizer.saveVisualization(outputImagePath);

        } catch (Exception e) {
            System.err.println("Error visualizing district solution: " + e.getMessage());
            e.printStackTrace();
        }
    }

    // Note: You'll need to implement these methods based on your output file format
    @SuppressWarnings("unchecked")
    private static ArrayList<Integer>[] parseZonesFromSolutionFile(String solutionFile, int k) throws IOException {
        ArrayList<Integer>[] zones = new ArrayList[k];
        for (int i = 0; i < k; i++) {
            zones[i] = new ArrayList<>();
        }

        BufferedReader reader = new BufferedReader(new FileReader(solutionFile));
        String line;
        int currentZone = -1;

        while ((line = reader.readLine()) != null) {
            if (line.startsWith("center ID:")) {
                currentZone++;
            } else if (currentZone >= 0 && !line.startsWith("best") && !line.startsWith("程序")) {
                String[] parts = line.trim().split("\\s+");
                for (String part : parts) {
                    if (!part.isEmpty()) {
                        try {
                            int areaId = Integer.parseInt(part);
                            zones[currentZone].add(areaId);
                        } catch (NumberFormatException e) {
                            // Skip non-numeric entries
                        }
                    }
                }
            }
        }
        reader.close();

        return zones;
    }

    private static ArrayList<Area> parseCentersFromSolutionFile(String solutionFile, Instance instance) throws IOException {
        ArrayList<Area> centers = new ArrayList<>();

        BufferedReader reader = new BufferedReader(new FileReader(solutionFile));
        String line;

        while ((line = reader.readLine()) != null) {
            if (line.startsWith("center ID:")) {
                String idStr = line.substring("center ID:".length()).trim();
                try {
                    int centerId = Integer.parseInt(idStr);
                    centers.add(instance.getAreas()[centerId]);
                } catch (NumberFormatException | ArrayIndexOutOfBoundsException e) {
                    System.err.println("Error parsing center ID: " + idStr);
                }
            }
        }
        reader.close();

        return centers;
    }
}