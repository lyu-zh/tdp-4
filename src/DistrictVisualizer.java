import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class DistrictVisualizer {
    private Instance instance;
    private ArrayList<Integer>[] zones;
    private ArrayList<Area> centers;
    private int width = 1200;
    private int height = 800;
    private int mainX = 10; // 设置主区域的起始X坐标（不是0）
    private int mainY = 50;
    private int mainWidth = 980; // 调整宽度以适应新的X起始点
    private int mainHeight = 700;
    private int padding = 80;
    private double minX, maxX, minY, maxY;
    private double scaleX, scaleY;
    private Map<Integer, Point> areaPositions;

    public DistrictVisualizer(Instance instance, ArrayList<Integer>[] zones, ArrayList<Area> centers) {
        this.instance = instance;
        this.zones = zones;
        this.centers = centers;
        this.areaPositions = new HashMap<>();

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

        // Calculate scaling factors for initial positions
        scaleX = (mainWidth - 2 * padding) / (maxX - minX);
        scaleY = (mainHeight - 2 * padding) / (maxY - minY);

        // Generate initial positions - 添加mainX
        for (Area area : areas) {
            int x = (int) (mainX + padding + (area.getX() - minX) * scaleX);
            int y = (int) (mainY + padding + (area.getY() - minY) * scaleY);
            areaPositions.put(area.getId(), new Point(x, y));
        }

        // Adjust positions to avoid overlaps
        adjustPositionsToAvoidOverlaps();

        // Ensure all nodes are inside borders
        ensureNodesInsideBorders();
    }

    // Ensure all nodes are inside borders
    private void ensureNodesInsideBorders() {
        int buffer = 25; // Minimum distance from edge

        for (Point p : areaPositions.values()) {
            p.x = Math.max(mainX + buffer, Math.min(mainX + mainWidth - buffer, p.x));
            p.y = Math.max(mainY + buffer, Math.min(mainY + mainHeight - buffer, p.y));
        }
    }

    private void adjustPositionsToAvoidOverlaps() {
        Area[] areas = instance.getAreas();
        int nodeRadius = 8;
        int minDistance = nodeRadius * 4; // Minimum distance between nodes

        // Use force-directed positioning to separate overlapping nodes
        boolean hasOverlap = true;
        int iterations = 0;
        int maxIterations = 150; // Increased iterations for better spacing

        while (hasOverlap && iterations < maxIterations) {
            hasOverlap = false;
            iterations++;

            for (int i = 0; i < areas.length; i++) {
                int id1 = areas[i].getId();
                Point p1 = areaPositions.get(id1);

                for (int j = i + 1; j < areas.length; j++) {
                    int id2 = areas[j].getId();
                    Point p2 = areaPositions.get(id2);

                    // Calculate distance between nodes
                    double dx = p2.x - p1.x;
                    double dy = p2.y - p1.y;
                    double distance = Math.sqrt(dx * dx + dy * dy);

                    // If nodes are too close, push them apart
                    if (distance < minDistance) {
                        hasOverlap = true;

                        // Calculate repulsion force
                        double force = minDistance - distance;
                        double fx = force * dx / distance;
                        double fy = force * dy / distance;

                        // Apply force to both nodes in opposite directions
                        p2.x += (int) (fx / 2);
                        p2.y += (int) (fy / 2);
                        p1.x -= (int) (fx / 2);
                        p1.y -= (int) (fy / 2);

                        // Keep nodes within the main area
                        p1.x = Math.max(mainX + padding, Math.min(mainX + mainWidth - padding, p1.x));
                        p1.y = Math.max(mainY + padding, Math.min(mainY + mainHeight - padding, p1.y));
                        p2.x = Math.max(mainX + padding, Math.min(mainX + mainWidth - padding, p2.x));
                        p2.y = Math.max(mainY + padding, Math.min(mainY + mainHeight - padding, p2.y));
                    }
                }
            }
        }

        // After separation, also apply attraction forces between connected nodes
        for (int iter = 0; iter < 60; iter++) {
            for (Area area : areas) {
                int id1 = area.getId();
                Point p1 = areaPositions.get(id1);

                for (int neighborId : area.getNeighbors()) {
                    Point p2 = areaPositions.get(neighborId);

                    // Calculate distance
                    double dx = p2.x - p1.x;
                    double dy = p2.y - p1.y;
                    double distance = Math.sqrt(dx * dx + dy * dy);

                    // If connected nodes are too far, pull them closer
                    if (distance > minDistance * 4) {
                        double force = (distance - minDistance * 3) / 10;
                        double fx = force * dx / distance;
                        double fy = force * dy / distance;

                        // Apply gentle attraction
                        p1.x += (int) (fx);
                        p1.y += (int) (fy);
                        p2.x -= (int) (fx);
                        p2.y -= (int) (fy);

                        // Keep nodes within bounds
                        p1.x = Math.max(mainX + padding, Math.min(mainX + mainWidth - padding, p1.x));
                        p1.y = Math.max(mainY + padding, Math.min(mainY + mainHeight - padding, p1.y));
                        p2.x = Math.max(mainX + padding, Math.min(mainX + mainWidth - padding, p2.x));
                        p2.y = Math.max(mainY + padding, Math.min(mainY + mainHeight - padding, p2.y));
                    }
                }
            }
        }
    }

    public void saveVisualization(String outputPath) {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = image.createGraphics();

        // Set rendering hints for better quality
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        // Fill background
        g2d.setColor(new Color(248, 248, 252)); // Very light blue-gray background
        g2d.fillRect(0, 0, width, height);

        // Draw main area with subtle gradient
        GradientPaint gradient = new GradientPaint(
                mainX, mainY, new Color(255, 255, 255),
                mainX, mainY + mainHeight, new Color(240, 242, 245)
        );
        g2d.setPaint(gradient);
        g2d.fillRect(mainX, mainY, mainWidth, mainHeight);

        // Generate colors for each district
        Color[] districtColors = generateDistrictColors(zones.length);

        // Draw title with shadow effect
        g2d.setColor(new Color(60, 60, 80));
        g2d.setFont(new Font("Arial", Font.BOLD, 20));

        String title = "District Visualization - " + zones.length + " Districts";
        FontMetrics fm = g2d.getFontMetrics();
        int titleWidth = fm.stringWidth(title);

        g2d.drawString(title, width / 2 - titleWidth / 2, 30);

        // Draw border around main diagram - 明确设置边框的所有边
        g2d.setStroke(new BasicStroke(2.0f));
        g2d.setColor(new Color(80, 80, 100));
        Rectangle2D mainRect = new Rectangle2D.Double(mainX, mainY, mainWidth, mainHeight);
        g2d.draw(mainRect);

        // Draw legend (outside the main diagram)
        drawLegend(g2d, districtColors);

        // Only draw within the main area from now on
        g2d.clip(mainRect);

        // Draw connections between adjacent areas
        drawConnections(g2d);

        // Draw areas colored by district
        drawAreas(g2d, districtColors);

        // Highlight centers
        drawCenters(g2d);

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

        // Enhanced color palette with more vibrant, distinguishable colors
        Color[] palette = {
                new Color(45, 125, 210),   // blue
                new Color(255, 112, 10),   // orange
                new Color(35, 170, 85),    // green
                new Color(220, 55, 60),    // red
                new Color(150, 100, 205),  // purple
                new Color(128, 80, 60),    // brown
                new Color(230, 90, 183),   // pink
                new Color(80, 80, 90),     // gray
                new Color(190, 180, 30),   // olive
                new Color(20, 180, 195),   // cyan
                new Color(90, 165, 255),   // light blue
                new Color(255, 165, 70),   // light orange
                new Color(160, 220, 60),   // light green
                new Color(255, 75, 150)    // light red
        };

        for (int i = 0; i < numDistricts; i++) {
            if (i < palette.length) {
                colors[i] = palette[i];
            } else {
                // Generate more colors if we run out of palette colors
                float hue = (float) i / numDistricts;
                colors[i] = Color.getHSBColor(hue, 0.85f, 0.9f);
            }
        }

        return colors;
    }

    private void drawConnections(Graphics2D g2d) {
        g2d.setColor(new Color(200, 200, 210, 180));  // Slightly transparent gray
        g2d.setStroke(new BasicStroke(0.9f));

        Area[] areas = instance.getAreas();

        for (int i = 0; i < areas.length; i++) {
            Area area = areas[i];
            Point p1 = areaPositions.get(area.getId());

            for (int neighborId : area.getNeighbors()) {
                // Only draw each connection once
                if (neighborId > i) {
                    Point p2 = areaPositions.get(neighborId);
                    g2d.draw(new Line2D.Double(p1.x, p1.y, p2.x, p2.y));
                }
            }
        }
    }

    private void drawAreas(Graphics2D g2d, Color[] districtColors) {
        Area[] areas = instance.getAreas();
        int nodeRadius = 9;  // Slightly larger nodes

        // First draw all area nodes
        for (Area area : areas) {
            Point p = areaPositions.get(area.getId());

            // Find which district this area belongs to
            int districtIndex = -1;
            for (int j = 0; j < zones.length; j++) {
                if (zones[j].contains(area.getId())) {
                    districtIndex = j;
                    break;
                }
            }

            // Check if this is a center
            boolean isCenter = false;
            for (Area center : centers) {
                if (center.getId() == area.getId()) {
                    isCenter = true;
                    break;
                }
            }

            // Draw area node with district color
            if (districtIndex >= 0) {
                // Fill with district color
                g2d.setColor(districtColors[districtIndex]);
                g2d.fill(new Ellipse2D.Double(p.x - nodeRadius, p.y - nodeRadius,
                        2 * nodeRadius, 2 * nodeRadius));

                // Draw outline - thicker for centers
                if (isCenter) {
                    g2d.setColor(new Color(40, 40, 40));
                    g2d.setStroke(new BasicStroke(2.0f));
                } else {
                    g2d.setColor(new Color(60, 60, 60));
                    g2d.setStroke(new BasicStroke(1.0f));
                }
                g2d.draw(new Ellipse2D.Double(p.x - nodeRadius, p.y - nodeRadius,
                        2 * nodeRadius, 2 * nodeRadius));
            }

            // Draw area ID labels (only for non-centers)
            if (!isCenter) {
                g2d.setColor(new Color(40, 40, 40));
                g2d.setFont(new Font("Arial", Font.PLAIN, 10));
                String idStr = String.valueOf(area.getId());
                FontMetrics fm = g2d.getFontMetrics();
                int textWidth = fm.stringWidth(idStr);
                g2d.drawString(idStr, p.x - textWidth / 2, p.y + nodeRadius + 12);
            }
        }
    }

    private void drawCenters(Graphics2D g2d) {
        int centerRadius = 9;  // Match node radius

        for (Area center : centers) {
            Point p = areaPositions.get(center.getId());

            // Draw center marker (just a cross inside the node)
            g2d.setColor(new Color(40, 40, 40));
            g2d.setStroke(new BasicStroke(2.0f));

            // Draw cross
            g2d.drawLine(p.x - centerRadius / 2, p.y, p.x + centerRadius / 2, p.y);
            g2d.drawLine(p.x, p.y - centerRadius / 2, p.x, p.y + centerRadius / 2);

            // Add center ID label
            g2d.setFont(new Font("Arial", Font.BOLD, 11));
            String idStr = "C" + center.getId();
            FontMetrics fm = g2d.getFontMetrics();
            int textWidth = fm.stringWidth(idStr);
            g2d.drawString(idStr, p.x - textWidth / 2, p.y - centerRadius - 5);
        }
    }

    private void drawLegend(Graphics2D g2d, Color[] districtColors) {
        int legendX = mainX + mainWidth + 20;
        int legendY = 70;
        int itemHeight = 25;
        int colorBoxSize = 16;

        // Draw legend background with rounded corners
        g2d.setColor(new Color(248, 248, 252));
        g2d.fillRoundRect(legendX - 10, legendY - 20, width - (mainX + mainWidth) - 30,
                zones.length * itemHeight + 120, 15, 15);

        g2d.setColor(new Color(220, 220, 230));
        g2d.setStroke(new BasicStroke(1.0f));
        g2d.drawRoundRect(legendX - 10, legendY - 20, width - (mainX + mainWidth) - 30,
                zones.length * itemHeight + 120, 15, 15);

        // Draw legend title
        g2d.setFont(new Font("Arial", Font.BOLD, 14));
        g2d.setColor(new Color(50, 50, 70));
        g2d.drawString("Districts", legendX, legendY);

        // Draw legend items
        for (int i = 0; i < zones.length; i++) {
            // Draw color box with rounded corners
            g2d.setColor(districtColors[i]);
            g2d.fillRoundRect(legendX, legendY + 10 + i * itemHeight, colorBoxSize, colorBoxSize, 4, 4);

            // Draw outline
            g2d.setColor(new Color(100, 100, 120));
            g2d.drawRoundRect(legendX, legendY + 10 + i * itemHeight, colorBoxSize, colorBoxSize, 4, 4);

            // Draw district label
            g2d.setFont(new Font("Arial", Font.PLAIN, 12));
            g2d.setColor(new Color(40, 40, 60));
            g2d.drawString("District " + (i + 1) + " (C" + centers.get(i).getId() + ")",
                    legendX + colorBoxSize + 8, legendY + 10 + i * itemHeight + 12);
        }

        // Draw divider line
        g2d.setColor(new Color(220, 220, 230));
        g2d.drawLine(legendX - 5, legendY + 15 + zones.length * itemHeight,
                legendX + width - (mainX + mainWidth) - 35, legendY + 15 + zones.length * itemHeight);

        // Draw node explanation
        int explanationY = legendY + 25 + zones.length * itemHeight;
        g2d.setFont(new Font("Arial", Font.BOLD, 14));
        g2d.setColor(new Color(50, 50, 70));
        g2d.drawString("Node Types", legendX, explanationY);

        // Regular node example
        g2d.setFont(new Font("Arial", Font.PLAIN, 12));
        int nodeRadius = 8;
        int nodeY = explanationY + 25;

        // Draw sample regular node
        g2d.setColor(new Color(180, 180, 210));
        g2d.fill(new Ellipse2D.Double(legendX + 8 - nodeRadius, nodeY - nodeRadius,
                2 * nodeRadius, 2 * nodeRadius));
        g2d.setColor(new Color(60, 60, 60));
        g2d.setStroke(new BasicStroke(1.0f));
        g2d.draw(new Ellipse2D.Double(legendX + 8 - nodeRadius, nodeY - nodeRadius,
                2 * nodeRadius, 2 * nodeRadius));

        g2d.setColor(new Color(40, 40, 60));
        g2d.drawString("Regular Area", legendX + 25, nodeY + 4);

        // Center node example
        int centerY = nodeY + 30;

        // Draw sample center node
        g2d.setColor(new Color(180, 180, 210));
        g2d.fill(new Ellipse2D.Double(legendX + 8 - nodeRadius, centerY - nodeRadius,
                2 * nodeRadius, 2 * nodeRadius));
        g2d.setColor(new Color(40, 40, 40));
        g2d.setStroke(new BasicStroke(2.0f));
        g2d.draw(new Ellipse2D.Double(legendX + 8 - nodeRadius, centerY - nodeRadius,
                2 * nodeRadius, 2 * nodeRadius));

        // Draw cross
        g2d.drawLine(legendX + 8 - nodeRadius / 2, centerY, legendX + 8 + nodeRadius / 2, centerY);
        g2d.drawLine(legendX + 8, centerY - nodeRadius / 2, legendX + 8, centerY + nodeRadius / 2);

        g2d.setColor(new Color(40, 40, 60));
        g2d.drawString("District Center", legendX + 25, centerY + 4);
    }

    // Static method to visualize from files
    public static void visualizeDistrictSolution(String instanceFile, String solutionFile, String outputImagePath) {
        try {
            // Load instance
            Instance instance = new Instance(instanceFile);

            // Parse solution file to extract zones and centers
            ArrayList<Integer>[] zones = parseZonesFromSolutionFile(solutionFile, instance.k);
            ArrayList<Area> centers = parseCentersFromSolutionFile(solutionFile, instance);

            // Create visualizer and generate image
            DistrictVisualizer visualizer = new DistrictVisualizer(instance, zones, centers);
            visualizer.saveVisualization(outputImagePath);

        } catch (Exception e) {
            System.err.println("Error visualizing district solution: " + e.getMessage());
            e.printStackTrace();
        }
    }

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
                if (currentZone >= k) break; // Safety check
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