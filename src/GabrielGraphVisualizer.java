import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

/**
 * Gabriel Graph 可视化工具
 * 读取 dat 文件并生成图的可视化图片
 */
public class GabrielGraphVisualizer {
    private int n; // 节点数量
    private int k; // 区域数量
    private ArrayList<NodePoint> nodes; // 节点坐标
    private ArrayList<int[]> edges; // 边列表
    private int width = 1200;
    private int height = 1000;
    private double padding = 50;
    private double minX, maxX, minY, maxY;
    private double scaleX, scaleY;

    /**
     * 从 dat 文件读取图数据
     */
    public GabrielGraphVisualizer(String datFilePath) throws FileNotFoundException {
        File file = new File(datFilePath);
        Scanner sc = new Scanner(file);

        // 读取节点数量
        n = sc.nextInt();
        
        // 读取下一行，判断是 k 还是节点数据
        // 新格式：节点数 -> k -> 节点数据
        // 旧格式：节点数 -> 节点数据
        String nextLine = sc.nextLine().trim();
        String[] tokens = nextLine.split("\\s+");
        
        boolean hasK = false;
        if (tokens.length == 1 && tokens[0].matches("^\\d+$")) {
            // 只有一个整数，可能是 k（新格式）
            // 检查再下一行是否是节点数据（3个值：id x y）
            if (sc.hasNextLine()) {
                String peekLine = sc.nextLine().trim();
                String[] peekTokens = peekLine.split("\\s+");
                if (peekTokens.length == 3) {
                    // 下一行是3个值，说明刚才读的是 k（新格式）
                    k = Integer.parseInt(tokens[0]);
                    hasK = true;
                    // 注意：我们已经读取了 peekLine，这是第一行节点数据
                    // 需要重新创建 Scanner 来正确读取所有数据
                    sc.close();
                    sc = new Scanner(file);
                    sc.nextInt(); // 跳过节点数
                    sc.nextInt(); // 跳过 k
                } else {
                    // 下一行不是3个值，说明刚才读的不是 k，需要回退
                    sc.close();
                    sc = new Scanner(file);
                    sc.nextInt(); // 跳过节点数
                    k = Math.max(4, n / 15); // 使用默认值
                }
            } else {
                k = Integer.parseInt(tokens[0]);
                hasK = true;
            }
        } else {
            // 不是单个整数，说明是节点数据（旧格式）
            sc.close();
            sc = new Scanner(file);
            sc.nextInt(); // 跳过节点数
            k = Math.max(4, n / 15); // 使用默认值
        }

        nodes = new ArrayList<>();
        edges = new ArrayList<>();

        // 读取节点坐标
        minX = Double.MAX_VALUE;
        maxX = Double.MIN_VALUE;
        minY = Double.MAX_VALUE;
        maxY = Double.MIN_VALUE;

        for (int i = 0; i < n; i++) {
            int id = sc.nextInt();
            double x = sc.nextDouble();
            double y = sc.nextDouble();
            nodes.add(new NodePoint(id, x, y));

            minX = Math.min(minX, x);
            maxX = Math.max(maxX, x);
            minY = Math.min(minY, y);
            maxY = Math.max(maxY, y);
        }

        // 计算缩放因子
        double rangeX = maxX - minX;
        double rangeY = maxY - minY;
        if (rangeX == 0) rangeX = 1;
        if (rangeY == 0) rangeY = 1;

        scaleX = (width - 2 * padding) / rangeX;
        scaleY = (height - 2 * padding) / rangeY;
        // 使用统一的缩放因子以保持纵横比
        double scale = Math.min(scaleX, scaleY);
        scaleX = scale;
        scaleY = scale;

        // 读取边
        int m = sc.nextInt();
        for (int i = 0; i < m; i++) {
            int a = sc.nextInt();
            int b = sc.nextInt();
            edges.add(new int[]{a, b});
        }

        sc.close();
    }

    /**
     * 将世界坐标转换为屏幕坐标
     */
    private java.awt.geom.Point2D.Double worldToScreen(double x, double y) {
        double screenX = padding + (x - minX) * scaleX;
        double screenY = padding + (y - minY) * scaleY;
        // 翻转 Y 轴（因为屏幕坐标 Y 轴向下）
        screenY = height - padding - (y - minY) * scaleY;
        return new java.awt.geom.Point2D.Double(screenX, screenY);
    }

    /**
     * 生成可视化图片
     */
    public void visualize(String outputPath) throws IOException {
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = image.createGraphics();

        // 设置渲染质量
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        // 绘制背景
        g2d.setColor(new Color(255, 255, 255));
        g2d.fillRect(0, 0, width, height);

        // 绘制标题
        g2d.setColor(new Color(40, 40, 40));
        g2d.setFont(new Font("Arial", Font.BOLD, 18));
        String title = String.format("Gabriel Graph - %d Nodes, %d Edges, %d Districts", n, edges.size(), k);
        FontMetrics fm = g2d.getFontMetrics();
        int titleWidth = fm.stringWidth(title);
        g2d.drawString(title, (width - titleWidth) / 2, 30);

        // 绘制统计信息
        g2d.setFont(new Font("Arial", Font.PLAIN, 12));
        String stats = String.format("Average Degree: %.2f", 2.0 * edges.size() / n);
        g2d.drawString(stats, 20, height - 30);

        // 绘制边（先绘制边，这样节点会显示在边上面）
        g2d.setColor(new Color(200, 200, 200));
        g2d.setStroke(new BasicStroke(1.0f));
        for (int[] edge : edges) {
            NodePoint nodeA = nodes.get(edge[0]);
            NodePoint nodeB = nodes.get(edge[1]);
            java.awt.geom.Point2D.Double screenA = worldToScreen(nodeA.x, nodeA.y);
            java.awt.geom.Point2D.Double screenB = worldToScreen(nodeB.x, nodeB.y);
            g2d.draw(new Line2D.Double(screenA.x, screenA.y, screenB.x, screenB.y));
        }

        // 绘制节点
        double nodeRadius = 4.0;
        for (NodePoint node : nodes) {
            java.awt.geom.Point2D.Double screenPos = worldToScreen(node.x, node.y);

            // 绘制节点（圆形）
            g2d.setColor(new Color(45, 125, 210)); // 蓝色
            g2d.fill(new Ellipse2D.Double(
                    screenPos.x - nodeRadius,
                    screenPos.y - nodeRadius,
                    2 * nodeRadius,
                    2 * nodeRadius
            ));

            // 绘制节点边框
            g2d.setColor(new Color(20, 80, 150));
            g2d.setStroke(new BasicStroke(1.5f));
            g2d.draw(new Ellipse2D.Double(
                    screenPos.x - nodeRadius,
                    screenPos.y - nodeRadius,
                    2 * nodeRadius,
                    2 * nodeRadius
            ));
        }

        // 可选：绘制节点标签（对于小图）
        if (n <= 100) {
            g2d.setColor(new Color(60, 60, 60));
            g2d.setFont(new Font("Arial", Font.PLAIN, 9));
            for (NodePoint node : nodes) {
                java.awt.geom.Point2D.Double screenPos = worldToScreen(node.x, node.y);
                String label = String.valueOf(node.id);
                FontMetrics labelFm = g2d.getFontMetrics();
                int labelWidth = labelFm.stringWidth(label);
                g2d.drawString(label, (float)(screenPos.x - labelWidth / 2), (float)(screenPos.y - nodeRadius - 3));
            }
        }

        g2d.dispose();

        // 保存图片
        File outputFile = new File(outputPath);
        ImageIO.write(image, "png", outputFile);
        System.out.println("可视化图片已保存到: " + outputPath);
    }

    /**
     * 简单的节点类
     */
    private static class NodePoint {
        int id;
        double x, y;

        NodePoint(int id, double x, double y) {
            this.id = id;
            this.x = x;
            this.y = y;
        }
    }

    /**
     * 主函数 - 用于测试
     */
    public static void main(String[] args) {
        // 如果没有提供参数，可以使用默认文件名（方便在 IDEA 中直接运行）
        // 取消下面这行的注释并设置默认文件名即可：
         if (args.length == 0) args = new String[]{"GG500-10-3.dat"};
        
        if (args.length < 1) {
            System.out.println("用法: java GabrielGraphVisualizer <文件名或路径> [输出图片路径]");
            System.out.println();
            System.out.println("示例:");
            System.out.println("  1. 指定完整路径:");
            System.out.println("     java GabrielGraphVisualizer Instances_new/GG100-6-1.dat");
            System.out.println("  2. 只指定文件名（自动在 Instances_new 目录中查找）:");
            System.out.println("     java GabrielGraphVisualizer GG100-6-1.dat");
            System.out.println("  3. 指定输出路径:");
            System.out.println("     java GabrielGraphVisualizer GG100-6-1.dat output/my_graph.png");
            System.out.println();
            System.out.println("提示: 在 IDEA 中运行时，可以在 Run Configuration 的 Program arguments 中输入文件名");
            return;
        }

        String inputPath = args[0];
        String datFilePath;
        
        // 如果输入路径不包含目录分隔符，自动在 Instances_new 目录中查找
        if (!inputPath.contains("/") && !inputPath.contains("\\")) {
            datFilePath = "Instances_new/" + inputPath;
            System.out.println("自动查找文件: " + datFilePath);
        } else {
            datFilePath = inputPath;
        }
        
        // 检查文件是否存在
        File datFile = new File(datFilePath);
        if (!datFile.exists()) {
            System.err.println("错误: 文件不存在: " + datFilePath);
            System.err.println("请检查文件路径是否正确");
            return;
        }

        String outputPath;
        if (args.length >= 2) {
            outputPath = args[1];
        } else {
            // 默认输出路径
            String fileName = datFile.getName();
            String baseName = fileName.substring(0, fileName.lastIndexOf('.'));
            outputPath = "output/" + baseName + "_visualization.png";
        }

        try {
            // 确保输出目录存在
            File outputDir = new File(outputPath).getParentFile();
            if (outputDir != null && !outputDir.exists()) {
                outputDir.mkdirs();
            }

            System.out.println("正在读取文件: " + datFilePath);
            GabrielGraphVisualizer visualizer = new GabrielGraphVisualizer(datFilePath);
            visualizer.visualize(outputPath);
            System.out.println("可视化完成！");
        } catch (Exception e) {
            System.err.println("错误: " + e.getMessage());
            e.printStackTrace();
        }
    }
}

