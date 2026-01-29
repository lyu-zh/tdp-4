import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class Instance {
    private int n; // 区域个数
    private Area[] areas; // 区域信息
    private int[][] edges; // 边集合

    public double average1;
    public double average2;
    public int k; // 需要分成的大区域个数
    double[][] dist;
    
    // 地球半径（公里），用于Haversine公式计算经纬度距离
    private static final double EARTH_RADIUS_KM = 6371.0;



    public Instance(String filepath) throws FileNotFoundException {
        File file = new File(filepath);
        Scanner sc = new Scanner(file);

        n = sc.nextInt(); // 节点个数
        
        // 读取区域数量 k
        k = sc.nextInt();

        areas = new Area[n]; // 存储区域信息
        double sum = 0.0;
        double sum2 = 0.0;
        
        // 读取所有节点数据（新格式：只有 id x y）
        for (int i = 0; i < n; i++) {
            int id = sc.nextInt();
            double x = sc.nextDouble();
            double y = sc.nextDouble();

            // 新格式没有 activeness 数据，使用默认值
            double[] activeness = new double[3];
            activeness[0] = 500.0; // 默认基础需求
            activeness[1] = 1500.0; // 默认值
            activeness[2] = 1500.0; // 默认值
            
            sum += activeness[0];
            sum2 += activeness[1];
            areas[i] = new Area(id, x, y, activeness);
        }

        edges = new int[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                edges[i][j] = 0;
            }
        }

        // 读取边数
        int m = sc.nextInt();
        // 读取边数据
        for (int i = 0; i < m; i++) {
            int a = sc.nextInt();
            int b = sc.nextInt();
            areas[a].addNeighbor(b);
            areas[b].addNeighbor(a);
            edges[a][b] = 1;
            edges[b][a] = 1;
        }
        
        average1 = sum / k;
        average2 = sum2 / k;
        sc.close();

        // 生成距离矩阵
        // 判断是否为经纬度格式：如果文件路径包含"unique_coordinates_list_filtered"或"filtered"，则使用Haversine公式
        boolean useHaversine = filepath.contains("unique_coordinates_list_filtered") || 
                               filepath.contains("filtered");
        
        dist = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double distance;
                if (useHaversine) {
                    // 使用Haversine公式计算经纬度距离（单位：公里）
                    // x是经度(longitude)，y是纬度(latitude)
                    distance = haversineDistance(areas[i].getX(), areas[i].getY(), 
                                                 areas[j].getX(), areas[j].getY());
                } else {
                    // 使用欧几里得距离（适用于平面坐标）
                    distance = Math.sqrt(Math.pow(areas[i].getX() - areas[j].getX(), 2) + 
                                         Math.pow(areas[i].getY() - areas[j].getY(), 2));
                }
                dist[i][j] = distance;
                dist[j][i] = dist[i][j];
            }
        }
        
        if (useHaversine) {
            System.out.println("已使用Haversine公式计算经纬度距离（单位：公里）");
        }
    }

    public Instance() {

    }

    /**
     * 创建当前实例的深度克隆
     * @return 克隆的实例对象
     */
    public Instance clone() {
        Instance cloned = new Instance();
        cloned.n = this.n;
        cloned.k = this.k;
        cloned.average1 = this.average1;
        cloned.average2 = this.average2;

        // 复制区域信息
        cloned.areas = new Area[this.n];
        for (int i = 0; i < this.n; i++) {
            cloned.areas[i] = new Area();
            cloned.areas[i].setId(this.areas[i].getId());
            cloned.areas[i].setX(this.areas[i].getX());
            cloned.areas[i].setY(this.areas[i].getY());

            // 复制活跃度指标
            double[] originalActiveness = this.areas[i].getActiveness();
            double[] newActiveness = new double[originalActiveness.length];
            System.arraycopy(originalActiveness, 0, newActiveness, 0, originalActiveness.length);
            cloned.areas[i].setActiveness(newActiveness);

            // 复制中心标识
            cloned.areas[i].setCenter(this.areas[i].isCenter());

            // 复制邻居列表
            for (int neighbor : this.areas[i].getNeighbors()) {
                cloned.areas[i].addNeighbor(neighbor);
            }
        }

        // 复制边信息
        cloned.edges = new int[this.n][this.n];
        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                cloned.edges[i][j] = this.edges[i][j];
            }
        }

        // 复制距离矩阵
        cloned.dist = new double[this.n][this.n];
        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                cloned.dist[i][j] = this.dist[i][j];
            }
        }

        return cloned;
    }

    /**
     * 基于原始实例创建新实例，使用特定场景的需求
     * @param original 原始实例
     * @param scenarioDemands 场景需求数组
     */
    public Instance(Instance original, int[] scenarioDemands) {
        this.n = original.n;
        this.k = original.k;
        this.average1 = original.average1;
        this.average2 = original.average2;

        // 复制区域信息，但使用新的需求值
        this.areas = new Area[this.n];
        double sum1 = 0.0;
        double sum2 = 0.0;

        for (int i = 0; i < this.n; i++) {
            // 使用带参数的构造函数，确保neighbors被初始化
            double[] originalActiveness = original.areas[i].getActiveness(); //拷贝第i个基本单元的原始活跃度数组
            double[] newActiveness = new double[originalActiveness.length]; //创建新的活跃度数组
            System.arraycopy(originalActiveness, 0, newActiveness, 0, originalActiveness.length);

            // 更新第一个活跃度指标为场景需求
            newActiveness[0] = scenarioDemands[i];

            // 使用带参数的构造函数，而不是无参构造函数
            this.areas[i] = new Area(original.areas[i].getId(),
                    original.areas[i].getX(),
                    original.areas[i].getY(),
                    newActiveness);

            sum1 += newActiveness[0];
            sum2 += newActiveness[1];

            // 复制中心标识
            this.areas[i].setCenter(original.areas[i].isCenter());

            // 复制邻居列表 - neighbors现在已经初始化了
            for (int neighbor : original.areas[i].getNeighbors()) {
                this.areas[i].addNeighbor(neighbor);
            }
        }

        // 更新平均值
        this.average1 = sum1 / this.k;
        this.average2 = sum2 / this.k;

        // 复制边信息
        this.edges = new int[this.n][this.n];
        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                this.edges[i][j] = original.edges[i][j];
            }
        }

        // 复制距离矩阵
        this.dist = new double[this.n][this.n];
        for (int i = 0; i < this.n; i++) {
            for (int j = 0; j < this.n; j++) {
                this.dist[i][j] = original.dist[i][j];
            }
        }
    }

    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }

    public Area[] getAreas() {
        return areas;
    }

    public void setAreas(Area[] areas) {
        this.areas = areas;
    }

    public int[][] getEdges() {
        return edges;
    }

    public void output() {
        System.out.println("区域个数: " + n);
        System.out.println("区域信息:");
        for (Area area : areas) {
            System.out.println("区域编号：" + area.getId() + "，横坐标：" + area.getX() + "，纵坐标：" + area.getY() +
                    "，活跃度指标：{" + area.getActiveness()[0] + ", " + area.getActiveness()[1] + ", " + area.getActiveness()[2] + "}");
        }
        System.out.println("边集合:");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print(edges[i][j] + " ");
            }
            System.out.println();
        }
    }


    public void setEdges(int[][] edges) {
        this.edges = edges;
    }
    
    /**
     * 使用Haversine公式计算两点间距离（公里）
     * @param lon1 第一个点的经度
     * @param lat1 第一个点的纬度
     * @param lon2 第二个点的经度
     * @param lat2 第二个点的纬度
     * @return 两点间的距离（公里）
     */
    private static double haversineDistance(double lon1, double lat1, double lon2, double lat2) {
        double lat1Rad = Math.toRadians(lat1);
        double lat2Rad = Math.toRadians(lat2);
        double deltaLatRad = Math.toRadians(lat2 - lat1);
        double deltaLonRad = Math.toRadians(lon2 - lon1);
        
        double a = Math.sin(deltaLatRad / 2) * Math.sin(deltaLatRad / 2) +
                   Math.cos(lat1Rad) * Math.cos(lat2Rad) *
                   Math.sin(deltaLonRad / 2) * Math.sin(deltaLonRad / 2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        
        return EARTH_RADIUS_KM * c;
    }
}


class Area {
    private int id; // 区域编号
    private double x; // 横坐标
    private double y; // 纵坐标
    private double[] activeness; // 活跃度指标，是一个数组，每一个代表不同的活跃度
    private boolean isCenter = false; // 是否是大区域中心
    private ArrayList<Integer> neighbors; // 存储所有相邻区域的编号

    public Area() {
        this.neighbors = new ArrayList<>();
    }

    public Area(int id, double x, double y, double[] activeness) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.activeness = activeness;
        this.isCenter = false; // 默认不是大区域中心
        this.neighbors = new ArrayList<>();
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public double[] getActiveness() {
        return activeness;
    }

    public void setActiveness(double[] activeness) {
        this.activeness = activeness;
    }

    public boolean isCenter() {
        return isCenter;
    }

    public void setCenter(boolean center) {
        isCenter = center;
    }

    public ArrayList<Integer> getNeighbors() {
        return neighbors;
    }

    public void addNeighbor(int neighborId) {
        neighbors.add(neighborId);
    }
}


