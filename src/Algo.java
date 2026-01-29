import com.gurobi.gurobi.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Algo {
    private Instance inst;
    private ArrayList<Area> centers; // 存储所有大区域中心的坐标
    private ArrayList<Integer>[] zones;
    private double r;
    private int timeLimit = Integer.MAX_VALUE; // 默认无时间限制
    private double demandUpperBound; // 新增变量存储U
    private double demandLowerBound; // 新增变量存储L

    // 改进模型相关字段
    private boolean useImprovedModel = false; // 是否使用改进模型
    private double[][] shortestPathDist = null; // 最短路径距离矩阵（改进模型使用）
    private double[] meanVector = null; // 均值向量（改进模型使用）
    private boolean useRelativeBalance = false; // 是否使用相对平衡性约束

    // 设置时间限制
    public void setTimeLimit(int seconds) {
        this.timeLimit = seconds;
    }

    // 设置改进模型相关参数
    public void setImprovedModelParams(boolean useImprovedModel, double[][] shortestPathDist, double[] meanVector) {
        this.useImprovedModel = useImprovedModel;
        this.shortestPathDist = shortestPathDist;
        this.meanVector = meanVector;
    }

    // 设置边界值（用于改进模型，从外部传入已计算好的边界值）
    public void setBounds(double lowerBound, double upperBound) {
        this.demandLowerBound = lowerBound;
        this.demandUpperBound = upperBound;
    }

    // 为改进模型自动计算并设置基于工作量的边界值（使用实例的实际需求）
    // 根据新方法：只考虑每个基本单元的 k 个最近中心计算平均配送距离
    public void calculateAndSetWorkloadBounds(double[][] shortestPathDist, double r) {
        int n = inst.getN();
        int p = inst.k;
        int k = (int) Math.ceil((double) n / p); // k = ⌈n/p⌉
        k = Math.min(k, n); // 确保 k 不超过总节点数

        double totalWorkload = 0.0;

        // 计算每个基本单元的平均配送距离（只考虑 k 个最近中心）
        for (int i = 0; i < n; i++) {
            // 找到基本单元 i 到所有潜在中心的距离，并排序
            double[] distances = new double[n];
            for (int j = 0; j < n; j++) {
                distances[j] = shortestPathDist[i][j];
            }

            // 使用部分排序找到 k 个最小的距离
            // 简单方法：排序并取前 k 个
            java.util.Arrays.sort(distances);

            // 计算 k 个最近中心的平均距离
            double sum = 0.0;
            for (int idx = 0; idx < k; idx++) {
                sum += distances[idx];
            }
            double avgDist = sum / k;

            // 使用实例的实际需求计算工作量
            double demand = inst.getAreas()[i].getActiveness()[0];
            totalWorkload += demand * avgDist;
        }

        // 计算每个区域的平均工作量并设置边界值
        double avgWorkload = totalWorkload / p;
        this.demandLowerBound = (1 - r) * avgWorkload;
        this.demandUpperBound = (1 + r) * avgWorkload;
    }

    public Algo(Instance instance) {
        this.inst = instance;
        this.zones = new ArrayList[inst.k];
        this.r = 0.1;

        // 计算总需求并设置上下限
        double totalDemand = 0;
        for (int i = 0; i < inst.getN(); i++) {
            totalDemand += inst.getAreas()[i].getActiveness()[0];
        }
        this.demandUpperBound = (1 + r) * (totalDemand / inst.k);
        this.demandLowerBound = (1 - r) * (totalDemand / inst.k);
    }


    public void run(String filename) throws GRBException, IOException {
        long startTime = System.currentTimeMillis(); // 获取开始时间
        double beta = 0.4;
        double Best = Double.MAX_VALUE;
        int iter = 0;
        int MaxIter = 1;
        double alpha = 0;
        double delta = 0.01;
        ArrayList<Integer>[] BestZones = new ArrayList[inst.k];
        Random rand = new Random();
//        rand.setSeed(1097879689);
        while (iter < MaxIter) {
            centers = new ArrayList<>();
            int startId = rand.nextInt(inst.getN()); // 随机选择一个起始点
            // 将随机选择的起始点作为第一个大区域的中心
            centers.add(inst.getAreas()[startId]);
            inst.getAreas()[startId].setCenter(true);

            // 贪心搜索区域中心点
            while (centers.size() < inst.k) {
                //随机贪心挑选作为中心的区域
                // 记录所有不是中心点的区域到中心点区域的最小距离
                double[] minDistances = new double[inst.getN()];
                for (int i = 0; i < inst.getN(); i++) {
                    if (!inst.getAreas()[i].isCenter()) { // 只考虑未被选为中心的点
                        double minDist = Double.MAX_VALUE;
                        for (Area center : centers) {
                            if (inst.dist[center.getId()][i] < minDist) {
                                minDist = inst.dist[center.getId()][i];
                            }
                        }
                        minDistances[i] = minDist;
                    }
                }

                // 找到最大和最小距离
                double maxMinDist = -1;
                double minMinDist = Double.MAX_VALUE;
                for (Area area : inst.getAreas()) {
                    if (!area.isCenter()) { // 只考虑未被选为中心的点
                        double minDist = minDistances[area.getId()];
                        if (minDist > maxMinDist) {
                            maxMinDist = minDist;
                        }
                        if (minDist < minMinDist) {
                            minMinDist = minDist;
                        }
                    }
                }
                double Thre = maxMinDist - alpha * (maxMinDist - minMinDist);

                List<Area> candidates = new ArrayList<>();
                for (Area area : inst.getAreas()) {
                    if (!area.isCenter() && minDistances[area.getId()] >= Thre) {
                        candidates.add(area);
                    }
                }

                int nextId = rand.nextInt(candidates.size());
                centers.add(candidates.get(nextId));
                inst.getAreas()[nextId].setCenter(true);

            }

            boolean change = true;
            double cur_value = 0.0;
            while (change) {
                change = false;

                GRBEnv env = new GRBEnv(true);  // Create the env with manual start mode

// Set logging parameters BEFORE starting the environment
                env.set(GRB.IntParam.OutputFlag, 0);        // Suppress all output
                env.set(GRB.IntParam.LogToConsole, 0);      // Disable console logging
                env.set(GRB.StringParam.LogFile, "");       // Empty log file path
                env.set(GRB.IntParam.Seed, 42);
// Now start the environment
                env.start();

                // 禁用控制台输出


                GRBModel model = new GRBModel(env);
                model.set(GRB.IntParam.OutputFlag, 0); // 关闭标准输出
                model.set(GRB.IntParam.LogToConsole, 0); // 禁用控制台输出

                // 创建变量
                GRBVar[] x = new GRBVar[inst.getN() * inst.k];
                for (int i = 0; i < centers.size(); i++) {
                    for (int j = 0; j < inst.getN(); j++) {
                        x[i * inst.getN() + j] = model.addVar(0, 1, 0, GRB.BINARY, "x_" + centers.get(i).getId() + "_" + inst.getAreas()[j].getId());
                        if (centers.get(i).getId() == inst.getAreas()[j].getId()) {
                            //作为区域中心的变量直接设置为1
                            x[i * inst.getN() + j].set(GRB.DoubleAttr.LB, 1.0);
                            x[i * inst.getN() + j].set(GRB.DoubleAttr.UB, 1.0);
                        }
                    }
                }

                //每个点都必须属于一个区域
                for (int j = 0; j < inst.getN(); j++) {
                    GRBLinExpr expr = new GRBLinExpr();
                    for (int i = 0; i < centers.size(); i++) {
                        expr.addTerm(1.0, x[i * inst.getN() + j]);
                    }
                    model.addConstr(expr, GRB.EQUAL, 1.0, "c1_" + j);
                }


                // 添加约束：一个大区域内部所有小区域的一号活跃度之和加起来在[a,b]之间
                double a = (1 - r) * inst.average1;
                double b = (1 + r) * inst.average1;

                double a2 = (1 - r) * inst.average2;
                double b2 = (1 + r) * inst.average2;

                for (int i = 0; i < centers.size(); i++) {
                    GRBLinExpr expr = new GRBLinExpr();
                    GRBLinExpr expr2 = new GRBLinExpr();
                    for (int j = 0; j < inst.getN(); j++) {
                        expr.addTerm(inst.getAreas()[j].getActiveness()[0], x[i * inst.getN() + j]);
                        expr2.addTerm(inst.getAreas()[j].getActiveness()[1], x[i * inst.getN() + j]);
                    }
                    model.addConstr(expr, GRB.GREATER_EQUAL, a, "c2a_" + i);
                    model.addConstr(expr, GRB.LESS_EQUAL, b, "c2b_" + i);
                    model.addConstr(expr2, GRB.GREATER_EQUAL, a2, "c2a__" + i);
                    model.addConstr(expr2, GRB.LESS_EQUAL, b2, "c2b__" + i);
                }

                GRBLinExpr objExpr = new GRBLinExpr();
                for (int i = 0; i < centers.size(); i++)
                    for (int j = 0; j < inst.getN(); j++) {
                        objExpr.addTerm(inst.dist[centers.get(i).getId()][j], x[i * inst.getN() + j]);
                    }
//                model.write("facility_location.lp");
                model.setObjective(objExpr, GRB.MINIMIZE);
                model.optimize();
                double objVal = model.get(GRB.DoubleAttr.ObjVal);
                cur_value = objVal;


//                if (model.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL) {
//                    double objVal = model.get(GRB.DoubleAttr.ObjVal);
//                    cur_value = objVal;
//                    System.out.println("Optimal objective value: " + objVal);
//                } else {
//                    System.out.println("Model status is not optimal");
//                }

                // 检测每个变量是否处于最优解中

                for (int i = 0; i < centers.size(); i++) {
                    ArrayList<Integer> beatList = new ArrayList<>();
                    for (int j = 0; j < inst.getN(); j++) {
                        if (Math.abs(x[i * inst.getN() + j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                            beatList.add(j);
                        }
                    }
                    zones[i] = beatList;
                }


                for (int z = 0; z < zones.length; z++) {
                    ArrayList<Integer> zone = zones[z];
                    int oldCenter = centers.get(z).getId();
                    double minDist = Double.POSITIVE_INFINITY;
                    int newCenter = -1;
                    // 找到使得总距离最小的新中心
                    for (int i = 0; i < zone.size(); i++) {
                        int beat = zone.get(i);
                        double sumDist = 0.0;
                        for (int j = 0; j < zone.size(); j++) {
                            if (beat != zone.get(j)) {
                                sumDist += inst.dist[beat][zone.get(j)];
                            }
                        }
                        if (sumDist < minDist) {
                            minDist = sumDist;
                            newCenter = beat;
                        }
                    }
                    // 如果新中心与旧中心不同，则更新centers列表中的中心
                    if (newCenter != oldCenter) {
                        change = true;
                        centers.set(z, inst.getAreas()[newCenter]);
                    }
                }

                model.dispose();
                env.dispose();

            }

            FileWriter writer = new FileWriter("violated1.txt");
            BufferedWriter buffer = new BufferedWriter(writer);

            for (int io = 0; io < zones.length; io++) {
                buffer.write("center ID: " + centers.get(io).getId() + "\n");
                for (int jo = 0; jo < zones[io].size(); jo++) {
                    buffer.write(zones[io].get(jo) + " ");
                }
                buffer.newLine();
            }
            buffer.close();


            if (cur_value < Best) {
                boolean vio = true;
                //新建一个模型
                GRBEnv env = new GRBEnv(true);  // Create the env with manual start mode

// Set logging parameters BEFORE starting the environment
                env.set(GRB.IntParam.OutputFlag, 0);        // Suppress all output
                env.set(GRB.IntParam.LogToConsole, 0);      // Disable console logging
                env.set(GRB.StringParam.LogFile, "");       // Empty log file path
                env.set(GRB.IntParam.Seed, 42);
// Now start the environment
                env.start();

                // 禁用控制台输出

                GRBModel model = new GRBModel(env);
                model.set(GRB.IntParam.OutputFlag, 0); // 关闭标准输出
                model.set(GRB.IntParam.LogToConsole, 0); // 禁用控制台输出

                // 创建变量
                GRBVar[] x = new GRBVar[inst.getN() * inst.k];
                for (int i = 0; i < centers.size(); i++) {
                    for (int j = 0; j < inst.getN(); j++) {
                        x[i * inst.getN() + j] = model.addVar(0, 1, 0, GRB.BINARY, "x_" + centers.get(i).getId() + "_" + inst.getAreas()[j].getId());
                        if (centers.get(i).getId() == inst.getAreas()[j].getId()) {
                            //作为区域中心的变量直接设置为1
                            x[i * inst.getN() + j].set(GRB.DoubleAttr.LB, 1.0);
                            x[i * inst.getN() + j].set(GRB.DoubleAttr.UB, 1.0);
                        }
                    }
                }

                //每个点都必须属于一个区域
                for (int j = 0; j < inst.getN(); j++) {
                    GRBLinExpr expr = new GRBLinExpr();
                    for (int i = 0; i < centers.size(); i++) {
                        expr.addTerm(1.0, x[i * inst.getN() + j]);
                    }
                    model.addConstr(expr, GRB.EQUAL, 1.0, "c1_" + j);
                }


                // 添加约束：一个大区域内部所有小区域的一号活跃度之和加起来在[a,b]之间
                double a = (1 - r) * inst.average1;
                double b = (1 + r) * inst.average1;

                double a2 = (1 - r) * inst.average2;
                double b2 = (1 + r) * inst.average2;

                for (int i = 0; i < centers.size(); i++) {
                    GRBLinExpr expr = new GRBLinExpr();
                    GRBLinExpr expr2 = new GRBLinExpr();
                    for (int j = 0; j < inst.getN(); j++) {
                        expr.addTerm(inst.getAreas()[j].getActiveness()[0], x[i * inst.getN() + j]);
                        expr2.addTerm(inst.getAreas()[j].getActiveness()[1], x[i * inst.getN() + j]);
                    }
                    model.addConstr(expr, GRB.GREATER_EQUAL, a, "c2a_" + i);
                    model.addConstr(expr, GRB.LESS_EQUAL, b, "c2b_" + i);
                    model.addConstr(expr2, GRB.GREATER_EQUAL, a2, "c2a__" + i);
                    model.addConstr(expr2, GRB.LESS_EQUAL, b2, "c2b__" + i);
                }

                GRBLinExpr objExpr = new GRBLinExpr();
                for (int i = 0; i < centers.size(); i++)
                    for (int j = 0; j < inst.getN(); j++) {
                        objExpr.addTerm(inst.dist[centers.get(i).getId()][j], x[i * inst.getN() + j]);
                    }

                int maxConnectivityIterations = 100; // 最大连通性改善迭代次数
                int connectivityIteration = 0;
                while (vio && connectivityIteration < maxConnectivityIterations) {
                    connectivityIteration++;
                    vio = false;


                    for (int i = 0; i < zones.length; i++) {
                        ArrayList<ArrayList<Integer>> CA = new ArrayList<>();
                        ArrayList<Integer> Tm = new ArrayList<>();
                        Tm.addAll(zones[i]);
                        ArrayList<ArrayList<Integer>> CC = ConnectComponents(Tm);
                        if (CC.size() == 1)
                            continue;
//                        for (Integer kk : CC.get(0))
//                            System.out.print(kk + " ");

                        System.out.println();
                        System.out.println("连续部分有" + CC.size() + "个");
                        for (int j = 0; j < CC.size(); j++) {
                            if (!CC.get(j).contains(centers.get(i).getId())) {
                                CA.add(CC.get(j));
                                GRBLinExpr expr = new GRBLinExpr();
                                Set<Integer> neighbor = new LinkedHashSet<>();
                                for (Integer vv : CC.get(j)) {
                                    neighbor.addAll(inst.getAreas()[vv].getNeighbors());
                                }
                                neighbor.removeAll(CC.get(j));
                                for (Integer nei : neighbor) {
                                    expr.addTerm(1, x[i * inst.getN() + nei]);
                                }
//                                System.out.println("违反元素包含：");
                                for (Integer now : CC.get(j)) {
//                                    System.out.println(now);
                                    expr.addTerm(-1, x[i * inst.getN() + now]);
                                }
                                int m = 1 - CC.get(j).size();
                                model.addConstr(expr, GRB.GREATER_EQUAL, m, "violated_" + i);
                                vio = true;
//                                System.out.println("违反啦！" + "此时区域中心为" + centers.get(i).getId());

                            }
                        }
                    }
                    if (vio) {

                        System.out.println("违反啦");
//                        model.write("facility_location_violated.lp");
                        model.setObjective(objExpr, GRB.MINIMIZE);
                        model.optimize();
                        //加入新的约束之后，相应的区域会改变，但是centers不变
                        for (int im = 0; im < centers.size(); im++) {
                            ArrayList<Integer> beatList = new ArrayList<>();
                            for (int jm = 0; jm < inst.getN(); jm++) {
                                if (Math.abs(x[im * inst.getN() + jm].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                                    beatList.add(jm);
                                }
                            }
                            zones[im] = beatList;
                        }
                        cur_value = model.get(GRB.DoubleAttr.ObjVal);
                    } else {
                        if (cur_value < Best) {
                            Best = cur_value;
                            for (int z = 0; z < inst.k; z++) {
                                BestZones[z] = new ArrayList<>();
                                BestZones[z].addAll(zones[z]);
                            }
                        }

//                        System.out.println("没有违反约束");
//                        System.out.println("最终结果为" + Best);
                    }


                }

                model.dispose();
                env.dispose();

            }


            iter = iter + 1;
            System.out.println("第" + iter + "次迭代最好结果为" + Best);
            if (alpha < beta)
                alpha = alpha + delta;
            else
                alpha = 0;

            for (Area a : inst.getAreas()) {
                a.setCenter(false);
            }


        }

        long endTime = System.currentTimeMillis(); // 获取结束时间
        double timeSpentInSeconds = (endTime - startTime) / 1000.0;


        String outputFilePath = "./output/" + filename.replace(".dat", ".txt");
        FileWriter writer = new FileWriter(outputFilePath);
        BufferedWriter buffer = new BufferedWriter(writer);

        for (int io = 0; io < BestZones.length; io++) {
            buffer.write("center ID: " + centers.get(io).getId() + "\n");
            for (int jo = 0; jo < BestZones[io].size(); jo++) {
                buffer.write(BestZones[io].get(jo) + " ");
            }
            buffer.newLine();
        }

        String result = String.format("%.2f", Best);
        System.out.println(result);

        buffer.write("best objective: " + result + "\n");
        buffer.write("程序运行时间为：" + timeSpentInSeconds + "s" + "\n");
        buffer.close();

    }

    /**
     * 返回解决方案中的区域中心列表，不写入文件
     *
     * @return 区域中心ID列表
     */
    public ArrayList<Integer> getCorrectSolutionCenters() throws GRBException, IOException {
        return getCorrectSolutionCenters(false, null, null);
    }

    /**
     * 返回解决方案中的区域中心列表，不写入文件（支持改进模型）
     *
     * @param useImprovedModel 是否使用改进模型
     * @param shortestPathDist 最短路径距离矩阵（改进模型使用）
     * @param meanVector 均值向量（改进模型使用）
     * @return 区域中心ID列表
     */
    public ArrayList<Integer> getCorrectSolutionCenters(boolean useImprovedModel, double[][] shortestPathDist, double[] meanVector) throws GRBException, IOException {
        return getCorrectSolutionCenters(useImprovedModel, shortestPathDist, meanVector, false);
    }
    
    /**
     * 返回解决方案中的区域中心列表，不写入文件（支持改进模型和相对平衡性约束）
     *
     * @param useImprovedModel 是否使用改进模型
     * @param shortestPathDist 最短路径距离矩阵（改进模型使用）
     * @param meanVector 均值向量（改进模型使用）
     * @param useRelativeBalance 是否使用相对平衡性约束
     * @return 区域中心ID列表
     */
    public ArrayList<Integer> getCorrectSolutionCenters(boolean useImprovedModel, double[][] shortestPathDist, double[] meanVector, boolean useRelativeBalance) throws GRBException, IOException {
        // 设置改进模型参数
        this.useImprovedModel = useImprovedModel;
        this.shortestPathDist = shortestPathDist;
        this.meanVector = meanVector;
        this.useRelativeBalance = useRelativeBalance;

        return getCorrectSolutionCentersInternal();
    }

    /**
     * 内部方法：返回解决方案中的区域中心列表
     * 每次完整的算法流程总时间限制为60秒
     */
    private ArrayList<Integer> getCorrectSolutionCentersInternal() throws GRBException, IOException {
        long algorithmStartTime = System.currentTimeMillis(); // 记录算法开始时间
        long algorithmTimeLimit = 60 * 1000; // 每次完整算法流程总时间最多60秒（毫秒）
        
        long startTime = System.currentTimeMillis(); // 获取开始时间（用于其他用途）
        double beta = 0.4;
        double Best = Double.MAX_VALUE;
        ArrayList<Integer> bestCenters = new ArrayList<>();
        int iter = 0;
        int MaxIter = 1;
        double alpha = 0;
        double delta = 0.01;
        ArrayList<Integer>[] BestZones = new ArrayList[inst.k];
        Random rand = new Random();
        rand.setSeed(42);
        while (iter < MaxIter) {
            // 检查算法时间限制
            if (System.currentTimeMillis() - algorithmStartTime > algorithmTimeLimit) {
                System.out.println("【Algo】算法流程时间限制（60秒）已超过，返回当前最佳中心");
                // 如果有最佳中心，返回它；否则返回空列表
                if (!bestCenters.isEmpty()) {
                    return bestCenters;
                } else {
                    // 如果没有最佳中心，返回当前中心
                    ArrayList<Integer> currentCenters = new ArrayList<>();
                    for (Area center : centers) {
                        currentCenters.add(center.getId());
                    }
                    return currentCenters;
                }
            }
            centers = new ArrayList<>();
            int startId = rand.nextInt(inst.getN()); // 随机选择一个起始点
            // 将随机选择的起始点作为第一个大区域的中心
            centers.add(inst.getAreas()[startId]);
            inst.getAreas()[startId].setCenter(true);

            // 贪心搜索区域中心点
            while (centers.size() < inst.k) {
                //随机贪心挑选作为中心的区域
                // 记录所有不是中心点的区域到中心点区域的最小距离
                // 贪心选择中心点时：统一使用欧氏距离
                double[] minDistances = new double[inst.getN()];
                for (int i = 0; i < inst.getN(); i++) {
                    if (!inst.getAreas()[i].isCenter()) { // 只考虑未被选为中心的点
                        double minDist = Double.MAX_VALUE;
                        for (Area center : centers) {
                            // 统一使用欧氏距离
                            double dist = inst.dist[center.getId()][i];
                            if (dist < minDist) {
                                minDist = dist;
                            }
                        }
                        minDistances[i] = minDist;
                    }
                }

                // 找到最大和最小距离
                double maxMinDist = -1;
                double minMinDist = Double.MAX_VALUE;
                for (Area area : inst.getAreas()) {
                    if (!area.isCenter()) { // 只考虑未被选为中心的点
                        double minDist = minDistances[area.getId()];
                        if (minDist > maxMinDist) {
                            maxMinDist = minDist;
                        }
                        if (minDist < minMinDist) {
                            minMinDist = minDist;
                        }
                    }
                }
                double Thre = maxMinDist - alpha * (maxMinDist - minMinDist);

                List<Area> candidates = new ArrayList<>();
                for (Area area : inst.getAreas()) {
                    if (!area.isCenter() && minDistances[area.getId()] >= Thre) {
                        candidates.add(area);
                    }
                }

                int nextId = rand.nextInt(candidates.size());
                centers.add(candidates.get(nextId));
                inst.getAreas()[nextId].setCenter(true);

            }
            boolean change = true;
            double cur_value = 0.0;
            int maxCenterIterations = 100; // 最大中心点改善迭代次数
            int centerIteration = 0;
            while (change && centerIteration < maxCenterIterations) {
                // 检查算法时间限制
                if (System.currentTimeMillis() - algorithmStartTime > algorithmTimeLimit) {
                    System.out.println("【Algo】算法流程时间限制（60秒）已超过，返回当前最佳中心");
                    if (!bestCenters.isEmpty()) {
                        return bestCenters;
                    } else {
                        ArrayList<Integer> currentCenters = new ArrayList<>();
                        for (Area center : centers) {
                            currentCenters.add(center.getId());
                        }
                        return currentCenters;
                    }
                }
                
                centerIteration++;
                change = false;

                GRBEnv env = new GRBEnv(true);  // Create the env with manual start mode

// Set logging parameters BEFORE starting the environment
                env.set(GRB.IntParam.OutputFlag, 0);        // Suppress all output
                env.set(GRB.IntParam.LogToConsole, 0);      // Disable console logging
                env.set(GRB.StringParam.LogFile, "");       // Empty log file path
                env.set(GRB.IntParam.Seed, 42);
// Now start the environment
                env.start();

                GRBModel model = new GRBModel(env);
                model.set(GRB.IntParam.OutputFlag, 0); // 关闭标准输出
                model.set(GRB.IntParam.LogToConsole, 0); // 禁用控制台输出

                // 创建变量
                GRBVar[] x = new GRBVar[inst.getN() * inst.k];
                for (int i = 0; i < centers.size(); i++) {
                    for (int j = 0; j < inst.getN(); j++) {
                        x[i * inst.getN() + j] = model.addVar(0, 1, 0, GRB.BINARY, "x_" + centers.get(i).getId() + "_" + inst.getAreas()[j].getId());
                        if (centers.get(i).getId() == inst.getAreas()[j].getId()) {
                            //作为区域中心的变量直接设置为1
                            x[i * inst.getN() + j].set(GRB.DoubleAttr.LB, 1.0);
                            x[i * inst.getN() + j].set(GRB.DoubleAttr.UB, 1.0);
                        }
                    }
                }

                //每个点都必须属于一个区域
                for (int j = 0; j < inst.getN(); j++) {
                    GRBLinExpr expr = new GRBLinExpr();
                    for (int i = 0; i < centers.size(); i++) {
                        expr.addTerm(1.0, x[i * inst.getN() + j]);
                    }
                    model.addConstr(expr, GRB.EQUAL, 1.0, "c1_" + j);
                }


                // 添加约束：根据模型类型和约束类型选择不同的约束
                if (useRelativeBalance) {
                    // 相对平衡性约束：对于每个区域 j，需求比例应该在 [1/p*(1-α), 1/p*(1+α)] 范围内
                    // 等价于：(1-α)/p * sum_k d^T * x_k - d^T * x_j <= 0
                    //        d^T * x_j - (1+α)/p * sum_k d^T * x_k <= 0
                    int p = centers.size();
                    double coeff = (1.0 - r) / p; // (1-r)/p
                    double coeffUpper = (1.0 + r) / p; // (1+r)/p
                    
                    if (useImprovedModel && shortestPathDist != null) {
                        // 改进模型：使用工作量约束 w_j = r_j ⊙ x_j
                        // 为每个区域添加相对平衡性约束
                        for (int i = 0; i < centers.size(); i++) {
                            int centerId = centers.get(i).getId();
                            
                            // 下界约束：(1-α)/p * sum_k w^T * d_k - w^T * d_i <= 0
                            GRBLinExpr lowerExpr = new GRBLinExpr();
                            // 添加 (1-α)/p * sum_k w^T * d_k = coeff * sum_k sum_j d_j * r_jk * x_kj
                            for (int j = 0; j < inst.getN(); j++) {
                                for (int k = 0; k < centers.size(); k++) {
                                    int centerId_k = centers.get(k).getId();
                                    double r_jk = shortestPathDist[j][centerId_k];
                                    lowerExpr.addTerm(coeff * inst.getAreas()[j].getActiveness()[0] * r_jk, x[k * inst.getN() + j]);
                                }
                            }
                            
                            // 减去 w^T * d_i = sum_j d_j * r_ij * x_ij
                            for (int j = 0; j < inst.getN(); j++) {
                                double r_ij = shortestPathDist[j][centerId];
                                lowerExpr.addTerm(-inst.getAreas()[j].getActiveness()[0] * r_ij, x[i * inst.getN() + j]);
                            }
                            model.addConstr(lowerExpr, GRB.LESS_EQUAL, 0, "rel_balance_lower_" + i);
                            
                            // 上界约束：w^T * d_i - (1+α)/p * sum_k w^T * d_k <= 0
                            GRBLinExpr upperExpr = new GRBLinExpr();
                            // 添加 w^T * d_i = sum_j d_j * r_ij * x_ij
                            for (int j = 0; j < inst.getN(); j++) {
                                double r_ij = shortestPathDist[j][centerId];
                                upperExpr.addTerm(inst.getAreas()[j].getActiveness()[0] * r_ij, x[i * inst.getN() + j]);
                            }
                            
                            // 减去 (1+α)/p * sum_k w^T * d_k = -coeffUpper * sum_k sum_j d_j * r_jk * x_kj
                            for (int j = 0; j < inst.getN(); j++) {
                                for (int k = 0; k < centers.size(); k++) {
                                    int centerId_k = centers.get(k).getId();
                                    double r_jk = shortestPathDist[j][centerId_k];
                                    upperExpr.addTerm(-coeffUpper * inst.getAreas()[j].getActiveness()[0] * r_jk, x[k * inst.getN() + j]);
                                }
                            }
                            model.addConstr(upperExpr, GRB.LESS_EQUAL, 0, "rel_balance_upper_" + i);
                        }
                    } else {
                        // 原始模型：使用需求约束
                        // 为每个区域添加相对平衡性约束
                        for (int i = 0; i < centers.size(); i++) {
                            // 下界约束：(1-α)/p * sum_k d^T * x_k - d^T * x_i <= 0
                            GRBLinExpr lowerExpr = new GRBLinExpr();
                            // 添加 (1-α)/p * sum_k d^T * x_k = coeff * sum_k sum_j d_j * x_kj
                            for (int j = 0; j < inst.getN(); j++) {
                                for (int k = 0; k < centers.size(); k++) {
                                    lowerExpr.addTerm(coeff * inst.getAreas()[j].getActiveness()[0], x[k * inst.getN() + j]);
                                }
                            }
                            
                            // 减去 d^T * x_i = sum_j d_j * x_ij
                            for (int j = 0; j < inst.getN(); j++) {
                                lowerExpr.addTerm(-inst.getAreas()[j].getActiveness()[0], x[i * inst.getN() + j]);
                            }
                            model.addConstr(lowerExpr, GRB.LESS_EQUAL, 0, "rel_balance_lower_" + i);
                            
                            // 上界约束：d^T * x_i - (1+α)/p * sum_k d^T * x_k <= 0
                            GRBLinExpr upperExpr = new GRBLinExpr();
                            // 添加 d^T * x_i = sum_j d_j * x_ij
                            for (int j = 0; j < inst.getN(); j++) {
                                upperExpr.addTerm(inst.getAreas()[j].getActiveness()[0], x[i * inst.getN() + j]);
                            }
                            
                            // 减去 (1+α)/p * sum_k d^T * x_k = -coeffUpper * sum_k sum_j d_j * x_kj
                            for (int j = 0; j < inst.getN(); j++) {
                                for (int k = 0; k < centers.size(); k++) {
                                    upperExpr.addTerm(-coeffUpper * inst.getAreas()[j].getActiveness()[0], x[k * inst.getN() + j]);
                                }
                            }
                            model.addConstr(upperExpr, GRB.LESS_EQUAL, 0, "rel_balance_upper_" + i);
                        }
                    }
                } else {
                    // 绝对平衡性约束：使用上下界约束
                    if (useImprovedModel && shortestPathDist != null) {
                        // 改进模型：使用工作量约束 w_j = r_j ⊙ x_j
                        for (int i = 0; i < centers.size(); i++) {
                            int centerId = centers.get(i).getId();
                            // 创建 w_j 变量
                            // GRBVar[] w_j = new GRBVar[inst.getN()];
                            GRBLinExpr workloadExpr = new GRBLinExpr();
                            for (int j = 0; j < inst.getN(); j++) {
                                double r_ij = shortestPathDist[j][centerId];
                                // w_j[j] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "w_" + i + "_" + j);
                                // 约束：w_j[j] = r_ij * x[i*inst.getN() + j]
                                // GRBLinExpr wExpr = new GRBLinExpr();
                                // wExpr.addTerm(r_ij, x[i * inst.getN() + j]);
                                // model.addConstr(w_j[j], GRB.EQUAL, wExpr, "w_def_" + i + "_" + j);
                                // 工作量：d_j * w_j[j] = d_j * r_ij * x_ij
                                workloadExpr.addTerm(inst.getAreas()[j].getActiveness()[0] * r_ij, x[i * inst.getN() + j]);
                            }
                            model.addConstr(workloadExpr, GRB.GREATER_EQUAL, demandLowerBound, "c2a_" + i);
                            model.addConstr(workloadExpr, GRB.LESS_EQUAL, demandUpperBound, "c2b_" + i);
                        }
                    } else {
                        // 原始模型：使用需求约束
                        for (int i = 0; i < centers.size(); i++) {
                            GRBLinExpr expr = new GRBLinExpr();
                            for (int j = 0; j < inst.getN(); j++) {
                                expr.addTerm(inst.getAreas()[j].getActiveness()[0], x[i * inst.getN() + j]);
                            }
                            model.addConstr(expr, GRB.GREATER_EQUAL, demandLowerBound, "c2a_" + i);
                            model.addConstr(expr, GRB.LESS_EQUAL, demandUpperBound, "c2b_" + i);
                        }
                    }
                }

                // 设置目标函数：根据模型类型选择不同的目标函数
                GRBLinExpr objExpr = new GRBLinExpr();
                if (useImprovedModel && shortestPathDist != null) {
                    // 改进模型：min ∑_j ∑_i d_i * r_ij * x_ij（使用场景实例的实际需求值）
                    for (int i = 0; i < centers.size(); i++) {
                        int centerId = centers.get(i).getId();
                        for (int j = 0; j < inst.getN(); j++) {
                            double r_ij = shortestPathDist[j][centerId];
                            // 使用场景实例的实际需求值，而不是均值
                            double demand_j = inst.getAreas()[j].getActiveness()[0];
                            objExpr.addTerm(demand_j * r_ij, x[i * inst.getN() + j]);
                        }
                    }
                } else {
                    // 原始模型：min ∑_j ∑_i dist[i][j] * x_ij
                    for (int i = 0; i < centers.size(); i++)
                        for (int j = 0; j < inst.getN(); j++) {
                            objExpr.addTerm(inst.dist[centers.get(i).getId()][j], x[i * inst.getN() + j]);
                        }
                }
//                model.write("facility_location.lp");
                model.setObjective(objExpr, GRB.MINIMIZE);
                model.optimize();
                double objVal = model.get(GRB.DoubleAttr.ObjVal);
                cur_value = objVal;
//                if (model.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL) {
//                    double objVal = model.get(GRB.DoubleAttr.ObjVal);
//                    cur_value = objVal;
//                    System.out.println("Optimal objective value: " + objVal);
//                } else {
//                    System.out.println("Model status is not optimal");
//                }

                // 检测每个变量是否处于最优解中
                for (int i = 0; i < centers.size(); i++) {
                    ArrayList<Integer> beatList = new ArrayList<>();
                    for (int j = 0; j < inst.getN(); j++) {
                        if (Math.abs(x[i * inst.getN() + j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                            beatList.add(j);
                        }
                    }
                    zones[i] = beatList;
                }

                for (int z = 0; z < zones.length; z++) {
                    ArrayList<Integer> zone = zones[z];
                    int oldCenter = centers.get(z).getId();//获得区域中心集合里，第z个位置的中心的id，oldcenter代表真正的id
                    double minDist = Double.POSITIVE_INFINITY;
                    int newCenter = -1;
                    // 找到使得总距离最小的新中心
                    for (int i = 0; i < zone.size(); i++) {
                        int beat = zone.get(i);//beat就是基本单元的id
                        double sumDist = 0.0;
                        for (int j = 0; j < zone.size(); j++) {
                            if (beat != zone.get(j)) {
                                if (useImprovedModel && shortestPathDist != null && meanVector != null) {
                                    // 改进模型：使用最短路距离乘以需求 r_ij * d_j
                                    double r_ij = shortestPathDist[beat][zone.get(j)];
                                    double demand_j = meanVector[zone.get(j)];
                                    sumDist += r_ij * demand_j;
                                } else {
                                    // 原始模型：使用欧氏距离
                                    sumDist += inst.dist[beat][zone.get(j)];
                                }
                            }
                        }
                        if (sumDist < minDist) {
                            minDist = sumDist;
                            newCenter = beat; //newCenter一定会发生更新
                        }
                    }
                    // 如果新中心与旧中心不同，则更新centers列表中的中心
                    if (newCenter != oldCenter) {
                        change = true;
                        centers.set(z, inst.getAreas()[newCenter]);
                    }
                }

                model.dispose();
                env.dispose();
                
                // 如果达到最大中心点迭代次数仍未收敛，给出警告
                if (centerIteration >= maxCenterIterations && change) {
                    System.out.println("警告：达到最大中心点改善迭代次数（" + maxCenterIterations + "），可能未完全收敛");
                }

            }

            if (cur_value < Best) {
                boolean vio = true;
                //新建一个模型
                // Create the environment
                GRBEnv env = new GRBEnv(true);  // Create the env with manual start mode

// Set logging parameters BEFORE starting the environment
                env.set(GRB.IntParam.OutputFlag, 0);        // Suppress all output
                env.set(GRB.IntParam.LogToConsole, 0);      // Disable console logging
                env.set(GRB.StringParam.LogFile, "");       // Empty log file path
                env.set(GRB.IntParam.Seed, 42);
// Now start the environment
                env.start();


                // 禁用控制台输出

                GRBModel model = new GRBModel(env);
                model.set(GRB.IntParam.OutputFlag, 0); // 关闭标准输出
                model.set(GRB.IntParam.LogToConsole, 0); // 禁用控制台输出

                // 创建变量
                GRBVar[] x = new GRBVar[inst.getN() * inst.k];
                for (int i = 0; i < centers.size(); i++) {
                    for (int j = 0; j < inst.getN(); j++) {
                        x[i * inst.getN() + j] = model.addVar(0, 1, 0, GRB.BINARY, "x_" + centers.get(i).getId() + "_" + inst.getAreas()[j].getId());
                        if (centers.get(i).getId() == inst.getAreas()[j].getId()) {
                            //作为区域中心的变量直接设置为1
                            x[i * inst.getN() + j].set(GRB.DoubleAttr.LB, 1.0);
                            x[i * inst.getN() + j].set(GRB.DoubleAttr.UB, 1.0);
                        }
                    }
                }

                //每个点都必须属于一个区域
                for (int j = 0; j < inst.getN(); j++) {
                    GRBLinExpr expr = new GRBLinExpr();
                    for (int i = 0; i < centers.size(); i++) {
                        expr.addTerm(1.0, x[i * inst.getN() + j]);
                    }
                    model.addConstr(expr, GRB.EQUAL, 1.0, "c1_" + j);
                }


                // 添加约束：根据模型类型和约束类型选择不同的约束
                if (useRelativeBalance) {
                    // 相对平衡性约束：对于每个区域 j，需求比例应该在 [1/p*(1-α), 1/p*(1+α)] 范围内
                    // 等价于：(1-α)/p * sum_k d^T * x_k - d^T * x_j <= 0
                    //        d^T * x_j - (1+α)/p * sum_k d^T * x_k <= 0
                    int p = centers.size();
                    double coeff = (1.0 - r) / p; // (1-r)/p
                    double coeffUpper = (1.0 + r) / p; // (1+r)/p
                    
                    if (useImprovedModel && shortestPathDist != null) {
                        // 改进模型：使用工作量约束 w_j = r_j ⊙ x_j
                        // 为每个区域添加相对平衡性约束
                        for (int i = 0; i < centers.size(); i++) {
                            int centerId = centers.get(i).getId();
                            
                            // 下界约束：(1-α)/p * sum_k w^T * d_k - w^T * d_i <= 0
                            GRBLinExpr lowerExpr = new GRBLinExpr();
                            // 添加 (1-α)/p * sum_k w^T * d_k = coeff * sum_k sum_j d_j * r_jk * x_kj
                            for (int j = 0; j < inst.getN(); j++) {
                                for (int k = 0; k < centers.size(); k++) {
                                    int centerId_k = centers.get(k).getId();
                                    double r_jk = shortestPathDist[j][centerId_k];
                                    lowerExpr.addTerm(coeff * inst.getAreas()[j].getActiveness()[0] * r_jk, x[k * inst.getN() + j]);
                                }
                            }
                            
                            // 减去 w^T * d_i = sum_j d_j * r_ij * x_ij
                            for (int j = 0; j < inst.getN(); j++) {
                                double r_ij = shortestPathDist[j][centerId];
                                lowerExpr.addTerm(-inst.getAreas()[j].getActiveness()[0] * r_ij, x[i * inst.getN() + j]);
                            }
                            model.addConstr(lowerExpr, GRB.LESS_EQUAL, 0, "rel_balance_lower_" + i);
                            
                            // 上界约束：w^T * d_i - (1+α)/p * sum_k w^T * d_k <= 0
                            GRBLinExpr upperExpr = new GRBLinExpr();
                            // 添加 w^T * d_i = sum_j d_j * r_ij * x_ij
                            for (int j = 0; j < inst.getN(); j++) {
                                double r_ij = shortestPathDist[j][centerId];
                                upperExpr.addTerm(inst.getAreas()[j].getActiveness()[0] * r_ij, x[i * inst.getN() + j]);
                            }
                            
                            // 减去 (1+α)/p * sum_k w^T * d_k = -coeffUpper * sum_k sum_j d_j * r_jk * x_kj
                            for (int j = 0; j < inst.getN(); j++) {
                                for (int k = 0; k < centers.size(); k++) {
                                    int centerId_k = centers.get(k).getId();
                                    double r_jk = shortestPathDist[j][centerId_k];
                                    upperExpr.addTerm(-coeffUpper * inst.getAreas()[j].getActiveness()[0] * r_jk, x[k * inst.getN() + j]);
                                }
                            }
                            model.addConstr(upperExpr, GRB.LESS_EQUAL, 0, "rel_balance_upper_" + i);
                        }
                    } else {
                        // 原始模型：使用需求约束
                        // 为每个区域添加相对平衡性约束
                        for (int i = 0; i < centers.size(); i++) {
                            // 下界约束：(1-α)/p * sum_k d^T * x_k - d^T * x_i <= 0
                            GRBLinExpr lowerExpr = new GRBLinExpr();
                            // 添加 (1-α)/p * sum_k d^T * x_k = coeff * sum_k sum_j d_j * x_kj
                            for (int j = 0; j < inst.getN(); j++) {
                                for (int k = 0; k < centers.size(); k++) {
                                    lowerExpr.addTerm(coeff * inst.getAreas()[j].getActiveness()[0], x[k * inst.getN() + j]);
                                }
                            }
                            
                            // 减去 d^T * x_i = sum_j d_j * x_ij
                            for (int j = 0; j < inst.getN(); j++) {
                                lowerExpr.addTerm(-inst.getAreas()[j].getActiveness()[0], x[i * inst.getN() + j]);
                            }
                            model.addConstr(lowerExpr, GRB.LESS_EQUAL, 0, "rel_balance_lower_" + i);
                            
                            // 上界约束：d^T * x_i - (1+α)/p * sum_k d^T * x_k <= 0
                            GRBLinExpr upperExpr = new GRBLinExpr();
                            // 添加 d^T * x_i = sum_j d_j * x_ij
                            for (int j = 0; j < inst.getN(); j++) {
                                upperExpr.addTerm(inst.getAreas()[j].getActiveness()[0], x[i * inst.getN() + j]);
                            }
                            
                            // 减去 (1+α)/p * sum_k d^T * x_k = -coeffUpper * sum_k sum_j d_j * x_kj
                            for (int j = 0; j < inst.getN(); j++) {
                                for (int k = 0; k < centers.size(); k++) {
                                    upperExpr.addTerm(-coeffUpper * inst.getAreas()[j].getActiveness()[0], x[k * inst.getN() + j]);
                                }
                            }
                            model.addConstr(upperExpr, GRB.LESS_EQUAL, 0, "rel_balance_upper_" + i);
                        }
                    }
                } else {
                    // 绝对平衡性约束：使用上下界约束
                    if (useImprovedModel && shortestPathDist != null) {
                        // 改进模型：使用工作量约束 w_j = r_j ⊙ x_j
                        for (int i = 0; i < centers.size(); i++) {
                            int centerId = centers.get(i).getId();
                            // 创建 w_j 变量
                            // GRBVar[] w_j = new GRBVar[inst.getN()];
                            GRBLinExpr workloadExpr = new GRBLinExpr();
                            for (int j = 0; j < inst.getN(); j++) {
                                double r_ij = shortestPathDist[j][centerId];
                                // w_j[j] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "w_" + i + "_" + j);
                                // 约束：w_j[j] = r_ij * x[i*inst.getN() + j]
                                // GRBLinExpr wExpr = new GRBLinExpr();
                                // wExpr.addTerm(r_ij, x[i * inst.getN() + j]);
                                // model.addConstr(w_j[j], GRB.EQUAL, wExpr, "w_def_" + i + "_" + j);
                                // 工作量：d_j * w_j[j] = d_j * r_ij * x_ij
                                workloadExpr.addTerm(inst.getAreas()[j].getActiveness()[0] * r_ij, x[i * inst.getN() + j]);
                            }
                            model.addConstr(workloadExpr, GRB.GREATER_EQUAL, demandLowerBound, "c2a_" + i);
                            model.addConstr(workloadExpr, GRB.LESS_EQUAL, demandUpperBound, "c2b_" + i);
                        }
                    } else {
                        // 原始模型：使用需求约束
                        for (int i = 0; i < centers.size(); i++) {
                            GRBLinExpr expr = new GRBLinExpr();
                            for (int j = 0; j < inst.getN(); j++) {
                                expr.addTerm(inst.getAreas()[j].getActiveness()[0], x[i * inst.getN() + j]);
                            }
                            model.addConstr(expr, GRB.GREATER_EQUAL, demandLowerBound, "c2a_" + i);
                            model.addConstr(expr, GRB.LESS_EQUAL, demandUpperBound, "c2b_" + i);
                        }
                    }
                }

                // 设置目标函数：根据模型类型选择不同的目标函数
                GRBLinExpr objExpr = new GRBLinExpr();
                if (useImprovedModel && shortestPathDist != null) {
                    // 改进模型：min ∑_j ∑_i d_i * r_ij * x_ij（使用实例的实际需求值）
                    for (int i = 0; i < centers.size(); i++) {
                        int centerId = centers.get(i).getId();
                        for (int j = 0; j < inst.getN(); j++) {
                            double r_ij = shortestPathDist[j][centerId];
                            // 使用实例的实际需求值，而不是均值
                            double demand_j = inst.getAreas()[j].getActiveness()[0];
                            objExpr.addTerm(demand_j * r_ij, x[i * inst.getN() + j]);
                        }
                    }
                } else {
                    // 原始模型：min ∑_j ∑_i dist[i][j] * x_ij
                    for (int i = 0; i < centers.size(); i++)
                        for (int j = 0; j < inst.getN(); j++) {
                            objExpr.addTerm(inst.dist[centers.get(i).getId()][j], x[i * inst.getN() + j]);
                        }
                }

                int maxConnectivityIterations = 100; // 最大连通性改善迭代次数
                int connectivityIteration = 0;
                while (vio && connectivityIteration < maxConnectivityIterations) {
                    // 检查算法时间限制
                    if (System.currentTimeMillis() - algorithmStartTime > algorithmTimeLimit) {
                        System.out.println("【Algo】算法流程时间限制（60秒）已超过，返回当前最佳中心");
                        if (!bestCenters.isEmpty()) {
                            model.dispose();
                            env.dispose();
                            return bestCenters;
                        } else {
                            ArrayList<Integer> currentCenters = new ArrayList<>();
                            for (Area center : centers) {
                                currentCenters.add(center.getId());
                            }
                            model.dispose();
                            env.dispose();
                            return currentCenters;
                        }
                    }
                    
                    connectivityIteration++;
                    vio = false;


                    for (int i = 0; i < zones.length; i++) {
                        ArrayList<ArrayList<Integer>> CA = new ArrayList<>();
                        ArrayList<Integer> Tm = new ArrayList<>();
                        Tm.addAll(zones[i]);
                        ArrayList<ArrayList<Integer>> CC = ConnectComponents(Tm);
                        if (CC.size() == 1)
                            continue;
//                        System.out.println();
//                        System.out.println("连续部分有" + CC.size() + "个");
                        for (int j = 0; j < CC.size(); j++) {
                            if (!CC.get(j).contains(centers.get(i).getId())) {
                                CA.add(CC.get(j));
                                GRBLinExpr expr = new GRBLinExpr();
                                Set<Integer> neighbor = new LinkedHashSet<>();
                                for (Integer vv : CC.get(j)) {
                                    neighbor.addAll(inst.getAreas()[vv].getNeighbors());
                                }
                                neighbor.removeAll(CC.get(j));
                                for (Integer nei : neighbor) {
                                    expr.addTerm(1, x[i * inst.getN() + nei]);
                                }
//                                System.out.println("违反元素包含：");
                                for (Integer now : CC.get(j)) {
//                                    System.out.println(now);
                                    expr.addTerm(-1, x[i * inst.getN() + now]);
                                }
                                int m = 1 - CC.get(j).size();
                                model.addConstr(expr, GRB.GREATER_EQUAL, m, "violated_" + i);
                                vio = true;
//                                System.out.println("违反啦！" + "此时区域中心为" + centers.get(i).getId());

                            }
                        }
                    }
                    if (vio) {
//                        System.out.println("-------违反连通性啦----------");
//                        model.write("facility_location_violated.lp");
                        model.setObjective(objExpr, GRB.MINIMIZE);
                        model.optimize();
                        //加入新的约束之后，相应的区域会改变，但是centers不变
                        for (int im = 0; im < centers.size(); im++) {
                            ArrayList<Integer> beatList = new ArrayList<>();
                            for (int jm = 0; jm < inst.getN(); jm++) {
                                if (Math.abs(x[im * inst.getN() + jm].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                                    beatList.add(jm);
                                }
                            }
                            zones[im] = beatList;
                        }
                        cur_value = model.get(GRB.DoubleAttr.ObjVal);
                    } else {
                        if (cur_value < Best) {
                            Best = cur_value;
                            for (int z = 0; z < inst.k; z++) {
                                BestZones[z] = new ArrayList<>();
                                BestZones[z].addAll(zones[z]);
                            }
                            // 清空并更新最佳中心列表
                            bestCenters.clear();
                            for (int kk = 0; kk < centers.size(); kk++) {
                                bestCenters.add(centers.get(kk).getId());
                            }
                        }

                    }
                }
                
                // 如果达到最大连通性迭代次数仍未满足连通性，给出警告
                if (connectivityIteration >= maxConnectivityIterations && vio) {
                    System.out.println("警告：达到最大连通性迭代次数（" + maxConnectivityIterations + "），可能未完全满足连通性约束");
                }

                model.dispose();
                env.dispose();
            }
            
            // 最终检查算法时间限制
            if (System.currentTimeMillis() - algorithmStartTime > algorithmTimeLimit) {
                System.out.println("【Algo】算法流程时间限制（60秒）已超过，返回当前最佳中心");
                if (!bestCenters.isEmpty()) {
                    return bestCenters;
                } else {
                    ArrayList<Integer> currentCenters = new ArrayList<>();
                    for (Area center : centers) {
                        currentCenters.add(center.getId());
                    }
                    return currentCenters;
                }
            }
            
            iter = iter + 1;
//            System.out.println("第" + iter + "次迭代最好结果为" + Best);
//            if (alpha < beta)
//                alpha = alpha + delta;
//            else
//                alpha = 0;
//            for (Area a : inst.getAreas()) {
//                a.setCenter(false);
//            }

        }

        long endTime = System.currentTimeMillis(); // 获取结束时间
        double timeSpentInSeconds = (endTime - startTime) / 1000.0;
        return bestCenters;
    }


    public ArrayList<ArrayList<Integer>> ViolatedCC(ArrayList<Integer> zone) {
        ArrayList<ArrayList<Integer>> CA = new ArrayList<>();
        ArrayList<Integer> Tm = new ArrayList<>();
        Tm.addAll(zone);
        ArrayList<ArrayList<Integer>> CC = ConnectComponents(Tm);
        for (int j = 0; j < CC.size(); j++) {

        }

        return CA;
    }

    public ArrayList<ArrayList<Integer>> ConnectComponents(ArrayList<Integer> zone) {
        //  zone本身就是区域的编号
        ArrayList<ArrayList<Integer>> CC = new ArrayList<>();
        boolean[] visited = new boolean[zone.size()];
        for (int i = 0; i < zone.size(); i++) {
            visited[i] = false;
        }
        int cnt = 0;
        for (int i = 0; i < zone.size(); i++) {
            if (!visited[i]) {
                ArrayList<Integer> CompConex = new ArrayList<>();
                Queue<Integer> L = new LinkedList();
                L.add(zone.get(i));
                cnt = cnt + 1;
                visited[i] = true;
                while (!L.isEmpty()) {
                    int cur = L.remove();
                    CompConex.add(cur);
                    for (int m = 0; m < inst.getAreas()[cur].getNeighbors().size(); m++) {
                        int neighbor = inst.getAreas()[cur].getNeighbors().get(m);
                        if (zone.contains(neighbor) && !visited[zone.indexOf(neighbor)]) {
                            L.add(neighbor);
                            visited[zone.indexOf(neighbor)] = true;
                        }
                    }
                }
                CC.add(CompConex);
            }
        }
//        for (int i = 0; i < CC.size(); i++) {
//            System.out.println("第" + (i + 1) + "个连续部分");
//            for (int j = 0; j < CC.get(i).size(); j++) {
//                System.out.println(CC.get(i).get(j));
//            }
//        }

        return CC;
    }

}


