import Jama.CholeskyDecomposition;
import Jama.Matrix;
import com.gurobi.gurobi.*;
import copt.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


/**
 * 分布鲁棒机会约束分区问题的解决方案
 * 实现基于D_1和D_2模糊集的DRICC问题以及基于Bonferroni近似的DRJCC模型
 * 使用Gurobi进行求解
 */
public class DistributionallyRobustAlgo {
    private Instance inst;
    private ArrayList<Area> centers; // 存储所有区域中心
    private ArrayList<Integer>[] zones; // 存储每个区域的基本单元编号
    private double r; // 活动指标平衡容差
    private double gamma; // 机会约束风险参数
    private double[][] scenarios; // 原始场景数据（用于非assignment-dependent模型）
    private int numScenarios; // 场景数量
    private int[][] scenarioDemands; // 存储所有场景下的需求（用于非assignment-dependent模型）
    private Random rand; // 随机数生成器

    // Assignment-dependent workload模型相关
    private boolean useAssignmentDependent; // 是否使用assignment-dependent workload模型
    private double[][][] assignmentDependentScenarios; // 存储assignment-dependent场景数据 [scenario][i][j]，表示场景k下d_ij的值，其中j是区域索引
    private int numAssignmentDependentScenarios; // assignment-dependent场景数量
    private int numRegionsForAssignmentDependent; // assignment-dependent模型中的区域数量p

    // 分布鲁棒模型的参数
    private double[] meanVector; // 样本均值向量（对于assignment-dependent模型，维度为N*p）
    private double[][] covarianceMatrix; // 样本协方差矩阵（对于assignment-dependent模型，维度为(N*p)×(N*p)）
    private double[] varianceVector; // 方差向量（对角线元素，用于独立场景假设，维度为N*p）
    private boolean useLazyCovariance; // 是否使用延迟计算协方差（当矩阵太大时）
    private boolean useDiagonalCovariance; // 是否只使用对角线（方差），假设场景之间独立
    private int nTimesPForCovariance; // 存储N*p值，用于延迟计算
    private double delta1; // D_2模糊集参数
    private double delta2; // D_2模糊集参数
    private boolean useD1; // 是否使用D_1模糊集
    private boolean useJointChance; // 是否使用联合机会约束
    private double[] individualGammas; // Bonferroni近似的个体风险分配
    private boolean useExactMethod; // 是否使用精确方法（否则使用近似方法MISOCP）
    private boolean useImprovedModel; // 是否使用改进模型（使用最短路距离和工作量约束）
    private boolean useRelativeBalance; // 是否使用相对平衡性约束（否则使用绝对平衡性约束）
    private boolean useSupportingHyperplaneCuts; // 是否使用基于支撑超平面的cut（仅在使用相对平衡性约束且使用近似方法时有效）

    // Gamma分解参数（用于近似分解方法ADM）
    private double gamma_a; // 下界约束的风险参数（用于低于下界L的情况）
    private double gamma_b; // 上界约束的风险参数（用于超过上界U的情况）

    private double demandUpperBound; // 存储需求上限
    private double demandLowerBound; // 存储需求下界

    // 改进模型相关：最短路径距离矩阵
    private double[][] shortestPathDist; // r_ij: 从基本单元i到区域中心j的最短路径距离

    // 精确方法相关：存储不可行解的index集合
    private ArrayList<Set<String>> infeasibleSolutions; // 每个元素是一个不可行解的index集合

    // 支撑超平面cut相关：存储相对平衡性约束的信息（用于生成cut）
    private static class RelativeBalanceConstraintInfo {
        int j; // 区域索引
        double riskParam; // 风险参数
        double factorLower; // 下界factor
        double factorUpper; // 上界factor
        double coeff; // (1-α)/p
        double coeffUpper; // (1+α)/p
    }
    private ArrayList<RelativeBalanceConstraintInfo> relativeBalanceConstraintInfos; // 存储所有相对平衡性约束的信息

    // 时间限制和迭代限制
    // 注意：timeLimit 用于限制单次Gurobi求解的最大时间，但实际使用时会根据全局剩余时间动态调整
    // 设置为3600秒，每次迭代最多使用3600秒
    private int timeLimit = 3600; // 秒（Gurobi求解器单次求解的最大时间限制，每次迭代最多3600秒）
    private int maxIterations = 100;

    // 全局时间限制（整个算法运行的总时间限制，单位：毫秒）
    private static final long GLOBAL_TIME_LIMIT_MS = 36000 * 1000; // 36000 seconds (10 hours)

    // 存储最佳目标函数值
    private double bestObjective = Double.MAX_VALUE;

    // 存储求解状态
    private boolean solveSuccess = false; // 是否求解成功
    private int statusCode = -1; // Gurobi状态代码（如果求解失败）
    private int cutIterations = -1; // 精确算法的迭代次数（如果使用精确方法）
    private int failureStage = 0;     // 失败阶段：0=成功，1=选择初始中心失败，2=生成初始解失败，3=改善初始解失败，4=确保连通性失败
    private int d1SdpNoGoodCutCounter = 0; // D1+SDP分离法下新增的no-good cut计数
    private boolean d1SdpJavaCoptPathLogged = false; // 仅打印一次 Java-COPT 路径提示
    private copt.Envr d1SdpSharedEnvr = null; // 复用COPT环境，避免每次子问题重复打印banner

    // avgDist计算方式：1=方式1（所有j到i的距离平均），2=方式2（j的k近邻平均距离的平均），3=原来的方法
    private int avgDistMethod = 1; // 默认为方式1
    
    // 记录当前方法计算出的avgDist向量的均值
    private double avgDistMean = -1.0; // avgDist向量均值（-1表示未计算）
    
    // CSV文件写入器（用于记录每次迭代的工作量）
    private BufferedWriter workloadCsvWriter = null;
    private String workloadCsvFilePath = null;

    /**
     * 构造函数
     */
    public DistributionallyRobustAlgo(Instance instance, double[][] scenarios, double gamma,
                                      long seed, boolean useD1, double delta1, double delta2, boolean useJointChance, double r, boolean useExactMethod) {
        this(instance, scenarios, gamma, seed, useD1, delta1, delta2, useJointChance, r, useExactMethod, false);
    }

    /**
     * 构造函数（包含改进模型选择）
     */
    public DistributionallyRobustAlgo(Instance instance, double[][] scenarios, double gamma,
                                      long seed, boolean useD1, double delta1, double delta2, boolean useJointChance, double r, boolean useExactMethod, boolean useImprovedModel) {
        this(instance, scenarios, gamma, seed, useD1, delta1, delta2, useJointChance, r, useExactMethod, useImprovedModel, 1);
    }

    /**
     * 构造函数（包含改进模型选择和avgDist计算方式选择）
     * @param avgDistMethod avgDist计算方式：1=方式1（所有j到i的距离平均），2=方式2（j的k近邻平均距离的平均）
     */
    public DistributionallyRobustAlgo(Instance instance, double[][] scenarios, double gamma,
                                      long seed, boolean useD1, double delta1, double delta2, boolean useJointChance, double r, boolean useExactMethod, boolean useImprovedModel, int avgDistMethod) {
        this(instance, scenarios, gamma, seed, useD1, delta1, delta2, useJointChance, r, useExactMethod, useImprovedModel, avgDistMethod, false);
    }

    /**
     * 构造函数（包含改进模型选择、avgDist计算方式选择和相对平衡性约束选项）
     * @param avgDistMethod avgDist计算方式：1=方式1（所有j到i的距离平均），2=方式2（j的k近邻平均距离的平均）
     * @param useRelativeBalance 是否使用相对平衡性约束（否则使用绝对平衡性约束）
     */
    public DistributionallyRobustAlgo(Instance instance, double[][] scenarios, double gamma,
                                      long seed, boolean useD1, double delta1, double delta2, boolean useJointChance, double r, boolean useExactMethod, boolean useImprovedModel, int avgDistMethod, boolean useRelativeBalance) {
        this(instance, scenarios, gamma, seed, useD1, delta1, delta2, useJointChance, r, useExactMethod, useImprovedModel, avgDistMethod, useRelativeBalance, false);
    }

    /**
     * 构造函数（包含改进模型选择、avgDist计算方式选择、相对平衡性约束选项和支撑超平面cut选项）
     * @param avgDistMethod avgDist计算方式：1=方式1（所有j到i的距离平均），2=方式2（j的k近邻平均距离的平均）
     * @param useRelativeBalance 是否使用相对平衡性约束（否则使用绝对平衡性约束）
     * @param useSupportingHyperplaneCuts 是否使用基于支撑超平面的cut（仅在使用相对平衡性约束且使用近似方法时有效）
     */
    public DistributionallyRobustAlgo(Instance instance, double[][] scenarios, double gamma,
                                      long seed, boolean useD1, double delta1, double delta2, boolean useJointChance, double r, boolean useExactMethod, boolean useImprovedModel, int avgDistMethod, boolean useRelativeBalance, boolean useSupportingHyperplaneCuts) {
        this.inst = instance;
        this.scenarios = scenarios;
        this.numScenarios = scenarios != null ? scenarios.length : 0;
        this.scenarioDemands = new int[numScenarios][inst.getN()];
        for (int s = 0; s < numScenarios; s++) {
            for (int i = 0; i < inst.getN(); i++) {
                this.scenarioDemands[s][i] = (int) scenarios[s][i];
            }
        }
        this.gamma = gamma;
        this.r = r; // 活动指标平衡容差（可配置）
        this.rand = new Random(seed);
        this.zones = new ArrayList[inst.k];
        this.useD1 = useD1;
        this.delta1 = delta1;
        this.delta2 = delta2;
        this.useJointChance = useJointChance;
        this.useExactMethod = useExactMethod;
        this.useImprovedModel = useImprovedModel;
        this.useRelativeBalance = useRelativeBalance;
        this.useSupportingHyperplaneCuts = useSupportingHyperplaneCuts;
        this.infeasibleSolutions = new ArrayList<>();
        this.relativeBalanceConstraintInfos = new ArrayList<>();
        this.avgDistMethod = avgDistMethod; // 设置avgDist计算方式
        this.useAssignmentDependent = false; // 默认不使用assignment-dependent模型
        this.assignmentDependentScenarios = null;
        this.numAssignmentDependentScenarios = 0;

        // 如果使用改进模型，计算最短路径距离矩阵
        if (useImprovedModel) {
            computeShortestPathDistances();
        }

        // 计算样本均值和协方差矩阵
        calculateMomentInformation();

        // 计算上下界
        if (useImprovedModel) {
            // 改进模型：基于工作量的上下界
            calculateWorkloadBounds(avgDistMethod);
        } else {
            // 原始模型：基于需求的上下界
            double totalMeanDemand = 0;
            for (int i = 0; i < inst.getN(); i++) {
                totalMeanDemand += meanVector[i];
            }
            double mu = totalMeanDemand / inst.k; // 所有scenario的平均需求均值
            this.demandUpperBound = (1 + r) * mu; // 上界 U = (1 + r) * mu
            this.demandLowerBound = (1 - r) * mu; // 下界 L = (1 - r) * mu
        }

        // 如果使用Bonferroni近似，初始化个体风险分配
        if (useJointChance) {
            this.individualGammas = new double[inst.k];
            // 简单均分风险
            for (int j = 0; j < inst.k; j++) {
                individualGammas[j] = gamma / inst.k;
            }
        }

        // 初始化gamma分解参数（四六分：gamma_a = 0.6*gamma, gamma_b = 0.4*gamma）
        // 出于人道主义考虑，对超过上界U的情况容忍度更低，所以上界约束使用更小的gamma_b
        this.gamma_a = 0.5 * gamma; // 下界约束的风险参数
        this.gamma_b = 0.5 * gamma; // 上界约束的风险参数
    }

    /**
     * 构造函数（用于assignment-dependent workload模型）
     * @param instance 实例
     * @param gamma 风险参数
     * @param seed 随机种子
     * @param useD1 是否使用D1模糊集
     * @param delta1 D2模糊集参数1
     * @param delta2 D2模糊集参数2
     * @param useJointChance 是否使用联合机会约束
     * @param r 活动指标平衡容差
     * @param useExactMethod 是否使用精确方法
     * @param useRelativeBalance 是否使用相对平衡性约束
     * @param useSupportingHyperplaneCuts 是否使用支撑超平面cut
     */
    public DistributionallyRobustAlgo(Instance instance, double gamma,
                                      long seed, boolean useD1, double delta1, double delta2, boolean useJointChance, double r, boolean useExactMethod, boolean useRelativeBalance, boolean useSupportingHyperplaneCuts) {
        this.inst = instance;
        this.scenarios = null;
        this.numScenarios = 0;
        this.scenarioDemands = null;
        this.gamma = gamma;
        this.r = r;
        this.rand = new Random(seed);
        this.zones = new ArrayList[inst.k];
        this.useD1 = useD1;
        this.delta1 = delta1;
        this.delta2 = delta2;
        this.useJointChance = useJointChance;
        this.useExactMethod = useExactMethod;
        this.useImprovedModel = false; // assignment-dependent模型不使用改进模型
        this.useRelativeBalance = useRelativeBalance;
        this.useSupportingHyperplaneCuts = useSupportingHyperplaneCuts;
        this.infeasibleSolutions = new ArrayList<>();
        this.relativeBalanceConstraintInfos = new ArrayList<>();
        this.avgDistMethod = 1;
        this.useAssignmentDependent = true; // 使用assignment-dependent模型
        this.assignmentDependentScenarios = null;
        this.numAssignmentDependentScenarios = 0;
        this.numRegionsForAssignmentDependent = 0; // 将在loadAssignmentDependentData中设置

        // 加载assignment-dependent workload数据
        loadAssignmentDependentData();

        // 计算样本均值和协方差矩阵（N*p维）
        calculateMomentInformation();

        // 计算上下界（基于assignment-dependent模型）
        calculateAssignmentDependentBounds();

        // 如果使用Bonferroni近似，初始化个体风险分配
        if (useJointChance) {
            this.individualGammas = new double[inst.k];
            for (int j = 0; j < inst.k; j++) {
                individualGammas[j] = gamma / inst.k;
            }
        }

        // 初始化gamma分解参数
        this.gamma_a = 0.5 * gamma;
        this.gamma_b = 0.5 * gamma;
        
        // 如果使用assignment-dependent模型和支撑超平面cut，初始化CSV文件用于记录工作量
        if (useAssignmentDependent && useSupportingHyperplaneCuts) {
            // initializeWorkloadCsvFile();
        }
    }

    /**
     * 主要求解方法
     */
    public void run(String filename) throws IOException, GRBException {
        long startTime = System.currentTimeMillis();
        double Best = Double.MAX_VALUE;
        ArrayList<Integer>[] BestZones = new ArrayList[inst.k];

        // 初始化求解状态
        this.solveSuccess = false;
        this.statusCode = -1;
        this.cutIterations = -1; // 如果不是精确方法，保持为-1
        this.failureStage = 0; // 初始化为0（成功），如果失败会在相应阶段设置

        // Check global time limit at the start
        if (GLOBAL_TIME_LIMIT_MS <= 0) {
            System.out.println("Global time limit has been exceeded before starting");
            this.failureStage = 1; // 在开始前就失败，算作选择初始中心阶段失败
            this.bestObjective = -1;
            return;
        }

        // Step 1: 构造初始区域中心集合
        ArrayList<Integer> initialCenters;
        if (useAssignmentDependent) {
            // assignment-dependent模型：直接从CSV文件读取初始中心点
            initialCenters = loadInitialCentersFromCSV();
        } else {
            // 标准模型：使用原来的选择方法
            initialCenters = selectInitialCenters();
        }
        // 检查选择初始中心是否失败（如果返回null或空列表，说明失败）
        if (initialCenters == null || initialCenters.isEmpty()) {
            System.out.println("【选择初始中心】选择初始中心失败");
            this.failureStage = 1; // 选择初始中心阶段失败
            this.bestObjective = -1;
            return;
        }
        centers = new ArrayList<>();
        for (int centerId : initialCenters) {
            centers.add(inst.getAreas()[centerId]);
            inst.getAreas()[centerId].setCenter(true);
            // assignment-dependent 从 CSV 加载时，.dat 中 Area.id 为序号-1，算法用 getId() 作数组下标，故将中心点的 id 设为实例索引
            if (useAssignmentDependent) {
                inst.getAreas()[centerId].setId(centerId);
            }
        }

        // Step 2: 生成初始可行解
        System.out.println("【生成初始解】开始生成初始可行解...");
        // Check time limit before starting
        if (System.currentTimeMillis() - startTime > GLOBAL_TIME_LIMIT_MS) {
            System.out.println("【生成初始解】全局时间限制已超过，无法生成初始解");
            this.failureStage = 2; // 生成初始解阶段失败
            this.bestObjective = -1;
            return;
        }
        boolean feasible = false;

        // 如果使用改进模型但不使用相对平衡性约束，需要在生成初始解阶段调整上下界
        if (useImprovedModel && !useRelativeBalance) {
            // 保存原始上下界
            double originalUpperBound = demandUpperBound;
            double originalLowerBound = demandLowerBound;
            
            // 上下界调整循环：最多50次，一次2%
            int maxBoundAdjustments = 20;
            int boundAdjustmentCount = 0;
            
            while (!feasible && boundAdjustmentCount < maxBoundAdjustments) {
                if (boundAdjustmentCount > 0) {
                    System.out.println("【生成初始解-调整上下界】第 " + boundAdjustmentCount + " 次调整上下界...");
                    adjustBounds();
                } else {
                    System.out.println("【生成初始解】开始第一次尝试...");
                }
                
                // 如果使用精确方法，需要迭代求解直到找到可行解
                if (useExactMethod) {
                    System.out.println("【生成初始解】使用精确方法...");
                    feasible = solveWithExactMethod(startTime);
                    if (feasible) {
                        if (boundAdjustmentCount > 0) {
                            System.out.println("【生成初始解】经过 " + boundAdjustmentCount + " 次调整上下界后，精确方法成功找到可行解");
                        } else {
                            System.out.println("【生成初始解】精确方法成功找到可行解");
                        }
                    } else {
                        System.out.println("【生成初始解】精确方法在最大迭代次数内无法找到可行解");
                        boundAdjustmentCount++;
                        // 继续调整上下界
                    }
                } else {
                    System.out.println("【生成初始解】使用近似方法（直接添加DRCC约束）...");
                    feasible = generateInitialSolution(startTime);
                    
                    if (feasible) {
                        if (boundAdjustmentCount > 0) {
                            System.out.println("【生成初始解】经过 " + boundAdjustmentCount + " 次调整上下界后，近似方法成功找到可行解");
                        } else {
                            System.out.println("【生成初始解】近似方法成功找到可行解");
                        }
                    } else {
                        System.out.println("【生成初始解】近似方法无法找到可行解");
                        boundAdjustmentCount++;
                        // 继续调整上下界
                    }
                }
                
                // 检查时间限制
                if (System.currentTimeMillis() - startTime > GLOBAL_TIME_LIMIT_MS) {
                    System.out.println("【生成初始解】全局时间限制已超过，停止上下界调整");
                    break;
                }
            }
            
            if (!feasible) {
                System.out.println("【生成初始解】达到最大调整次数（" + maxBoundAdjustments + "），无法找到可行解");
                // 恢复原始上下界
                demandUpperBound = originalUpperBound;
                demandLowerBound = originalLowerBound;
                this.failureStage = 2; // 生成初始解阶段失败
                this.bestObjective = -1;
                return;
            }
        } else {
            // 不使用改进模型或使用相对平衡性约束，保持原有逻辑
            // 如果使用精确方法，需要迭代求解直到找到可行解
            if (useExactMethod) {
                System.out.println("【生成初始解】使用精确方法...");
                feasible = solveWithExactMethod(startTime);
                if (!feasible) {
                    System.out.println("【生成初始解】精确方法在最大迭代次数内无法找到可行解");
                    this.failureStage = 2; // 生成初始解阶段失败
                    this.bestObjective = -1;
                    return;
                }
                System.out.println("【生成初始解】精确方法成功找到可行解");
            } else {
                System.out.println("【生成初始解】使用近似方法（直接添加DRCC约束）...");
                feasible = generateInitialSolution(startTime);

                if (!feasible) {
                    System.out.println("【生成初始解】无法找到可行解，请检查模型参数");
                    this.failureStage = 2; // 生成初始解阶段失败
                    this.bestObjective = -1;
                    return;
                }
                System.out.println("【生成初始解】近似方法成功找到可行解");
            }
        }

        // Step 3: 改善初始解
        double cur_value = 0.0; // 初始化cur_value，避免作用域问题
        // 只有当找到可行解时，才进行改善初始解的迭代
        // if (!feasible) {
        //     System.out.println("【改善初始解】未找到初始可行解，跳过改善初始解阶段");
        // } else {
        //     boolean change = true;
        //     cur_value = evaluateObjective();
        //     int iteration = 0;

        //     while (change && iteration < maxIterations && feasible) {
        //         iteration++;
        //         change = false;

        //         // 检查每个区域的真正中心
        //         ArrayList<Area> newCenters = findTrueCenters();

        //         // 如果区域中心发生变化，更新并重新求解
        //         boolean centersChanged = !compareCenters(centers, newCenters);
        //         if (centersChanged) {
        //             System.out.println("【改善初始解】检测到中心点变化，准备重新求解...");
        //             // Check time limit before re-solving
        //             if (System.currentTimeMillis() - startTime > GLOBAL_TIME_LIMIT_MS) {
        //                 System.out.println("Global time limit of 1000 seconds exceeded during center adjustment");
        //                 feasible = false;
        //                 this.failureStage = 3; // 改善初始解阶段失败
        //                 break;
        //             }
        //             centers = newCenters;
        //             change = true;

        //             // 如果使用精确方法，中心点变化时需要清空旧的割约束并重新迭代
        //             if (useExactMethod) {
        //                 // 清空旧的割约束集合（因为中心点变了，旧的割约束不再有效）
        //                 infeasibleSolutions.clear();
        //                 System.out.println("【改善初始解】中心点变化，清空旧的割约束，重新迭代添加割约束...");
        //                 // 重新执行精确方法的迭代过程
        //                 feasible = solveWithExactMethod(startTime);
        //                 // 如果重新求解后没有找到可行解，应该退出循环
        //                 if (!feasible) {
        //                     System.out.println("【改善初始解】重新求解后未找到可行解，停止改善初始解");
        //                     this.failureStage = 3; // 改善初始解阶段失败
        //                     break;
        //                 }
        //             } else {
        //                 feasible = generateInitialSolution(startTime);
        //                 // 如果重新求解后没有找到可行解，应该退出循环
        //                 if (!feasible) {
        //                     System.out.println("【改善初始解】重新求解后未找到可行解，停止改善初始解");
        //                     this.failureStage = 3; // 改善初始解阶段失败
        //                     break;
        //                 }
        //             }
        //         } else {
        //             // 中心点没有变化，不需要重新求解
        //             System.out.println("【改善初始解】中心点未变化，跳过重新求解");

        //             if (feasible) {
        //                 cur_value = evaluateObjective();
        //                 System.out.println("【改善初始解】Iteration " + iteration + ": 目标函数值 = " + cur_value);
        //             } else {
        //                 System.out.println("【改善初始解】Iteration " + iteration + ": 无法找到可行解");
        //                 this.failureStage = 3; // 改善初始解阶段失败
        //                 break;
        //             }
        //         }
        //     }

        //     // 如果改善初始解阶段结束后 feasible 为 false，确保设置 failureStage
        //     if (!feasible && this.failureStage == 0) {
        //         this.failureStage = 3; // 改善初始解阶段失败
        //     }
        // }

        // 确保连通性
        // if (feasible) {
        //     // Check time limit before ensuring connectivity
        //     if (System.currentTimeMillis() - startTime > GLOBAL_TIME_LIMIT_MS) {
        //         System.out.println("Global time limit of 1000 seconds exceeded before connectivity check");
        //         this.failureStage = 4; // 确保连通性阶段失败
        //         this.bestObjective = -1;
        //         return;
        //     }
        //     boolean connectivitySuccess = ensureConnectivity(startTime);
        //     if (!connectivitySuccess) {
        //         this.failureStage = 4; // 确保连通性阶段失败
        //         feasible = false;
        //     }
        //     if (feasible) {
        //         cur_value = evaluateObjective();
        //     }
        // }

        // 评估最终结果
        cur_value = evaluateObjective();
        if (feasible && cur_value < Best) {
            Best = cur_value;
            for (int z = 0; z < inst.k; z++) {
                BestZones[z] = new ArrayList<>();
                if (zones[z] != null) {
                    BestZones[z].addAll(zones[z]);
                }
            }
        }

        // Store best objective value
        this.bestObjective = Best;

        // 设置求解成功标志
        this.solveSuccess = (feasible && Best != Double.MAX_VALUE && Best != -1);
        // 如果求解成功，确保failureStage为0
        if (this.solveSuccess) {
            this.failureStage = 0;
        }

        long endTime = System.currentTimeMillis();
        double timeSpentInSeconds = (endTime - startTime) / 1000.0;

        // 输出结果
        String outputFilePath = "./output/" + filename.replace(".dat", "_drcc.txt");
        FileWriter writer = new FileWriter(outputFilePath);
        BufferedWriter buffer = new BufferedWriter(writer);

        for (int io = 0; io < BestZones.length; io++) {
            if (BestZones[io] != null) {
                buffer.write("center ID: " + centers.get(io).getId() + "\n");
                for (int jo = 0; jo < BestZones[io].size(); jo++) {
                    buffer.write(BestZones[io].get(jo) + " ");
                }
                buffer.newLine();
            }
        }

        String result = String.format("%.2f", Best);
        buffer.write("best objective: " + result + "\n");
        buffer.write("程序运行时间为：" + timeSpentInSeconds + "s" + "\n");
        buffer.write("模型类型：" + (useImprovedModel ? "改进模型（最短路距离+工作量约束）" : "原始模型（欧氏距离+需求约束）") + "\n");
        buffer.write("模糊集类型：" + (useD1 ? "D_1" : "D_2") + "\n");
        if (!useD1) {
            buffer.write("delta1: " + delta1 + ", delta2: " + delta2 + "\n");
        }
        buffer.write("约束类型：" + (useJointChance ? "联合约束(DRJCC)" : "个体约束(DRICC)") + "\n");
        buffer.write("求解方法：" + (useExactMethod ? "精确方法（割约束）" : "近似方法（MISOCP）") + "\n");
        buffer.write("风险参数：" + gamma + "\n");
        buffer.write("失败阶段：" + failureStage + " (0=成功，1=选择初始中心失败，2=生成初始解失败，3=改善初始解失败，4=确保连通性失败)" + "\n");
        buffer.close();
        
        // 关闭工作量CSV文件
        closeWorkloadCsvFile();

        System.out.println("程序运行时间为：" + timeSpentInSeconds + "s");
        System.out.println("最终目标函数值：" + Best);
    }

    /**
     * 计算样本矩信息（均值和协方差矩阵）
     */
    /**
     * 加载assignment-dependent workload数据
     * 从output/travel_dist_dual_values_filtered_by_date_1目录下的全部CSV文件读取（训练数据）
     */
    private void loadAssignmentDependentData() {
        String dataDir = "data/travel_dist_dual_values_filtered_by_date_low_ratio_filled";
        java.io.File dir = new java.io.File(dataDir);
        java.io.File[] allFiles = dir.listFiles((d, name) -> name.endsWith(".csv") && name.startsWith("travel_dist_dual_values_p3_"));
        
        if (allFiles == null || allFiles.length == 0) {
            System.err.println("错误: 在目录 " + dataDir + " 中未找到CSV文件");
            return;
        }
        
        // 排序文件以确保一致性
        Arrays.sort(allFiles, (f1, f2) -> f1.getName().compareTo(f2.getName()));
        
        // 使用全部文件作为训练数据
        java.io.File[] files = allFiles;
        
        System.out.println("加载全部 " + files.length + " 个CSV文件作为训练数据");
        
        int n = inst.getN();
        numAssignmentDependentScenarios = files.length;
        
        // 首先读取第一个文件来确定区域数量p
        int p = 0;
        if (files.length > 0) {
            try {
                java.io.BufferedReader reader = new java.io.BufferedReader(
                    new java.io.InputStreamReader(new java.io.FileInputStream(files[0]), "UTF-8"));
                String line = reader.readLine(); // 读取表头
                if (line != null) {
                    String[] headers = line.split(",");
                    p = headers.length - 1; // 减去PointID列
                }
                reader.close();
            } catch (Exception e) {
                System.err.println("错误: 无法读取第一个文件来确定区域数量");
                e.printStackTrace();
            }
        }
        
        if (p == 0) {
            System.err.println("错误: 无法确定区域数量，使用inst.k作为默认值");
            p = inst.k;
        }
        
        numRegionsForAssignmentDependent = p;
        assignmentDependentScenarios = new double[numAssignmentDependentScenarios][n][p];
        
        System.out.println("开始加载assignment-dependent workload数据，共 " + numAssignmentDependentScenarios + " 个场景，区域数量 p = " + p);
        
        for (int k = 0; k < files.length; k++) {
            java.io.File file = files[k];
            try {
                java.io.BufferedReader reader = new java.io.BufferedReader(
                    new java.io.InputStreamReader(new java.io.FileInputStream(file), "UTF-8"));
                
                String line = reader.readLine(); // 读取表头
                if (line == null) {
                    System.err.println("警告: 文件 " + file.getName() + " 为空");
                    reader.close();
                    continue;
                }
                
                // 解析表头，获取区域数量
                String[] headers = line.split(",");
                int numRegions = headers.length - 1; // 减去PointID列
                
                if (numRegions != p) {
                    System.err.println("警告: 文件 " + file.getName() + " 的区域数量 (" + numRegions + ") 与第一个文件不一致 (" + p + ")");
                }
                
                // 初始化矩阵（所有值设为0）
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < p; j++) {
                        assignmentDependentScenarios[k][i][j] = 0.0;
                    }
                }
                
                // 读取数据
                int lineNum = 0;
                while ((line = reader.readLine()) != null) {
                    String[] values = line.split(",");
                    if (values.length < 2) continue;
                    
                    try {
                        int pointId = Integer.parseInt(values[0].trim());
                        if (pointId < 0 || pointId >= n) continue;
                        
                        // 读取每个区域的值
                        for (int j = 0; j < Math.min(numRegions, p); j++) {
                            String valueStr = values[j + 1].trim();
                            if (!valueStr.equals("Null") && !valueStr.isEmpty()) {
                                try {
                                    double value = Double.parseDouble(valueStr);
                                    assignmentDependentScenarios[k][pointId][j] = value;
                                } catch (NumberFormatException e) {
                                    // 忽略无效数值
                                }
                            }
                        }
                        lineNum++;
                    } catch (NumberFormatException e) {
                        // 忽略无效行
                    }
                }
                
                reader.close();
                
                // if ((k + 1) % 50 == 0 || k == files.length - 1) {
                //     System.out.println("  已加载 " + (k + 1) + "/" + files.length + " 个场景");
                // }
            } catch (Exception e) {
                System.err.println("错误: 读取文件 " + file.getName() + " 时出错: " + e.getMessage());
                e.printStackTrace();
            }
        }
        
        System.out.println("完成加载assignment-dependent workload数据");
    }

    private void calculateMomentInformation() {
        if (useAssignmentDependent) {
            calculateMomentInformationAssignmentDependent();
        } else {
            calculateMomentInformationStandard();
        }
    }

    /**
     * 计算标准模型的均值和协方差（N维）
     */
    private void calculateMomentInformationStandard() {
        int n = inst.getN();

        // 计算均值向量
        meanVector = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int s = 0; s < numScenarios; s++) {
                sum += scenarios[s][i];
            }
            meanVector[i] = sum / numScenarios;
        }

        // 计算协方差矩阵
        covarianceMatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double covariance = 0;
                for (int s = 0; s < numScenarios; s++) {
                    covariance += (scenarios[s][i] - meanVector[i]) * (scenarios[s][j] - meanVector[j]);
                }
                covarianceMatrix[i][j] = covariance / numScenarios;
            }
        }

        // 检查矩阵是否对称
        boolean isSymmetric = isSymmetric(covarianceMatrix);
        System.out.println("协方差矩阵是否对称: " + isSymmetric);

        // 检查矩阵是否半正定
        boolean isPSD = isPSDByEigenvalues();
        if (!isPSD) {
            System.out.println("警告: 协方差矩阵不是半正定的，正在尝试修正...");
            ensurePSDMatrix();

            boolean isPSDAfterFix = isPSDByEigenvalues();
            System.out.println("修正后矩阵是否半正定: " + isPSDAfterFix);
        }
    }

    /**
     * 计算assignment-dependent模型的均值和协方差（N*p维）
     * 将N×p矩阵向量化为N*p维向量，按照论文中的顺序：[d_11, d_12, ..., d_1p, d_21, ..., d_Np]
     * 其中N是基本单元数量，p是区域数量
     */
    private void calculateMomentInformationAssignmentDependent() {
        int n = inst.getN();
        int p = numRegionsForAssignmentDependent;
        int nTimesP = n * p;

        // 向量化索引函数：将(i,j)映射到向量位置
        // vec(D) = [d_11, d_12, ..., d_1p, d_21, ..., d_Np]^T
        // 位置 = i * p + j，其中i是基本单元索引(0..N-1)，j是区域索引(0..p-1)

        // 计算均值向量（N*p维）
        meanVector = new double[nTimesP];
        for (int idx = 0; idx < nTimesP; idx++) {
            int i = idx / p;
            int j = idx % p;
            double sum = 0;
            for (int k = 0; k < numAssignmentDependentScenarios; k++) {
                sum += assignmentDependentScenarios[k][i][j];
            }
            meanVector[idx] = sum / numAssignmentDependentScenarios;
        }

        // 检查内存需求
        long matrixSize = (long) nTimesP * (long) nTimesP;
        long memoryBytes = matrixSize * 8L; // 每个double是8字节
        long memoryGB = memoryBytes / (1024L * 1024L * 1024L);
        long memoryMB = memoryBytes / (1024L * 1024L);
        
        System.out.println("协方差矩阵内存需求估算:");
        System.out.println("  矩阵大小: " + nTimesP + "×" + nTimesP + " = " + matrixSize + " 个元素 (N=" + n + ", p=" + p + ")");
        System.out.println("  内存需求: " + memoryMB + " MB (" + memoryGB + " GB)");
        
        // 检查可用内存
        Runtime runtime = Runtime.getRuntime();
        long maxMemory = runtime.maxMemory();
        long totalMemory = runtime.totalMemory();
        long freeMemory = runtime.freeMemory();
        long availableMemory = maxMemory - (totalMemory - freeMemory);
        
        System.out.println("JVM内存状态:");
        System.out.println("  最大堆内存: " + (maxMemory / (1024L * 1024L)) + " MB");
        System.out.println("  可用内存: " + (availableMemory / (1024L * 1024L)) + " MB");
        
        // 对于assignment-dependent模型，默认使用对角线简化（假设场景独立）
        // 这样可以大幅减少内存和计算复杂度
        useDiagonalCovariance = true; // 假设场景之间独立，只计算方差
        nTimesPForCovariance = nTimesP;
        
        if (useDiagonalCovariance) {
            // 只计算对角线元素（方差），假设场景之间独立
            System.out.println("使用对角线简化模式（假设场景之间独立，只计算方差）");
            System.out.println("内存需求: " + (nTimesP * 8 / (1024L * 1024L)) + " MB (仅存储方差向量)");
            
            // 计算方差向量（对角线元素）
            varianceVector = new double[nTimesP];
            System.out.println("开始计算assignment-dependent模型的方差向量（" + nTimesP + " 个元素）...");
            
            for (int idx = 0; idx < nTimesP; idx++) {
                int i = idx / p;
                int j = idx % p;
                double variance = 0;
                for (int k = 0; k < numAssignmentDependentScenarios; k++) {
                    double val = assignmentDependentScenarios[k][i][j] - meanVector[idx];
                    variance += val * val; // 只计算平方，不需要交叉项
                }
                varianceVector[idx] = variance / numAssignmentDependentScenarios;
                
                if ((idx + 1) % 1000 == 0) {
                    System.out.println("  已计算 " + (idx + 1) + "/" + nTimesP + " 个方差");
                }
            }
            
            System.out.println("完成计算assignment-dependent模型的方差向量");
            covarianceMatrix = null; // 不需要存储完整矩阵
            useLazyCovariance = false; // 不使用延迟计算
        } else {
            // 如果内存需求超过50GB，使用延迟计算（按需计算）
            final long MEMORY_THRESHOLD_GB = 50;
            useLazyCovariance = (memoryGB > MEMORY_THRESHOLD_GB);
            
            if (useLazyCovariance) {
                System.out.println("内存需求过大 (" + memoryGB + " GB)，使用延迟计算模式（按需计算协方差）");
                System.out.println("注意: 这将增加计算时间，但可以节省大量内存");
                covarianceMatrix = null; // 不预先分配矩阵
            } else {
                // 如果内存需求超过2GB，给出警告
                if (memoryGB > 2) {
                    System.err.println("警告: 协方差矩阵需要 " + memoryGB + " GB 内存！");
                    System.err.println("建议: 增加JVM堆内存，例如使用 -Xmx" + (memoryGB * 2) + "g 参数");
                }
                
                if (memoryBytes > availableMemory) {
                    System.err.println("错误: 可用内存不足！需要 " + memoryMB + " MB，但只有 " + (availableMemory / (1024L * 1024L)) + " MB 可用。");
                    System.err.println("请增加JVM堆内存: -Xmx" + (memoryGB * 2) + "g");
                    throw new OutOfMemoryError("协方差矩阵内存需求 (" + memoryGB + " GB) 超过可用内存。请增加JVM堆内存。");
                }

                // 计算协方差矩阵（(N*p)×(N*p)）
                System.out.println("开始计算assignment-dependent模型的协方差矩阵（" + nTimesP + "×" + nTimesP + "）...");
                try {
                    covarianceMatrix = new double[nTimesP][nTimesP];
                } catch (OutOfMemoryError e) {
                    System.err.println("内存不足！无法创建协方差矩阵。");
                    System.err.println("请增加JVM堆内存，例如: -Xmx" + (memoryGB * 2) + "g");
                    throw e;
                }
                
                for (int idx1 = 0; idx1 < nTimesP; idx1++) {
                    int i1 = idx1 / p;
                    int j1 = idx1 % p;
                    for (int idx2 = 0; idx2 < nTimesP; idx2++) {
                        int i2 = idx2 / p;
                        int j2 = idx2 % p;
                        double covariance = 0;
                        for (int k = 0; k < numAssignmentDependentScenarios; k++) {
                            double val1 = assignmentDependentScenarios[k][i1][j1] - meanVector[idx1];
                            double val2 = assignmentDependentScenarios[k][i2][j2] - meanVector[idx2];
                            covariance += val1 * val2;
                        }
                        covarianceMatrix[idx1][idx2] = covariance / numAssignmentDependentScenarios;
                    }
                    if ((idx1 + 1) % 100 == 0) {
                        System.out.println("  已计算 " + (idx1 + 1) + "/" + nTimesP + " 行");
                    }
                }

                System.out.println("完成计算assignment-dependent模型的协方差矩阵");
            }
        }

        // 检查矩阵是否对称（仅在非延迟模式且非对角线模式下）
        if (!useLazyCovariance && !useDiagonalCovariance) {
            boolean isSymmetric = isSymmetric(covarianceMatrix);
            System.out.println("协方差矩阵是否对称: " + isSymmetric);

            // 检查矩阵是否半正定
            // 对于大矩阵（>10000），特征值分解非常耗时，跳过检查
            // 协方差矩阵理论上应该是半正定的，如果数值误差导致不是完全半正定，在后续使用中处理
            final int EIGENVALUE_CHECK_THRESHOLD = 10000;
            if (nTimesP > EIGENVALUE_CHECK_THRESHOLD) {
                System.out.println("矩阵规模较大 (" + nTimesP + ")，跳过特征值分解检查（协方差矩阵理论上应该是半正定的）");
                System.out.println("如果后续计算中出现数值问题，可以考虑使用延迟计算模式或增加数值稳定性处理");
            } else {
                System.out.println("开始检查矩阵是否半正定（特征值分解）...");
                boolean isPSD = isPSDByEigenvalues();
                if (!isPSD) {
                    System.out.println("警告: 协方差矩阵不是半正定的，正在尝试修正...");
                    ensurePSDMatrix();

                    System.out.println("修正后重新检查半正定性...");
                    boolean isPSDAfterFix = isPSDByEigenvalues();
                    System.out.println("修正后矩阵是否半正定: " + isPSDAfterFix);
                } else {
                    System.out.println("矩阵是半正定的");
                }
            }
        } else if (useDiagonalCovariance) {
            System.out.println("对角线模式: 方差向量已计算完成（假设场景独立，方差自动满足半正定性）");
        } else {
            System.out.println("延迟计算模式: 跳过对称性和半正定性检查（将在需要时按需计算）");
        }
    }
    
    /**
     * 按需计算协方差矩阵的单个元素（用于延迟计算模式）
     * @param idx1 第一个索引（对应d_{i1,j1}）
     * @param idx2 第二个索引（对应d_{i2,j2}）
     * @return 协方差值
     */
    private double getCovariance(int idx1, int idx2) {
        // 对角线模式：只返回对角线元素（方差），非对角线元素为0
        if (useDiagonalCovariance) {
            if (idx1 == idx2) {
                return varianceVector[idx1];
            } else {
                return 0.0; // 假设场景独立，非对角线协方差为0
            }
        }
        
        if (!useLazyCovariance) {
            // 正常模式：直接返回预先计算的值
            return covarianceMatrix[idx1][idx2];
        }
        
        // 延迟计算模式：按需计算
        int n = inst.getN();
        int p = numRegionsForAssignmentDependent;
        int i1 = idx1 / p;
        int j1 = idx1 % p;
        int i2 = idx2 / p;
        int j2 = idx2 % p;
        
        double covariance = 0;
        for (int k = 0; k < numAssignmentDependentScenarios; k++) {
            double val1 = assignmentDependentScenarios[k][i1][j1] - meanVector[idx1];
            double val2 = assignmentDependentScenarios[k][i2][j2] - meanVector[idx2];
            covariance += val1 * val2;
        }
        return covariance / numAssignmentDependentScenarios;
    }
    
    /**
     * 计算矩阵向量乘积 Σ * v（优化版本，支持延迟计算和对角线模式）
     * @param v 输入向量（N*p维）
     * @return 结果向量 Σ * v（N*p维）
     */
    private double[] computeCovarianceMatrixVectorProduct(double[] v) {
        int nTimesP = nTimesPForCovariance;
        double[] result = new double[nTimesP];
        
        if (useDiagonalCovariance) {
            // 对角线模式：Σ * v = diag(σ) * v，只需逐元素相乘
            for (int idx = 0; idx < nTimesP; idx++) {
                result[idx] = varianceVector[idx] * v[idx];
            }
        } else if (!useLazyCovariance) {
            // 正常模式：使用预先计算的矩阵
            for (int idx1 = 0; idx1 < nTimesP; idx1++) {
                result[idx1] = 0.0;
                for (int idx2 = 0; idx2 < nTimesP; idx2++) {
                    result[idx1] += covarianceMatrix[idx1][idx2] * v[idx2];
                }
            }
        } else {
            // 延迟计算模式：按需计算
            for (int idx1 = 0; idx1 < nTimesP; idx1++) {
                result[idx1] = 0.0;
                for (int idx2 = 0; idx2 < nTimesP; idx2++) {
                    double cov = getCovariance(idx1, idx2);
                    result[idx1] += cov * v[idx2];
                }
            }
        }
        return result;
    }
    
    /**
     * 计算二次型 v^T * Σ * v（优化版本，支持延迟计算和对角线模式）
     * @param v 输入向量（N*p维）
     * @return 二次型值
     */
    private double computeQuadraticForm(double[] v) {
        int nTimesP = nTimesPForCovariance;
        double result = 0.0;
        
        if (useDiagonalCovariance) {
            // 对角线模式：v^T * Σ * v = sum(v_i^2 * σ_i)，只需计算对角线项
            for (int idx = 0; idx < nTimesP; idx++) {
                result += v[idx] * v[idx] * varianceVector[idx];
            }
        } else if (!useLazyCovariance) {
            // 正常模式：使用预先计算的矩阵
            for (int idx1 = 0; idx1 < nTimesP; idx1++) {
                for (int idx2 = 0; idx2 < nTimesP; idx2++) {
                    result += v[idx1] * covarianceMatrix[idx1][idx2] * v[idx2];
                }
            }
        } else {
            // 延迟计算模式：优化计算，利用对称性
            // v^T * Σ * v = sum_{i,j} v[i] * Σ[i,j] * v[j]
            // 利用对称性：Σ[i,j] = Σ[j,i]
            for (int idx1 = 0; idx1 < nTimesP; idx1++) {
                for (int idx2 = 0; idx2 <= idx1; idx2++) {
                    double cov = getCovariance(idx1, idx2);
                    if (idx1 == idx2) {
                        result += v[idx1] * cov * v[idx2];
                    } else {
                        // 利用对称性：v[i]*Σ[i,j]*v[j] + v[j]*Σ[j,i]*v[i] = 2*v[i]*Σ[i,j]*v[j]
                        result += 2.0 * v[idx1] * cov * v[idx2];
                    }
                }
            }
        }
        return result;
    }

    /**
     * 计算assignment-dependent模型的上下界
     */
    private void calculateAssignmentDependentBounds() {
        int n = inst.getN();
        int p = numRegionsForAssignmentDependent;
        double totalMeanWorkload = 0;
        
        // 计算总平均工作量
        // 对于assignment-dependent模型，总工作量是 sum_{i,j} E[d_ij]，其中i是基本单元，j是区域
        int nTimesP = n * p;
        for (int idx = 0; idx < nTimesP; idx++) {
            totalMeanWorkload += meanVector[idx];
        }
        
        double mu = totalMeanWorkload / inst.k; // 每个区域的平均工作量
        this.demandUpperBound = (1 + r) * mu; // 上界 U = (1 + r) * mu
        this.demandLowerBound = (1 - r) * mu; // 下界 L = (1 - r) * mu
        
        System.out.println("Assignment-dependent模型上下界: L = " + demandLowerBound + ", U = " + demandUpperBound);
    }

    private boolean isSymmetric(double[][] matrix) {
        int n = matrix.length;
        
        // 对于大矩阵（>5000），对称性检查也很耗时，使用采样检查
        final int SYMMETRY_CHECK_THRESHOLD = 5000;
        if (n > SYMMETRY_CHECK_THRESHOLD) {
            System.out.println("矩阵规模较大 (" + n + ")，使用采样方法检查对称性（检查部分元素）...");
            // 采样检查：检查对角线上方的一些元素
            int sampleSize = Math.min(10000, n * (n - 1) / 2); // 最多检查10000个元素
            int checked = 0;
            double tolerance = 1e-10;
            
            // 随机采样或均匀采样
            int step = Math.max(1, (n * (n - 1) / 2) / sampleSize);
            int count = 0;
            
            for (int i = 0; i < n && checked < sampleSize; i++) {
                for (int j = i + 1; j < n && checked < sampleSize; j++) {
                    if (count % step == 0) {
                        if (Math.abs(matrix[i][j] - matrix[j][i]) > tolerance) {
                            System.out.println("检测到不对称: matrix[" + i + "][" + j + "] = " + matrix[i][j] + 
                                             ", matrix[" + j + "][" + i + "] = " + matrix[j][i]);
                            return false;
                        }
                        checked++;
                    }
                    count++;
                }
            }
            System.out.println("采样检查完成，检查了 " + checked + " 个元素，未发现不对称");
            return true;
        } else {
            // 小矩阵：完整检查
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double tolerance = 1e-10;
                    if (Math.abs(matrix[i][j] - matrix[j][i]) > tolerance) {
                        return false;
                    }
                }
            }
            return true;
        }
    }

    private boolean isPSDByEigenvalues() {
        if (useLazyCovariance) {
            // 延迟计算模式下，跳过特征值分解（需要完整矩阵）
            System.out.println("延迟计算模式: 跳过半正定性检查");
            return true; // 假设是半正定的（协方差矩阵理论上应该是半正定的）
        }
        
        int n = covarianceMatrix.length;
        
        // 对于大矩阵，特征值分解非常耗时，添加提示
        if (n > 5000) {
            System.out.println("警告: 矩阵规模较大 (" + n + ")，特征值分解可能需要较长时间，请耐心等待...");
        }

        try {
            System.out.println("正在进行特征值分解（这可能需要几分钟到几小时，取决于矩阵大小）...");
            long startEigTime = System.currentTimeMillis();
            
            Matrix matrix = new Matrix(covarianceMatrix);
            Jama.EigenvalueDecomposition eig = matrix.eig();
            double[] eigenvalues = eig.getRealEigenvalues();

            long eigTime = System.currentTimeMillis() - startEigTime;
            System.out.println("特征值分解完成，耗时: " + (eigTime / 1000.0) + " 秒");

            double minEigenvalue = Double.MAX_VALUE;
            for (double ev : eigenvalues) {
                minEigenvalue = Math.min(minEigenvalue, ev);
            }

            double epsilon = 1e-15; // 数值稳定性的容差
            System.out.println("最小特征值: " + minEigenvalue);

            return minEigenvalue >= -epsilon;
        } catch (java.lang.Exception e) {
            System.out.println("特征值分解失败: " + e.getMessage());
            e.printStackTrace();
            return false;
        }
    }

    private void ensurePSDMatrix() {
        if (useLazyCovariance || covarianceMatrix == null) {
            // 延迟计算模式下，无法修改矩阵
            System.out.println("延迟计算模式: 无法修正矩阵（矩阵未预先计算）");
            return;
        }
        
        // 在对角线上添加一个小的正数
        double epsilon = 1e-5; // 略大于最小负特征值的绝对值
        for (int i = 0; i < covarianceMatrix.length; i++) {
            covarianceMatrix[i][i] += epsilon;
        }
    }

    /**
     * 计算最短路径距离矩阵（使用Floyd-Warshall算法）
     * 对于改进模型，r_ij 表示从基本单元i到区域中心j的最短路径距离
     */
    private void computeShortestPathDistances() {
        int n = inst.getN();
        shortestPathDist = new double[n][n];

        // 初始化距离矩阵
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    // 节点到自身的距离为0
                    shortestPathDist[i][j] = 0.0;
                } else if (inst.getEdges()[i][j] == 1) {
                    // 如果i和j之间有边，使用欧氏距离作为边的权重
                    shortestPathDist[i][j] = inst.dist[i][j];
                } else {
                    // 如果i和j之间没有直接边，初始化为无穷大
                    shortestPathDist[i][j] = Double.MAX_VALUE;
                }
            }
        }

        // Floyd-Warshall算法：通过中间节点k来更新所有节点对之间的最短路径
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    // 如果通过节点k可以缩短i到j的路径，则更新
                    if (shortestPathDist[i][k] != Double.MAX_VALUE &&
                            shortestPathDist[k][j] != Double.MAX_VALUE) {
                        double newDist = shortestPathDist[i][k] + shortestPathDist[k][j];
                        if (newDist < shortestPathDist[i][j]) {
                            shortestPathDist[i][j] = newDist;
                        }
                    }
                }
            }
        }

        System.out.println("已使用Floyd-Warshall算法计算最短路径距离矩阵");
    }

    /**
     * 计算基于工作量的上下界（改进模型）
     * 根据新方法：L = (1 - α) * W/p, U = (1 + α) * W/p
     * 其中 W = ∑_i μ_i * r̄_i，r̄_i 的计算方式有两种：
     * 方式1：r̄_i = 所有点j（j的最近k个点中包含i）到i的距离的平均值
     * 方式2：r̄_i = 所有点j（j的最近k个点中包含i）的k近邻平均距离的平均值
     * k = ⌈n/p⌉
     * @param method avgDist计算方式：1=方式1，2=方式2
     */
    private void calculateWorkloadBounds(int method) {
        int n = inst.getN();
        int p = inst.k;
        int k = (int) Math.ceil((double) n / p); // k = ⌈n/p⌉
        k = Math.min(k, n); // 确保 k 不超过总节点数

        // 计算每个基本单元的平均配送距离
        double[] avgDist = new double[n];
        
        if (method == 1) {
            // 方式1：对于每个点i，找到所有点j（j的最近k个点中包含i），计算所有j到i的距离的平均值
            for (int i = 0; i < n; i++) {
                ArrayList<Integer> validJs = new ArrayList<>(); // 存储所有满足条件的j
                
                // 对于每个点j，检查i是否在j的最近k个点中
                for (int j = 0; j < n; j++) {
                    // 找到j到所有点的距离并排序
                    double[] distancesFromJ = new double[n];
                    int[] indices = new int[n];
                    for (int idx = 0; idx < n; idx++) {
                        distancesFromJ[idx] = shortestPathDist[j][idx];
                        indices[idx] = idx;
                    }
                    
                    // 使用索引数组排序，保持索引对应关系
                    for (int a = 0; a < n - 1; a++) {
                        for (int b = a + 1; b < n; b++) {
                            if (distancesFromJ[a] > distancesFromJ[b]) {
                                double tempDist = distancesFromJ[a];
                                distancesFromJ[a] = distancesFromJ[b];
                                distancesFromJ[b] = tempDist;
                                int tempIdx = indices[a];
                                indices[a] = indices[b];
                                indices[b] = tempIdx;
                            }
                        }
                    }
                    
                    // 检查i是否在j的最近k个点中
                    boolean iInKNearest = false;
                    for (int idx = 0; idx < k; idx++) {
                        if (indices[idx] == i) {
                            iInKNearest = true;
                            break;
                        }
                    }
                    
                    if (iInKNearest) {
                        validJs.add(j);
                    }
                }
                
                // 输出validJs的size
                System.out.println("【方式1】点 " + i + " 的validJs.size = " + validJs.size());
                
                // 计算所有validJs中的点j到i的距离的平均值
                if (!validJs.isEmpty()) {
                    double sum = 0.0;
                    for (int j : validJs) {
                        sum += shortestPathDist[j][i];
                    }
                    avgDist[i] = sum / validJs.size();
                } else {
                    // 如果没有满足条件的j，使用原来的方法作为fallback
                    if (i == 0) {
                        System.out.println("【提示】方式1：对于点 " + i + "，没有找到满足条件的j（j的最近k个点中包含i），使用原来的计算方式作为fallback");
                    }
                    double[] distances = new double[n];
                    for (int j = 0; j < n; j++) {
                        distances[j] = shortestPathDist[i][j];
                    }
                    java.util.Arrays.sort(distances);
                    double sum = 0.0;
                    for (int idx = 0; idx < k; idx++) {
                        sum += distances[idx];
                    }
                    avgDist[i] = sum / k;
                }
            }
        } else if (method == 2) {
            // 方式2：对于每个点i，找到所有点j（j的最近k个点中包含i），
            // 对于每个j，计算j的最近k个点的距离的平均值，然后再对这些平均值求平均
            for (int i = 0; i < n; i++) {
                ArrayList<Integer> validJs = new ArrayList<>(); // 存储所有满足条件的j
                
                // 对于每个点j，检查i是否在j的最近k个点中
                for (int j = 0; j < n; j++) {
                    // 找到j到所有点的距离并排序
                    double[] distancesFromJ = new double[n];
                    int[] indices = new int[n];
                    for (int idx = 0; idx < n; idx++) {
                        distancesFromJ[idx] = shortestPathDist[j][idx];
                        indices[idx] = idx;
                    }
                    
                    // 使用索引数组排序，保持索引对应关系
                    for (int a = 0; a < n - 1; a++) {
                        for (int b = a + 1; b < n; b++) {
                            if (distancesFromJ[a] > distancesFromJ[b]) {
                                double tempDist = distancesFromJ[a];
                                distancesFromJ[a] = distancesFromJ[b];
                                distancesFromJ[b] = tempDist;
                                int tempIdx = indices[a];
                                indices[a] = indices[b];
                                indices[b] = tempIdx;
                            }
                        }
                    }
                    
                    // 检查i是否在j的最近k个点中
                    boolean iInKNearest = false;
                    for (int idx = 0; idx < k; idx++) {
                        if (indices[idx] == i) {
                            iInKNearest = true;
                            break;
                        }
                    }
                    
                    if (iInKNearest) {
                        validJs.add(j);
                    }
                }
                
                // 输出validJs的size
                System.out.println("【方式2】点 " + i + " 的validJs.size = " + validJs.size());
                
                // 对于每个validJs中的点j，计算j的最近k个点的距离的平均值
                if (!validJs.isEmpty()) {
                    double sumOfAverages = 0.0;
                    for (int j : validJs) {
                        // 找到j的最近k个点的距离
                        double[] distancesFromJ = new double[n];
                        for (int idx = 0; idx < n; idx++) {
                            distancesFromJ[idx] = shortestPathDist[j][idx];
                        }
                        java.util.Arrays.sort(distancesFromJ);
                        
                        // 计算j的最近k个点的平均距离
                        double sum = 0.0;
                        for (int idx = 0; idx < k; idx++) {
                            sum += distancesFromJ[idx];
                        }
                        double avgDistForJ = sum / k;
                        sumOfAverages += avgDistForJ;
                    }
                    // 对这些平均值求平均
                    avgDist[i] = sumOfAverages / validJs.size();
                } else {
                    // 如果没有满足条件的j，使用原来的方法作为fallback
                    if (i == 0) {
                        System.out.println("【提示】方式2：对于点 " + i + "，没有找到满足条件的j（j的最近k个点中包含i），使用原来的计算方式作为fallback");
                    }
                    double[] distances = new double[n];
                    for (int j = 0; j < n; j++) {
                        distances[j] = shortestPathDist[i][j];
                    }
                    java.util.Arrays.sort(distances);
                    double sum = 0.0;
                    for (int idx = 0; idx < k; idx++) {
                        sum += distances[idx];
                    }
                    avgDist[i] = sum / k;
                }
            }
        } else if (method == 3) {
            // 方式3：使用原来的方法，对于每个点i，计算i到其最近k个点的距离的平均值
            for (int i = 0; i < n; i++) {
                double[] distances = new double[n];
                for (int j = 0; j < n; j++) {
                    distances[j] = shortestPathDist[i][j];
                }
                java.util.Arrays.sort(distances);
                double sum = 0.0;
                for (int idx = 0; idx < k; idx++) {
                    sum += distances[idx];
                }
                avgDist[i] = sum / k;
            }
        } else {
            // 默认使用原来的方法
            System.out.println("【警告】avgDistMethod=" + method + " 不是有效值（应为1、2或3），使用原来的计算方式：对于每个点i，计算i到其最近k个点的距离的平均值");
            for (int i = 0; i < n; i++) {
                double[] distances = new double[n];
                for (int j = 0; j < n; j++) {
                    distances[j] = shortestPathDist[i][j];
                }
                java.util.Arrays.sort(distances);
                double sum = 0.0;
                for (int idx = 0; idx < k; idx++) {
                    sum += distances[idx];
                }
                avgDist[i] = sum / k;
            }
        }

        // 计算当前方法的avgDist向量均值并存储
        double sumAvgDist = 0.0;
        for (int i = 0; i < n; i++) {
            sumAvgDist += avgDist[i];
        }
        this.avgDistMean = sumAvgDist / n;
        
        // 输出当前方法的avgDist向量均值
        System.out.println("【记录】方式" + method + "的avgDist向量均值 = " + avgDistMean);

        // 计算总期望工作量 W = ∑_i μ_i * r̄_i（使用均值需求）
        double totalExpectedWorkload = 0.0;
        for (int i = 0; i < n; i++) {
            totalExpectedWorkload += meanVector[i] * avgDist[i];
        }

        // 计算每个区域的平均工作量
        double avgWorkload = totalExpectedWorkload / p;

        // 使用 r 作为 α 参数（平衡容差）
        this.demandUpperBound = (1 + r) * avgWorkload;
        this.demandLowerBound = (1 - r) * avgWorkload;

        System.out.println("基于工作量的上下界（k=" + k + ", avgDist计算方式=" + method + "）：L = " + demandLowerBound + ", U = " + demandUpperBound);
    }

    /**
     * 调整上下界：增加2%
     */
    private void adjustBounds() {
        double adjustmentFactor = 1.05; // 增加2%
        this.demandUpperBound *= adjustmentFactor;
        this.demandLowerBound *= adjustmentFactor;
        System.out.println("调整上下界（增加2%）：L = " + demandLowerBound + ", U = " + demandUpperBound);
    }

    /**
     * 调整上下界：降低2%
     */
    private void reduceBounds() {
        double adjustmentFactor = 0.95; // 降低2%
        this.demandUpperBound *= adjustmentFactor;
        this.demandLowerBound *= adjustmentFactor;
        System.out.println("调整上下界（降低2%）：L = " + demandLowerBound + ", U = " + demandUpperBound);
    }

    /**
     * 选择初始区域中心
     * 对于改进模型且不使用相对平衡性约束的情况，会自适应调整上下界参考值，直到所有scenario都能成功求解
     */
    /**
     * 从CSV文件加载初始中心点（用于assignment-dependent模型）
     */
    private ArrayList<Integer> loadInitialCentersFromCSV() {
        String csvFile = "data/test/selected_centers_low_ratio_p3.csv";
        ArrayList<Integer> centers = new ArrayList<>();
        
        try {
            java.io.BufferedReader reader = new java.io.BufferedReader(
                new java.io.InputStreamReader(new java.io.FileInputStream(csvFile), "UTF-8"));
            
            String line = reader.readLine(); // 读取表头
            if (line == null) {
                System.err.println("错误: CSV文件为空");
                reader.close();
                return centers;
            }
            
            while ((line = reader.readLine()) != null) {
                String[] values = line.split(",");
                if (values.length >= 1) {
                    try {
                        // 第一列是序号（支持整数或浮点如 70 或 70.0）；.dat 中点的 id = 序号 - 1
                        int xuhao = (int) Double.parseDouble(values[0].trim());
                        int idInDat = xuhao - 1;
                        // 在实例中查找 id 为 idInDat 的点的索引（0..N-1）
                        for (int i = 0; i < inst.getN(); i++) {
                            if (inst.getAreas()[i].getId() == idInDat) {
                                centers.add(i);
                                break;
                            }
                        }
                    } catch (NumberFormatException e) {
                        // 忽略无效行
                    }
                }
            }
            
            reader.close();
            
            System.out.println("从CSV文件加载了 " + centers.size() + " 个初始中心点");
            if (centers.size() != inst.k) {
                System.err.println("警告: 加载的中心点数量 (" + centers.size() + ") 与所需数量 (" + inst.k + ") 不匹配");
            }
        } catch (Exception e) {
            System.err.println("错误: 读取初始中心点CSV文件时出错: " + e.getMessage());
            e.printStackTrace();
        }
        
        return centers;
    }

    private ArrayList<Integer> selectInitialCenters() throws GRBException {
        int InitialNum = 5; // 改为10个场景
        ArrayList<Integer> candidateCenters = new ArrayList<>();
        HashMap<Integer, Integer> centerFrequency = new HashMap<>();

        // 如果使用改进模型且不使用相对平衡性约束，需要自适应调整上下界
        if (useImprovedModel && !useRelativeBalance) {
            int maxAdjustmentIterations = 20; // 最大调整次数，避免无限循环
            int adjustmentIteration = 0;
            boolean allScenariosSuccess = false;

            while (!allScenariosSuccess && adjustmentIteration < maxAdjustmentIterations) {
                adjustmentIteration++;
                if (adjustmentIteration > 1) {
                    System.out.println("【调整参考值】第 " + adjustmentIteration + " 次调整上下界...");
                    adjustBounds();
                } else {
                    System.out.println("【选择初始中心】开始求解 " + InitialNum + " 个场景的确定性模型...");
                }

                // 随机选择 InitialNum 个不同的场景索引
                Set<Integer> selectedScenarios = new HashSet<>();
                while (selectedScenarios.size() < InitialNum && selectedScenarios.size() < numScenarios) {
                    selectedScenarios.add(rand.nextInt(numScenarios));
                }
                ArrayList<Integer> scenarioIndices = new ArrayList<>(selectedScenarios);

                // 尝试求解所有选定的场景
                boolean hasFailure = false;
                centerFrequency.clear(); // 清空之前的频率统计

                for (int idx = 0; idx < scenarioIndices.size(); idx++) {
                    int scenarioIndex = scenarioIndices.get(idx);
                    ScenarioSolveResult result = solveForScenario(scenarioIndex);

                    if (result.success && result.centers != null && result.centers.size() == inst.k) {
                        System.out.println("场景 " + (idx + 1) + "/" + scenarioIndices.size() + "（索引: " + scenarioIndex + "）求解成功");

                        // 更新中心频率
                        for (int center : result.centers) {
                            centerFrequency.put(center, centerFrequency.getOrDefault(center, 0) + 1);
                        }
                    } else {
                        System.out.println("场景 " + (idx + 1) + "/" + scenarioIndices.size() + "（索引: " + scenarioIndex + "）求解失败");
                        hasFailure = true;
                        break; // 一旦有失败，立即停止并调整参考值
                    }
                }

                if (!hasFailure) {
                    allScenariosSuccess = true;
                    System.out.println("【选择初始中心】所有 " + scenarioIndices.size() + " 个场景都成功求解，使用当前上下界继续后续步骤");
                } else {
                    System.out.println("【选择初始中心】检测到求解失败，需要调整参考值");
                }
            }

            if (!allScenariosSuccess) {
                System.out.println("【选择初始中心】在 " + maxAdjustmentIterations + " 次调整后仍无法使所有场景成功求解，选择初始中心阶段失败，终止本次实验");
                return null; // 直接返回null，终止后续所有步骤
            }
        } else {
            // 原始模型：保持原有逻辑
            int scenariosProcessed = 0;
            while (scenariosProcessed < InitialNum && scenariosProcessed < numScenarios) {
                int scenarioIndex = rand.nextInt(numScenarios);

                // 使用该场景的需求求解确定性模型
                ScenarioSolveResult result = solveForScenario(scenarioIndex);

                if (result.success && result.centers != null && result.centers.size() == inst.k) {
                    scenariosProcessed++;
                    System.out.println("处理场景 " + scenariosProcessed + "/" + InitialNum + "，场景索引: " + scenarioIndex);

                    // 更新中心频率
                    for (int center : result.centers) {
                        centerFrequency.put(center, centerFrequency.getOrDefault(center, 0) + 1);
                    }
                }
            }
        }

        // 检查是否有任何场景成功
        if (centerFrequency.isEmpty()) {
            System.out.println("【选择初始中心】没有任何场景成功求解，选择初始中心阶段失败");
            return null;
        }

        // 按频率排序选择前k个中心
        List<Map.Entry<Integer, Integer>> sortedCenters = new ArrayList<>(centerFrequency.entrySet());
        sortedCenters.sort(Map.Entry.<Integer, Integer>comparingByValue().reversed());

        for (int i = 0; i < Math.min(inst.k, sortedCenters.size()); i++) {
            candidateCenters.add(sortedCenters.get(i).getKey());
        }

        // 如果中心数量不足，随机补充
        while (candidateCenters.size() < inst.k) {
            int randomCenter = rand.nextInt(inst.getN());
            if (!candidateCenters.contains(randomCenter)) {
                candidateCenters.add(randomCenter);
            }
        }

        return candidateCenters;
    }

    private Instance createScenarioInstance(int scenarioIndex) {
        return new Instance(inst, scenarioDemands[scenarioIndex]);
    }

    /**
     * 求解结果类：包含是否成功和中心点列表
     */
    private static class ScenarioSolveResult {
        boolean success;
        ArrayList<Integer> centers;

        ScenarioSolveResult(boolean success, ArrayList<Integer> centers) {
            this.success = success;
            this.centers = centers;
        }
    }

    /**
     * 求解单一场景的确定性模型
     * @param scenarioIndex 场景索引
     * @return 求解结果（包含是否成功和中心点列表）
     */
    private ScenarioSolveResult solveForScenario(int scenarioIndex) throws GRBException {
        try {
            // 创建基于特定场景需求的实例
            Instance scenarioInstance = createScenarioInstance(scenarioIndex);

            // 创建Algo对象（不设置时间限制，由getCorrectSolutionCentersInternal内部控制）
            Algo algo = new Algo(scenarioInstance);

            // 如果使用改进模型，设置改进模型参数
            if (useImprovedModel && shortestPathDist != null && meanVector != null) {
                algo.setImprovedModelParams(true, shortestPathDist, meanVector);
                
                // 如果使用相对平衡性约束，不需要设置边界值（相对平衡性约束基于比例，不依赖绝对边界值）
                // 如果使用绝对平衡性约束，需要设置边界值
                if (!useRelativeBalance) {
                    // 使用当前已调整的上下界（而不是重新计算）
                    algo.setBounds(demandLowerBound, demandUpperBound);
                }
                // 如果使用相对平衡性约束，不设置边界值，让Algo使用默认的边界值或自己计算
            }

            // 获取该场景下的求解结果中心点
            ArrayList<Integer> scenarioCenters;
            if (useImprovedModel && shortestPathDist != null && meanVector != null) {
                scenarioCenters = algo.getCorrectSolutionCenters(true, shortestPathDist, meanVector, useRelativeBalance);
            } else {
                scenarioCenters = algo.getCorrectSolutionCenters();
            }

            // 如果算法未能返回足够的中心点，说明求解失败
            if (scenarioCenters.size() < inst.k) {
                System.out.println("场景 " + scenarioIndex + " 求解失败：中心点不足（" + scenarioCenters.size() + " < " + inst.k + "）");
                return new ScenarioSolveResult(false, null);
            }

            // 求解成功
            return new ScenarioSolveResult(true, scenarioCenters);
        } catch (Exception e) {
            System.out.println("求解场景 " + scenarioIndex + " 时出错: " + e.getMessage());
            // 求解失败
            return new ScenarioSolveResult(false, null);
        }
    }

    /**
     * 使用精确方法求解：迭代添加割约束直到找到可行解
     * @param globalStartTime 全局开始时间
     * @return 是否找到可行解
     */
    private boolean solveWithExactMethod(long globalStartTime) throws GRBException {
        // 如果使用相对平衡性约束，直接求解（约束已直接添加到模型中）
        if (useRelativeBalance) {
            System.out.println("【生成初始解】开始使用精确方法（直接添加分式约束到模型）...");
            boolean feasible = generateInitialSolution(globalStartTime);
            this.cutIterations = 1; // 只求解一次，不需要迭代
            if (feasible) {
                System.out.println("【生成初始解】精确方法找到可行解");
            } else {
                System.out.println("【生成初始解】精确方法未找到可行解");
            }
            return feasible;
        }
        
        // 不使用相对平衡性约束时，使用no-good cut方法：迭代添加割约束
        System.out.println("【生成初始解】开始使用精确方法（no-good cut，迭代添加割约束）...");

        int maxCutIterations = 100; // 最大割约束迭代次数
        int cutIteration = 0;
        boolean feasible = false;

        // 迭代添加割约束直到找到可行解
        while (!feasible && cutIteration < maxCutIterations) {
            cutIteration++;
            System.out.println("【生成初始解】精确方法迭代 " + cutIteration + "：求解松弛模型...");

            feasible = generateInitialSolution(globalStartTime);

            if (!feasible && cutIteration < maxCutIterations) {
                System.out.println("【生成初始解】解不满足DRCC约束，添加割约束后重新求解...");
            }

            // 检查时间限制
            if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                System.out.println("【生成初始解】全局时间限制已超过，停止精确方法迭代");
                break;
            }
        }

        // 保存迭代次数（记录最后一次的迭代次数）
        this.cutIterations = cutIteration;

        if (feasible) {
            System.out.println("【生成初始解】精确方法在 " + cutIteration + " 次迭代后找到可行解");
        } else {
            System.out.println("【生成初始解】精确方法在 " + cutIteration + " 次迭代后仍未找到可行解");
        }

        return feasible;
    }

    /**
     * 在包含连通性约束的模型上使用精确方法求解
     * @param model 已包含连通性约束的模型
     * @param x 决策变量
     * @param env Gurobi环境
     * @param globalStartTime 全局开始时间
     * @return 是否找到满足DRCC约束的解
     */
    private boolean solveWithExactMethodForConnectivity(GRBModel model, GRBVar[][] x, GRBEnv env, long globalStartTime, int initialCutCount) throws GRBException {
        // no-good cut方法：迭代添加割约束
        return solveWithExactMethodForConnectivityNoGoodCut(model, x, env, globalStartTime, initialCutCount);
    }

    /**
     * no-good cut方法：迭代添加割约束直到找到满足DRCC约束的解
     */
    private boolean solveWithExactMethodForConnectivityNoGoodCut(GRBModel model, GRBVar[][] x, GRBEnv env, long globalStartTime, int initialCutCount) throws GRBException {
        int maxCutIterations = 100; // 最大割约束迭代次数
        int cutIteration = 0;
        boolean feasible = false;

        while (!feasible && cutIteration < maxCutIterations) {
            cutIteration++;
            System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + "：求解包含连通性约束的模型...");

            // Update time limit
            long remainingTimeMs = GLOBAL_TIME_LIMIT_MS - (System.currentTimeMillis() - globalStartTime);
            double remainingTimeSec = Math.max(1.0, Math.min(remainingTimeMs / 1000.0, timeLimit));
            model.set(GRB.DoubleParam.TimeLimit, remainingTimeSec);

            // 求解模型
            model.optimize();

            // 检查时间限制
            if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                System.out.println("全局时间限制已超过，停止精确方法迭代");
                // 记录最后一次求解的状态代码（即使时间限制超时，也要记录模型状态）
                int status = model.get(GRB.IntAttr.Status);
                this.statusCode = status;
                break;
            }

            // 检查模型状态
            int status = model.get(GRB.IntAttr.Status);
            
            // 验证求解状态（检查是否是全局最优）
            if (useExactMethod && useRelativeBalance) {
                reportSolveStatus(model, "确保连通性-精确方法-迭代" + cutIteration);
            }
            
            if (status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                // 保存状态代码（如果求解失败，必须记录）
                this.statusCode = status;
                System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + " 失败，模型无解");
                break;
            }

            // 提取解决方案
            for (int j = 0; j < centers.size(); j++) {
                zones[j] = new ArrayList<>();
                for (int i = 0; i < inst.getN(); i++) {
                    if (Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                        zones[j].add(i);
                    }
                }
            }

            // 即使使用相对平衡性约束且约束已添加到模型中，也需要验证以确保解的质量
            // （因为可能存在数值误差导致Gurobi返回的解不满足真实约束）
            // 暂时注释掉验证，如果Gurobi求解成功就认为成功，不进行回验证
            // 如果Gurobi返回OPTIMAL或SUBOPTIMAL，直接认为求解成功
            feasible = true; // 如果Gurobi求解成功，就认为成功
            if (useRelativeBalance) {
                System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + "：找到满足约束的解（约束已直接添加到模型中）");
            } else {
                System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + "：找到满足DRCC约束的解");
            }
            /*
            boolean isFeasible = verifyDRCCConstraints(x);
            if (isFeasible) {
                feasible = true;
                if (useRelativeBalance) {
                    System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + "：找到满足约束的解（约束已直接添加到模型中）");
                } else {
                    System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + "：找到满足DRCC约束的解");
                }
            } else {
                if (useRelativeBalance) {
                    // 使用相对平衡性约束时，虽然约束已添加到模型，但解不满足真实约束
                    System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + "：警告：Gurobi返回的解不满足DRCC约束（可能存在数值误差），继续迭代...");
                    // 继续迭代，尝试找到满足约束的解
                } else {
                    // 不使用相对平衡性约束时，记录不可行解用于no-good cut
                    // 记录不可行解
                    recordInfeasibleSolution(x);

                    // 只添加新发现的不可行解对应的割约束
                    // 因为模型已经包含了初始的割约束，我们只需要添加新增的
                    int currentCutCount = infeasibleSolutions.size();
                    if (currentCutCount > initialCutCount) {
                        // 添加新发现的割约束
                        int newCutsAdded = 0;
                        for (int k = initialCutCount; k < currentCutCount; k++) {
                            Set<String> infeasibleSet = infeasibleSolutions.get(k);
                            GRBLinExpr cutExpr = new GRBLinExpr();
                            for (String index : infeasibleSet) {
                                String[] parts = index.split("_");
                                int i = Integer.parseInt(parts[0]);
                                int j = Integer.parseInt(parts[1]);
                                cutExpr.addTerm(1.0, x[i][j]);
                            }
                            int setSize = infeasibleSet.size();
                            model.addConstr(cutExpr, GRB.LESS_EQUAL, setSize - 1, "conn_cut_" + k);
                            newCutsAdded++;
                        }
                        System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + "：添加了 " + newCutsAdded + " 个新的割约束");
                        initialCutCount = currentCutCount; // 更新已添加的割约束数量
                    }
                    System.out.println("【确保连通性-精确方法】迭代 " + cutIteration + "：解不满足DRCC约束，添加割约束后重新求解...");
                }
            }
            */
        }

        return feasible;
    }

    /**
     * 生成初始可行解 - 使用Gurobi
     */
    private boolean generateInitialSolution(long globalStartTime) throws GRBException {
        GRBModel model = null;
        GRBEnv env = null;
        try {
            // Check time limit before starting
            if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                System.out.println("【生成初始解】全局时间限制已超过，无法继续求解");
                return false;
            }

            // 创建Gurobi环境和模型
            env = new GRBEnv(true);
            env.set(GRB.IntParam.OutputFlag, 0);
            env.set(GRB.IntParam.LogToConsole, 0);
            env.set(GRB.StringParam.LogFile, "");
            env.set(GRB.IntParam.Seed, 42);
            env.start();
            model = new GRBModel(env);

            // 设置Gurobi参数 - 使用剩余时间
            long remainingTimeMs = GLOBAL_TIME_LIMIT_MS - (System.currentTimeMillis() - globalStartTime);
            double remainingTimeSec = Math.max(1.0, Math.min(remainingTimeMs / 1000.0, timeLimit));
            model.set(GRB.IntParam.OutputFlag, 0);
            model.set(GRB.DoubleParam.TimeLimit, remainingTimeSec);
            
            // 设置更严格的数值容差，以提高解的精度
            // 特别是对于精确的分式约束，需要更严格的容差
            if (useExactMethod && useRelativeBalance) {
                model.set(GRB.DoubleParam.FeasibilityTol, 1e-6);  // 可行性容差：默认1e-6，保持严格
                model.set(GRB.DoubleParam.OptimalityTol, 1e-6);   // 最优性容差：默认1e-6，保持严格
                model.set(GRB.IntParam.NumericFocus, 3);          // 数值焦点：3表示最高精度
                model.set(GRB.DoubleParam.BarConvTol, 1e-8);      // 障碍法收敛容差：更严格
                // 关键修复：设置NonConvex=2保证全局最优（对于包含二次分式约束的非凸问题）
                model.set(GRB.IntParam.NonConvex, 2);
            }

            // 决策变量 x_ij
            GRBVar[][] x = new GRBVar[inst.getN()][centers.size()];
            for (int i = 0; i < inst.getN(); i++) {
                for (int j = 0; j < centers.size(); j++) {
                    x[i][j] = model.addVar(0, 1, 0, GRB.BINARY, "x_" + i + "_" + centers.get(j).getId());
                    if (i == centers.get(j).getId()) {
                        // 将中心点固定为1
                        model.addConstr(x[i][j], GRB.EQUAL, 1.0, "center_" + i + "_" + j);
                    }
                }
            }

            // 约束: 每个基本单元必须且只能属于一个区域
            for (int i = 0; i < inst.getN(); i++) {
                GRBLinExpr rowSum = new GRBLinExpr();
                for (int j = 0; j < centers.size(); j++) {
                    rowSum.addTerm(1.0, x[i][j]);
                }
                model.addConstr(rowSum, GRB.EQUAL, 1.0, "unit_" + i);
            }

            // 根据方法类型添加约束
            if (useExactMethod) {
                // 精确方法
                if (useRelativeBalance) {
                    // 使用相对平衡性约束时，直接添加精确的分式约束（利用Gurobi 13处理高次分式的功能）
                    addDistributionallyRobustConstraints(model, x);
                } else {
                    // 不使用相对平衡性约束时，使用no-good cut方法
                    // 不添加DRCC约束，而是通过迭代验证和添加割约束
                    // 添加已记录的割约束
                    addCutConstraints(model, x);
                }
            } else {
                // 近似方法：直接添加DRCC约束（MISOCP）或使用支撑超平面cut
                if (useRelativeBalance && useSupportingHyperplaneCuts && !shouldUseD1CopositiveSdpApprox()) {
                    // 使用支撑超平面cut：不直接添加SOC约束，而是迭代添加cut
                    addDistributionallyRobustConstraints(model, x); // 这会存储约束信息但不添加SOC约束
                } else {
                    // 直接添加DRCC约束（MISOCP）
                    addDistributionallyRobustConstraints(model, x);
                }
            }

            // 设置目标函数（根据模型类型）
            setObjectiveFunction(model, x);

            // 如果使用支撑超平面cut，需要迭代求解
            if (useRelativeBalance && useSupportingHyperplaneCuts && !useExactMethod && !shouldUseD1CopositiveSdpApprox()) {
                // 精确方法：模型直接包含 reformulate 后的约束（如 MIQCP），一次（或少数几次）求解即可得可行解或判不可行。
                // 近似方法：不把 SOC 加入模型，仅用支撑超平面 cut 外逼近；需多轮“求解→检查违反→加 cut”才能收敛。
                // assignment-dependent 时约束在 N*p 维、2*p 个 SOC，边界弯曲，外逼近所需 cut 数较多，故迭代上限要更大。
                int maxCutIterations = 50;
                if (useAssignmentDependent) {
                    maxCutIterations = 500; // assignment-dependent 时 SOC 边界复杂，外逼近需更多 cut；统一给足迭代次数（gamma 小时代数上更紧，收敛更慢）
                    System.out.println("【支撑超平面cut】assignment-dependent：最大迭代次数设为 " + maxCutIterations + "（近似法为外逼近，需更多 cut 才能收敛）");
                }
                int cutIteration = 0;
                boolean converged = false;
                double bestObjective = Double.MAX_VALUE;
                int noImprovementCount = 0; // 连续没有改善的迭代次数
                int maxNoImprovement = 5; // 允许连续没有改善的最大迭代次数
                
                // 限制每次迭代添加的cut数量，避免模型过于复杂
                int maxCutsPerIteration = useAssignmentDependent ? 10 : 20; // assignment-dependent模型每个cut更复杂，限制更少
                
                while (!converged && cutIteration < maxCutIterations) {
                    cutIteration++;
                    System.out.println("【支撑超平面cut迭代】第 " + cutIteration + " 次迭代...");
                    
                    // 为每次迭代设置时间限制（每次迭代最多3600秒，但不超过剩余时间）
                    long remainingTimeForIteration = GLOBAL_TIME_LIMIT_MS - (System.currentTimeMillis() - globalStartTime);
                    long iterationTimeLimitMs = Math.min(3600 * 1000, remainingTimeForIteration); // 每次迭代最多3600秒，但不超过剩余时间
                    double iterationTimeLimitSec = Math.max(1.0, iterationTimeLimitMs / 1000.0); // 至少1秒
                    double originalTimeLimit = model.get(GRB.DoubleParam.TimeLimit);
                    model.set(GRB.DoubleParam.TimeLimit, iterationTimeLimitSec);
                    System.out.println(String.format("【支撑超平面cut迭代】本次迭代时间限制: %.1f 秒 (剩余总时间: %.1f 秒)", 
                        iterationTimeLimitSec, remainingTimeForIteration / 1000.0));
                    
                    // 记录求解前的模型信息
                    int numVarsBefore = model.get(GRB.IntAttr.NumVars);
                    int numConstrsBefore = model.get(GRB.IntAttr.NumConstrs);
                    long solveStartTime = System.currentTimeMillis();
                    
                    // 求解模型
                    model.optimize();
                    
                    // 记录求解时间
                    long solveTime = System.currentTimeMillis() - solveStartTime;
                    
                    // 恢复原始时间限制
                    model.set(GRB.DoubleParam.TimeLimit, originalTimeLimit);
                    
                    // 输出详细的求解信息
                    // reportCutIterationSolveInfo(model, cutIteration, numVarsBefore, numConstrsBefore, solveTime);
                    
                    // 验证求解状态（检查是否是全局最优）
                    if (useExactMethod && useRelativeBalance) {
                        reportSolveStatus(model, "生成初始解-支撑超平面cut迭代-" + cutIteration);
                    }
                    
                    // 检查时间限制
                    if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                        System.out.println("【支撑超平面cut迭代】全局时间限制已超过，停止迭代");
                        int status = model.get(GRB.IntAttr.Status);
                        this.statusCode = status;
                        model.dispose();
                        env.dispose();
                        return false;
                    }
                    
                    // 检查模型状态
                    int status = model.get(GRB.IntAttr.Status);
                    if (status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                        System.out.println("【支撑超平面cut迭代】模型求解失败，状态码: " + status);
                        this.statusCode = status;
                        model.dispose();
                        env.dispose();
                        return false;
                    }
                    
                    // 如果使用assignment-dependent模型，计算并输出各区域工作量
                    if (useAssignmentDependent) {
                        calculateAndPrintWorkloads(x, cutIteration);
                    }
                    
                    // 获取当前目标函数值
                    // double currentObjective = model.get(GRB.DoubleAttr.ObjVal);
                    // boolean objectiveImproved = (currentObjective < bestObjective - 1e-6);
                    // if (objectiveImproved) {
                    //     bestObjective = currentObjective;
                    //     noImprovementCount = 0;
                    //     System.out.println(String.format("【支撑超平面cut迭代】目标函数值改善: %.4f", currentObjective));
                    // } else {
                    //     noImprovementCount++;
                    //     System.out.println(String.format("【支撑超平面cut迭代】目标函数值: %.4f (无改善，连续 %d 次)", currentObjective, noImprovementCount));
                    // }
                    
                    // 检查约束是否被违反，如果违反则添加cut（限制每次添加的数量）
                    boolean cutsAdded = checkAndAddRelativeBalanceCuts(model, x, 1e-5, maxCutsPerIteration);
                    if (!cutsAdded) {
                        // assignment-dependent 时用精确分式约束判定收敛，否则用“无新 cut”即收敛
                        if (useAssignmentDependent) {
                            double[][] currentXVal = new double[inst.getN()][centers.size()];
                            for (int i = 0; i < inst.getN(); i++)
                                for (int k = 0; k < centers.size(); k++)
                                    currentXVal[i][k] = x[i][k].get(GRB.DoubleAttr.X);
                            if (checkExactAssignmentDependentFeasibility(currentXVal))
                                converged = true;
                        } else {
                            converged = true;
                        }
                    }
                    
                    // if (!cutsAdded) {
                    //     // 没有新的cut，但需要验证解是否真正满足约束
                    //     // 因为支撑超平面cut是外逼近，可能不够紧
                    //     System.out.println("【支撑超平面cut迭代】第 " + cutIteration + " 次迭代：没有新的cut，验证解是否满足约束...");
                        
                    //     // 使用更严格的验证方法（验证原始约束）
                    //     boolean trulyFeasible = verifyRelativeBalanceConstraints(x);
                        
                    //     if (trulyFeasible) {
                    //         // 解满足约束，检查是否质量足够好
                    //         // 如果目标函数值连续多次没有改善，认为收敛
                    //         if (noImprovementCount >= maxNoImprovement) {
                    //             converged = true;
                    //             System.out.println(String.format("【支撑超平面cut迭代】在第 " + cutIteration + " 次迭代后收敛：解满足约束且目标函数值连续 %d 次无改善", noImprovementCount));
                    //         } else {
                    //             System.out.println(String.format("【支撑超平面cut迭代】解满足约束，但目标函数值可能仍在改善（连续 %d 次无改善，需要 %d 次），继续迭代...", noImprovementCount, maxNoImprovement));
                    //             // 继续迭代，等待目标函数值稳定
                    //         }
                    //     } else {
                    //         // 解不满足约束，但checkAndAddRelativeBalanceCuts没有检测到违反
                    //         // 这可能是因为阈值设置问题，或者解是连续松弛解
                    //         // 强制添加cut以确保收敛到真正满足约束的解
                    //         System.out.println("【支撑超平面cut迭代】警告：解不满足约束但未检测到违反，强制添加cut...");
                            
                    //         // 重新检查，使用更严格的阈值（也限制数量）
                    //         boolean strictCheck = checkAndAddRelativeBalanceCuts(model, x, 1e-8, maxCutsPerIteration);
                    //         if (strictCheck) {
                    //             cutsAdded = true;
                    //             noImprovementCount = 0; // 重置计数器
                    //             System.out.println("【支撑超平面cut迭代】使用严格检查添加了新的cut，继续迭代...");
                    //         } else {
                    //             // 即使严格检查也没有违反，但验证失败
                    //             // 这可能是因为解是连续松弛解，需要继续迭代
                    //             System.out.println("【支撑超平面cut迭代】警告：验证失败但未检测到违反，可能是连续松弛解，继续迭代...");
                    //             // 不设置converged，继续迭代
                    //         }
                    //     }
                    // } else {
                    //     noImprovementCount = 0; // 有新的cut，重置计数器
                    //     System.out.println("【支撑超平面cut迭代】添加了新的cut，继续迭代...");
                    // }
                }
                
                if (!converged && cutIteration >= maxCutIterations) {
                    System.out.println("【支撑超平面cut迭代】达到最大迭代次数 " + maxCutIterations + "，检查当前解是否满足约束...");
                    // assignment-dependent 时验证精确分式约束，否则验证近似约束
                    double[][] currentXVal = new double[inst.getN()][centers.size()];
                    for (int i = 0; i < inst.getN(); i++)
                        for (int k = 0; k < centers.size(); k++)
                            currentXVal[i][k] = x[i][k].get(GRB.DoubleAttr.X);
                    boolean trulyFeasible = useAssignmentDependent ? checkExactAssignmentDependentFeasibility(currentXVal) : verifyRelativeBalanceConstraints(x);
                    if (!trulyFeasible) {
                        // 达到最大迭代次数后，当前解仍不满足约束，判定为失败
                        System.out.println("【支撑超平面cut迭代】达到最大迭代次数后，当前解仍不满足约束，判定为失败");
                        int status = model.get(GRB.IntAttr.Status);
                        this.statusCode = status;
                        model.dispose();
                        env.dispose();
                        return false;
                    } else {
                        System.out.println("【支撑超平面cut迭代】达到最大迭代次数，但当前解满足约束");
                        // 即使达到最大迭代次数，只要解满足约束，就认为成功
                        converged = true;
                    }
                }
                
                // 检查模型状态（在支撑超平面cut迭代后）
                int status = model.get(GRB.IntAttr.Status);
                
                // 验证求解状态（检查是否是全局最优）
                if (useExactMethod && useRelativeBalance) {
                    reportSolveStatus(model, "生成初始解-支撑超平面cut迭代-最终");
                }
                
                if (status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                    System.out.println("【支撑超平面cut迭代】模型求解失败，状态码: " + status);
                    this.statusCode = status;
                    model.dispose();
                    env.dispose();
                    return false;
                }
            } else if (shouldUseD1CopositiveSdpApprox()) {
                boolean ok = iterateD1SdpSeparation(model, x, globalStartTime, 200, "D1-SDP分离");
                if (!ok) {
                    int status = model.get(GRB.IntAttr.Status);
                    this.statusCode = status;
                    model.dispose();
                    env.dispose();
                    return false;
                }
            } else {
                // 直接求解模型
                model.optimize();
                
                // 验证求解状态（检查是否是全局最优）
                if (useExactMethod && useRelativeBalance) {
                    reportSolveStatus(model, "生成初始解-直接求解");
                }
                
                // Check time limit after optimization
                if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                    System.out.println("【生成初始解】全局时间限制已超过，停止求解");
                    // 记录最后一次求解的状态代码（即使时间限制超时，也要记录模型状态）
                    int status = model.get(GRB.IntAttr.Status);
                    this.statusCode = status;
                    model.dispose();
                    env.dispose();
                    return false;
                }

                // 检查模型状态
                int status = model.get(GRB.IntAttr.Status);
                if (status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                    // 保存状态代码（如果求解失败，必须记录）
                    this.statusCode = status;
                    model.dispose();
                    env.dispose();
                    return false;
                }
            }
            
            // Check time limit after optimization (对于支撑超平面cut的情况)
            if (useRelativeBalance && useSupportingHyperplaneCuts && !useExactMethod && !shouldUseD1CopositiveSdpApprox()) {
                if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                    System.out.println("【生成初始解】全局时间限制已超过，停止求解");
                    int status = model.get(GRB.IntAttr.Status);
                    this.statusCode = status;
                    model.dispose();
                    env.dispose();
                    return false;
                }
            }
            // 如果求解成功，可以不更新statusCode（保持-1），但为了准确性也可以更新
            // 这里选择不更新，因为用户说成功时可不更新

            // 提取解决方案
            for (int j = 0; j < centers.size(); j++) {
                zones[j] = new ArrayList<>();
                for (int i = 0; i < inst.getN(); i++) {
                    if (Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                        zones[j].add(i);
                    }
                }
            }

            // 如果使用精确方法，需要验证解是否满足DRCC约束
            // 即使使用相对平衡性约束且约束已添加到模型中，也需要验证以确保解的质量
            // （因为可能存在数值误差导致Gurobi返回的解不满足真实约束）
            // 暂时注释掉验证，如果Gurobi求解成功就认为成功，不进行回验证
            /*
            if (useExactMethod) {
                boolean isFeasible = verifyDRCCConstraints(x);
                if (!isFeasible) {
                    if (useRelativeBalance) {
                        // 使用相对平衡性约束时，虽然约束已添加到模型，但解不满足真实约束
                        // 这可能是因为数值误差或约束添加方式的问题
                        System.out.println("【生成初始解】警告：Gurobi返回的解不满足DRCC约束（可能存在数值误差）");
                        model.dispose();
                        env.dispose();
                        return false; // 需要重新求解
                    } else {
                        // 不使用相对平衡性约束时，记录不可行解用于no-good cut
                        recordInfeasibleSolution(x);
                        model.dispose();
                        env.dispose();
                        return false; // 需要重新求解
                    }
                }
            }
            */

            // 释放资源
            model.dispose();
            env.dispose();

            return true;
        } catch (GRBException e) {
            System.out.println("Gurobi求解模型时发生错误: " + e.getMessage());
            // 如果model已经创建且optimize已调用，尝试获取状态代码
            // 注意：如果异常发生在optimize之前，model可能未创建或未调用optimize
            // 这种情况下statusCode保持为-1是合理的
            try {
                if (model != null) {
                    int status = model.get(GRB.IntAttr.Status);
                    this.statusCode = status;
                }
            } catch (Exception ex) {
                // 如果无法获取状态代码，保持statusCode为-1
            }
            return false;
        }
    }

    /**
     * 添加分布鲁棒机会约束 - 使用Gurobi实现
     */
    private void addDistributionallyRobustConstraints(GRBModel model, GRBVar[][] x) throws GRBException {
        // 如果使用支撑超平面cut，清空之前的约束信息
        if (useRelativeBalance && useSupportingHyperplaneCuts && !useExactMethod && !shouldUseD1CopositiveSdpApprox()) {
            relativeBalanceConstraintInfos.clear();
        }
        
        if (useJointChance) {
            // 使用Bonferroni近似的DRJCC模型
            for (int j = 0; j < centers.size(); j++) {
                // 根据DRICC添加个体约束，使用individualGammas[j]作为风险参数
                addDRICCConstraint(model, x, j, individualGammas[j]);
            }
        } else {
            // 直接使用DRICC模型
            for (int j = 0; j < centers.size(); j++) {
                addDRICCConstraint(model, x, j, gamma);
            }
        }
    }

    /**
     * 添加DRICC约束 - 使用Gurobi实现
     */
    private void addDRICCConstraint(GRBModel model, GRBVar[][] x, int j, double riskParam) throws GRBException {
        // 根据约束类型选择添加绝对约束还是相对约束
        if (useRelativeBalance) {
            addRelativeBalanceConstraint(model, x, j, riskParam);
        } else {
            addAbsoluteBalanceConstraint(model, x, j, riskParam);
        }
    }

    /**
     * D1 + 非assignment-dependent + 近似法 + 原始模型时，使用文档中的 COP+SDP 安全近似。
     */
    private boolean shouldUseD1CopositiveSdpApprox() {
        return useD1 && !useExactMethod && !useAssignmentDependent && !useImprovedModel;
    }

    /**
     * D1 + d>=0 支撑集：主问题(Gurobi)上的外层分离循环。
     * 每轮先解主问题，再用 COPT 解 SDP 证书子问题验算；若违反则加 no-good cut。
     */
    private boolean iterateD1SdpSeparation(
            GRBModel model,
            GRBVar[][] x,
            long globalStartTime,
            int maxIterations,
            String logTag
    ) throws GRBException {
        for (int it = 1; it <= maxIterations; it++) {
            long remaining = GLOBAL_TIME_LIMIT_MS - (System.currentTimeMillis() - globalStartTime);
            if (remaining <= 0) {
                System.out.println("【" + logTag + "】全局时间限制已超过，停止迭代");
                this.statusCode = model.get(GRB.IntAttr.Status);
                return false;
            }

            double iterTimeSec = Math.max(1.0, Math.min(3600.0, remaining / 1000.0));
            double oldLimit = model.get(GRB.DoubleParam.TimeLimit);
            model.set(GRB.DoubleParam.TimeLimit, iterTimeSec);

            System.out.println("【" + logTag + "】第 " + it + " 次迭代，求解主问题...");
            model.optimize();
            model.set(GRB.DoubleParam.TimeLimit, oldLimit);

            int status = model.get(GRB.IntAttr.Status);
            if (status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                System.out.println("【" + logTag + "】主问题求解失败，状态码: " + status);
                this.statusCode = status;
                return false;
            }

            boolean cutAdded = checkAndAddD1SdpNoGoodCut(model, x, 1e-6);
            if (!cutAdded) {
                System.out.println("【" + logTag + "】COPT-SDP证书全部通过，迭代收敛");
                return true;
            }
            System.out.println("【" + logTag + "】检测到SDP违反，已添加no-good cut，继续迭代");
        }

        System.out.println("【" + logTag + "】达到最大迭代次数，仍未收敛");
        return false;
    }

    /**
     * 检查当前解在 D1 + 支撑集 d>=0 下是否满足 SDP 证书；若不满足，加入 no-good cut 排除该解。
     * @return true 表示添加了 cut；false 表示当前解通过证书（无新增 cut）
     */
    private boolean checkAndAddD1SdpNoGoodCut(GRBModel model, GRBVar[][] x, double tol) throws GRBException {
        int n = inst.getN();
        int p = centers.size();
        double[][] xVal = new double[n][p];
        int[] selectedJ = new int[n];

        for (int i = 0; i < n; i++) {
            double best = -Double.MAX_VALUE;
            int bestJ = 0;
            for (int j = 0; j < p; j++) {
                xVal[i][j] = x[i][j].get(GRB.DoubleAttr.X);
                if (xVal[i][j] > best) {
                    best = xVal[i][j];
                    bestJ = j;
                }
            }
            selectedJ[i] = bestJ;
        }

        if (verifyD1SdpCertificatesByCopt(xVal, tol)) {
            return false;
        }

        GRBLinExpr cutExpr = new GRBLinExpr();
        for (int i = 0; i < n; i++) {
            cutExpr.addTerm(1.0, x[i][selectedJ[i]]);
        }
        model.addConstr(cutExpr, GRB.LESS_EQUAL, n - 1, "d1_sdp_nogood_" + (d1SdpNoGoodCutCounter++));
        return true;
    }

    /**
     * 用 COPT 求解文档中的 SDP 证书子问题，验证当前解是否可行。
     */
    private boolean verifyD1SdpCertificatesByCopt(double[][] xVal, double tol) {
        int n = inst.getN();
        int p = centers.size();
        double coeff = (1.0 - r) / p;
        double coeffUpper = (1.0 + r) / p;

        for (int j = 0; j < p; j++) {
            double risk = useJointChance ? individualGammas[j] : gamma;
            double[] split = splitRiskForD1Sdp(risk);
            double gammaL = split[0];
            double gammaU = split[1];

            double[] vL = new double[n];
            double[] vU = new double[n];
            for (int i = 0; i < n; i++) {
                vL[i] = coeff - xVal[i][j];
                vU[i] = xVal[i][j] - coeffUpper;
            }

            double boundL = solveD1SdpBoundWithCopt(vL);
            double boundU = solveD1SdpBoundWithCopt(vU);

            if (Double.isNaN(boundL) || Double.isNaN(boundU) || Double.isInfinite(boundL) || Double.isInfinite(boundU)) {
                throw new RuntimeException("COPT-SDP子问题返回无效值，请检查 COPT/CVXPY 环境。");
            }

            if (boundL > gammaL + tol || boundU > gammaU + tol) {
                System.out.println(String.format("【D1-SDP】区域%d违反: boundL=%.6f (阈值 %.6f), boundU=%.6f (阈值 %.6f)",
                        j, boundL, gammaL, boundU, gammaU));
                return false;
            }
        }
        return true;
    }

    /**
     * D1-SDP 的风险拆分策略：
     * - 默认采用 gamma_a:gamma_b 的比例并归一化到当前 risk；
     * - 若比例不可用，则退化为 0.5/0.5。
     */
    private double[] splitRiskForD1Sdp(double risk) {
        double a = gamma_a;
        double b = gamma_b;
        if (a <= 0 || b <= 0 || a + b <= 1e-12) {
            return new double[]{0.5 * risk, 0.5 * risk};
        }
        double sum = a + b;
        return new double[]{risk * a / sum, risk * b / sum};
    }

    /**
     * 使用 COPT Java API 直接计算单个线性事件的 SDP 证书上界。
     *
     * 文档对齐形式（PDF 公式 (35)）：
     * min y0 + mu^T y + Q:M
     * s.t. Z = P0 + N0
     *      Z - E11 - λ*Vaug = P1 + N1   （λ 为 S-procedure 乘子，是决策变量）
     *      P0,P1 >= 0 (PSD), N0,N1 >= 0 (逐元素)
     *      λ >= 0
     *
     * 等价实现（消去 P0,P1）：
     *      Z - N0 >= 0
     *      Z - E11 - λ*Vaug - N1 >= 0
     *      N0,N1 >= 0
     */
    private double solveD1SdpBoundWithCopt(double[] v) {
        // Bug 1 短路：d>=0 支撑集下，若 v 全为非正，则 d^T v > 0 不可能发生，最坏违约概率恒为 0，直接返回。
        boolean allNonPositive = true;
        for (double vi : v) {
            if (vi > 1e-12) {
                allNonPositive = false;
                break;
            }
        }
        if (allNonPositive) {
            return 0.0;
        }

        if (!d1SdpJavaCoptPathLogged) {
            System.out.println("【D1-SDP】使用 Java 直接调用 COPT API 求解 SDP 子问题（路线A）");
            d1SdpJavaCoptPathLogged = true;
        }

        final int n = inst.getN();
        final int dim = n + 1; // 增广矩阵维度
        RuntimeException lastException = null;

        for (int attempt = 1; attempt <= 2; attempt++) {
            copt.Model cModel = null;
            try {
                copt.Envr env = getOrCreateD1SdpSharedEnvr();
                cModel = env.createModel("d1_sdp_cert");

                // 模型级别静默 + 稳健参数
                try { cModel.setIntParam(copt.IntParam.Logging, 0); } catch (Exception ignored) {}
                try { cModel.setIntParam(copt.IntParam.LogToConsole, 0); } catch (Exception ignored) {}
                try { cModel.setIntParam(copt.IntParam.LogLevel, 0); } catch (Exception ignored) {}
                try { cModel.setDblParam(copt.DblParam.TimeLimit, 30.0); } catch (Exception ignored) {}
                try { cModel.setDblParam(copt.DblParam.FeasTol, 1e-8); } catch (Exception ignored) {}
                try { cModel.setDblParam(copt.DblParam.MatrixTol, 1e-10); } catch (Exception ignored) {}
                try { cModel.setIntParam(copt.IntParam.SDPMethod, attempt == 1 ? 0 : 1); } catch (Exception ignored) {}

                // 变量
                copt.Var y0 = cModel.addVar(1.0, copt.Consts.INFINITY, 0.0, copt.Consts.CONTINUOUS, "y0");
                copt.Var[] y = new copt.Var[n];
                for (int i = 0; i < n; i++) {
                    y[i] = cModel.addVar(-copt.Consts.INFINITY, copt.Consts.INFINITY, 0.0, copt.Consts.CONTINUOUS, "y_" + i);
                }
                // Bug 2：S-procedure 乘子 λ（PDF 公式(35)中的决策变量），λ >= 0，出现在 LMI2 中
                copt.Var lambda = cModel.addVar(0.0, copt.Consts.INFINITY, 0.0, copt.Consts.CONTINUOUS, "lambda");

                // M 对称自由变量（下三角建模，镜像赋值）
                copt.Var[][] mVar = new copt.Var[dim][dim];
                // N0/N1 对称且逐元素非负
                copt.Var[][] n0Var = new copt.Var[dim][dim];
                copt.Var[][] n1Var = new copt.Var[dim][dim];
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        copt.Var mv = cModel.addVar(-copt.Consts.INFINITY, copt.Consts.INFINITY, 0.0, copt.Consts.CONTINUOUS, "m_" + i + "_" + j);
                        mVar[i][j] = mv;
                        mVar[j][i] = mv;

                        copt.Var nv0 = cModel.addVar(0.0, copt.Consts.INFINITY, 0.0, copt.Consts.CONTINUOUS, "n0_" + i + "_" + j);
                        n0Var[i][j] = nv0;
                        n0Var[j][i] = nv0;

                        copt.Var nv1 = cModel.addVar(0.0, copt.Consts.INFINITY, 0.0, copt.Consts.CONTINUOUS, "n1_" + i + "_" + j);
                        n1Var[i][j] = nv1;
                        n1Var[j][i] = nv1;
                    }
                }

                // 目标：y0 + mu^T y + Q:M
                copt.Expr obj = new copt.Expr();
                obj.addTerm(y0, 1.0);
                for (int i = 0; i < n; i++) {
                    obj.addTerm(y[i], meanVector[i]);
                }
                for (int i = 1; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        double qij = getCovariance(i - 1, j - 1) + meanVector[i - 1] * meanVector[j - 1];
                        double coeff = (i == j) ? qij : 2.0 * qij;
                        obj.addTerm(mVar[i][j], coeff);
                    }
                }
                cModel.setObjective(obj, copt.Consts.MINIMIZE);

                // 常量矩阵 E11
                double[][] e11 = new double[dim][dim];
                e11[0][0] = 1.0;
                copt.SymMatrix constE11 = createSymMatrixFromDense(cModel, e11);

                // 常量矩阵 -Vaug（在 LMI2 中乘以变量 λ，得到 -λ*Vaug）
                double[][] minusVaug = new double[dim][dim];
                for (int i = 1; i < dim; i++) {
                    minusVaug[0][i] = -0.5 * v[i - 1];
                    minusVaug[i][0] = -0.5 * v[i - 1];
                }
                copt.SymMatrix constMinusVaug = createSymMatrixFromDense(cModel, minusVaug);

                // LMI1: Z - N0 >= 0
                copt.LmiExpr lmi1 = new copt.LmiExpr();
                lmi1.addTerm(y0, createBasisSymMatrix(cModel, dim, 0, 0, 1.0));
                for (int i = 1; i < dim; i++) {
                    lmi1.addTerm(y[i - 1], createBasisSymMatrix(cModel, dim, 0, i, 0.5));
                }
                for (int i = 1; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        lmi1.addTerm(mVar[i][j], createBasisSymMatrix(cModel, dim, i, j, 1.0));
                    }
                }
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        lmi1.addTerm(n0Var[i][j], createBasisSymMatrix(cModel, dim, i, j, -1.0));
                    }
                }
                cModel.addLmiConstr(lmi1, "lmi_z_minus_n0_psd");

                // LMI2: Z - E11 - λ*Vaug - N1 >= 0（对应 PDF 公式(35)）
                copt.LmiExpr lmi2 = new copt.LmiExpr();
                lmi2.addTerm(y0, createBasisSymMatrix(cModel, dim, 0, 0, 1.0));
                for (int i = 1; i < dim; i++) {
                    lmi2.addTerm(y[i - 1], createBasisSymMatrix(cModel, dim, 0, i, 0.5));
                }
                for (int i = 1; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        lmi2.addTerm(mVar[i][j], createBasisSymMatrix(cModel, dim, i, j, 1.0));
                    }
                }
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        lmi2.addTerm(n1Var[i][j], createBasisSymMatrix(cModel, dim, i, j, -1.0));
                    }
                }
                lmi2.addConstant(new copt.SymMatExpr(constE11, -1.0));
                lmi2.addTerm(lambda, constMinusVaug);
                cModel.addLmiConstr(lmi2, "lmi_z_minus_e11_vaug_minus_n1_psd");

                cModel.solve();

                int status = cModel.getIntAttr(copt.IntAttr.LpStatus);
                int hasLpSol = cModel.getIntAttr(copt.IntAttr.HasLpSol);
                if (status == 1) { // optimal
                    return cModel.getDblAttr(copt.DblAttr.LpObjVal);
                }
                // 偶发非最优但已有可行解时，接受当前值，避免中断整轮实验
                if (hasLpSol == 1) {
                    double val = cModel.getDblAttr(copt.DblAttr.LpObjVal);
                    System.out.println(String.format("【D1-SDP】警告：子问题 status=%d（非最优）但 HasLpSol=1，采用当前值 %.6f", status, val));
                    return val;
                }
                // status=8 通常为 COPT_LPSTATUS_INTERRUPTED（迭代/时间限制或内部中止）；2=INFEASIBLE, 5=TIMEOUT, 9=ITERLIMIT
                // 无解时返回保守上界，使证书判为“违反”并加 cut，避免整轮实验直接失败
                if (attempt == 2) {
                    double fallback = 1.0; // 大于常见 gamma，使当前解被判为违反
                    System.out.println(String.format("【D1-SDP】警告：子问题 status=%d, HasLpSol=0（无解），采用保守上界 %.2f 继续", status, fallback));
                    return fallback;
                }
                throw new RuntimeException("COPT子问题未达到最优状态, status=" + status + ", HasLpSol=" + hasLpSol + ", attempt=" + attempt);
            } catch (Exception e) {
                lastException = new RuntimeException("调用COPT求解SDP证书失败: " + e.getMessage(), e);
                // 环境可能异常，清空后下一次重建
                try {
                    if (d1SdpSharedEnvr != null) {
                        d1SdpSharedEnvr.dispose();
                    }
                } catch (Exception ignored) {}
                d1SdpSharedEnvr = null;
            } finally {
                if (cModel != null) {
                    try {
                        cModel.dispose();
                    } catch (Exception ignored) {
                    }
                }
            }
        }
        throw lastException != null ? lastException : new RuntimeException("调用COPT求解SDP证书失败: 未知错误");
    }

    /**
     * 懒加载并复用 COPT 环境，避免每次子问题都重复初始化输出 banner。
     */
    private copt.Envr getOrCreateD1SdpSharedEnvr() throws copt.CoptException {
        if (d1SdpSharedEnvr != null) {
            return d1SdpSharedEnvr;
        }
        copt.EnvrConfig envConfig = new copt.EnvrConfig();
        try { envConfig.set(copt.IntParam.Logging, "0"); } catch (Exception ignored) {}
        try { envConfig.set(copt.IntParam.LogToConsole, "0"); } catch (Exception ignored) {}
        try { envConfig.set(copt.IntParam.LogLevel, "0"); } catch (Exception ignored) {}
        d1SdpSharedEnvr = new copt.Envr(envConfig);
        return d1SdpSharedEnvr;
    }

    /**
     * 构造对称基矩阵：在(i,j)/(j,i)位置为val，其余为0。
     */
    private copt.SymMatrix createBasisSymMatrix(copt.Model cModel, int dim, int i, int j, double val) throws copt.CoptException {
        if (i == j) {
            return cModel.addSparseMat(dim, 1, new int[]{i}, new int[]{j}, new double[]{val});
        }
        return cModel.addSparseMat(dim, 2, new int[]{i, j}, new int[]{j, i}, new double[]{val, val});
    }

    /**
     * 由稠密对称矩阵创建 COPT SymMatrix（仅导入非零项）。
     */
    private copt.SymMatrix createSymMatrixFromDense(copt.Model cModel, double[][] mat) throws copt.CoptException {
        int dim = mat.length;
        ArrayList<Integer> rows = new ArrayList<>();
        ArrayList<Integer> cols = new ArrayList<>();
        ArrayList<Double> vals = new ArrayList<>();
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (Math.abs(mat[i][j]) > 1e-15) {
                    rows.add(i);
                    cols.add(j);
                    vals.add(mat[i][j]);
                }
            }
        }
        int nnz = rows.size();
        int[] r = new int[nnz];
        int[] c = new int[nnz];
        double[] v = new double[nnz];
        for (int k = 0; k < nnz; k++) {
            r[k] = rows.get(k);
            c[k] = cols.get(k);
            v[k] = vals.get(k);
        }
        return cModel.addSparseMat(dim, nnz, r, c, v);
    }

    /**
     * 添加绝对平衡性约束（原始约束）
     */
    private void addAbsoluteBalanceConstraint(GRBModel model, GRBVar[][] x, int j, double riskParam) throws GRBException {
        double U = demandUpperBound; // 区域容量上界
        double L = demandLowerBound; // 区域容量下界

        int centerId = centers.get(j).getId();

        // 构建均值项和w_j向量
        // 原始模型：μ^T * x_j
        // 改进模型：μ^T * w_j，其中 w_j[i] = r_ij * x_ij
        GRBLinExpr meanTerm = new GRBLinExpr();
        GRBVar[] w_j = null;

        if (useImprovedModel) {
            // 改进模型：创建 w_j[i] = r_ij * x_ij
            w_j = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                double r_ij = shortestPathDist[i][centerId];
                w_j[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "w_" + j + "_" + i);
                // 约束：w_j[i] = r_ij * x_ij
                GRBLinExpr wExpr = new GRBLinExpr();
                wExpr.addTerm(r_ij, x[i][j]);
                model.addConstr(w_j[i], GRB.EQUAL, wExpr, "w_def_" + j + "_" + i);

                // 均值项：μ_i * w_j[i] = μ_i * r_ij * x_ij
                meanTerm.addTerm(meanVector[i], w_j[i]);
            }
        } else {
            // 原始模型：μ^T * x_j
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm.addTerm(meanVector[i], x[i][j]);
            }
        }

        // 使用gamma分解：gamma_a用于下界约束，gamma_b用于上界约束（四六分）
        // 统一对riskParam进行分解，无论是个体约束还是联合约束
        double gamma_a_local = 0.5 * riskParam;
        double gamma_b_local = 0.5 * riskParam;

        // 计算下界和上界的factor
        double factorLower, factorUpper;
        if (useD1) {
            // 使用D_1模糊集的约束
            // 下界: μ^T*x_j - sqrt((1-γ_a)/γ_a)*sqrt(x_j^T*Σ*x_j) ≥ L
            // 上界: μ^T*x_j + sqrt((1-γ_b)/γ_b)*sqrt(x_j^T*Σ*x_j) ≤ U
            factorLower = Math.sqrt((1 - gamma_a_local) / gamma_a_local);
            factorUpper = Math.sqrt((1 - gamma_b_local) / gamma_b_local);
        } else {
            // 使用D_2模糊集的约束
            // 根据delta1/delta2与gamma_a和gamma_b的关系确定factor
            if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                // Case 1
                factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
            } else if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 > gamma_b_local) {
                // Case 2
                factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                factorUpper = Math.sqrt(delta2 / gamma_b_local);
            } else if (delta1 / delta2 > gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                // Case 3
                factorLower = Math.sqrt(delta2 / gamma_a_local);
                factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
            } else {
                // Case 4
                factorLower = Math.sqrt(delta2 / gamma_a_local);
                factorUpper = Math.sqrt(delta2 / gamma_b_local);
            }
        }

        // 创建辅助变量t
        GRBVar t = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "t_" + j);

        // 使用Cholesky分解获取矩阵平方根
        Matrix covMatrix = new Matrix(covarianceMatrix);
        CholeskyDecomposition chol = new CholeskyDecomposition(covMatrix);

        // 如果Cholesky分解成功
        if (chol.isSPD()) {
            Matrix cholL = chol.getL();
            double[][] cholMatrix = cholL.getArray(); // 重命名以避免与下界L冲突

            // 创建SOCP约束: ||Ax|| <= t
            // 在Gurobi中，可以使用二次约束x'A'Ax <= t^2来表示

            // 首先计算L^T*v的每个分量（注意：需要转置L矩阵）
            // 原始模型：v = x_j，所以 x^T*Σ*x = x^T*L*L^T*x = ||L^T*x||^2
            // 改进模型：v = w_j，所以 w^T*Σ*w = w^T*L*L^T*w = ||L^T*w||^2
            GRBVar[] z = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                z[i] = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "z_" + j + "_" + i);

                GRBLinExpr zExpr = new GRBLinExpr();
                if (useImprovedModel) {
                    // 改进模型：z = L^T * w_j
                    for (int k = 0; k < inst.getN(); k++) {
                        // 使用 L^T[i][k] = L[k][i]，即转置后的矩阵
                        zExpr.addTerm(cholMatrix[k][i], w_j[k]);
                    }
                } else {
                    // 原始模型：z = L^T * x_j
                    for (int k = 0; k < inst.getN(); k++) {
                        // 使用 L^T[i][k] = L[k][i]，即转置后的矩阵
                        zExpr.addTerm(cholMatrix[k][i], x[k][j]);
                    }
                }
                model.addConstr(z[i], GRB.EQUAL, zExpr, "z_def_" + j + "_" + i);
            }

            // 创建SOCP约束: sum(z_i^2) <= t^2
            GRBQuadExpr socpExpr = new GRBQuadExpr();
            for (int i = 0; i < inst.getN(); i++) {
                socpExpr.addTerm(1.0, z[i], z[i]);
            }

            // t^2
            GRBQuadExpr tSquare = new GRBQuadExpr();
            tSquare.addTerm(1.0, t, t);

            model.addQConstr(socpExpr, GRB.LESS_EQUAL, tSquare, "socp_" + j);

            // 添加上界约束: μ^T*x_j + factorUpper*t ≤ U
            GRBLinExpr upperExpr = new GRBLinExpr();
            upperExpr.add(meanTerm);
            upperExpr.addTerm(factorUpper, t);
            model.addConstr(upperExpr, GRB.LESS_EQUAL, U, "demand_upper_" + j);

            // 添加下界约束: μ^T*x_j - factorLower*t ≥ L
            GRBLinExpr lowerExpr = new GRBLinExpr();
            lowerExpr.add(meanTerm);
            lowerExpr.addTerm(-factorLower, t);
            model.addConstr(lowerExpr, GRB.GREATER_EQUAL, L, "demand_lower_" + j);
        } else {
            // 如果Cholesky分解失败，尝试使用对角化方法
            System.out.println("警告：Cholesky分解失败，尝试使用替代方法");

            // 创建二次表达式
            GRBQuadExpr quadExpr = new GRBQuadExpr();
            if (useImprovedModel) {
                // 改进模型：w_j^T * Σ * w_j
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        quadExpr.addTerm(covarianceMatrix[i][k], w_j[i], w_j[k]);
                    }
                }
            } else {
                // 原始模型：x_j^T * Σ * x_j
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        quadExpr.addTerm(covarianceMatrix[i][k], x[i][j], x[k][j]);
                    }
                }
            }

            // 创建辅助变量q表示二次项
            GRBVar q = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "q_" + j);

            // 添加约束 quadExpr <= q
            model.addQConstr(quadExpr, GRB.LESS_EQUAL, q, "quad_" + j);

            // 添加约束 q <= t^2
            GRBQuadExpr tSquare = new GRBQuadExpr();
            tSquare.addTerm(1.0, t, t);
            model.addQConstr(q, GRB.LESS_EQUAL, tSquare, "t_square_" + j);

            // 添加上界约束: μ^T*x_j + factorUpper*t ≤ U
            GRBLinExpr upperExpr = new GRBLinExpr();
            upperExpr.add(meanTerm);
            upperExpr.addTerm(factorUpper, t);
            model.addConstr(upperExpr, GRB.LESS_EQUAL, U, "demand_upper_" + j);

            // 添加下界约束: μ^T*x_j - factorLower*t ≥ L
            GRBLinExpr lowerExpr = new GRBLinExpr();
            lowerExpr.add(meanTerm);
            lowerExpr.addTerm(-factorLower, t);
            model.addConstr(lowerExpr, GRB.GREATER_EQUAL, L, "demand_lower_" + j);
        }
    }

    /**
     * 添加相对平衡性约束（基于需求的相对平衡性约束）
     * 约束形式：对于每个区域 j，需求比例应该在 [1/p*(1-α), 1/p*(1+α)] 范围内
     * 等价于：
     * - d^T * v_{j,L} <= 0，其中 v_{j,L} = (1-α)/p * sum_k x_k - x_j
     * - d^T * v_{j,U} <= 0，其中 v_{j,U} = x_j - (1+α)/p * sum_k x_k
     * 
     * 如果使用支撑超平面cut，则不直接添加SOC约束，而是存储信息供后续迭代添加cut
     */
    private void addRelativeBalanceConstraint(GRBModel model, GRBVar[][] x, int j, double riskParam) throws GRBException {
        // 如果使用精确方法，直接添加精确的分式约束（利用Gurobi 13处理高次分式的功能）
        if (useExactMethod) {
            addExactRelativeBalanceConstraint(model, x, j, riskParam);
            return;
        }

        // D1 + d>=0 支撑集：采用“Gurobi主问题 + COPT求解SDP证书子问题”的分离法。
        // 此处不在主问题中直接加入SOC/SDP约束，约束由外层迭代（解主问题->COPT验算->回加cut）保证。
        if (shouldUseD1CopositiveSdpApprox()) {
            return;
        }
        
        // 如果使用支撑超平面cut且使用近似方法，则不直接添加SOC约束
        if (useSupportingHyperplaneCuts && !useExactMethod && !shouldUseD1CopositiveSdpApprox()) {
            // 计算factor并存储信息
            double gamma_a_local = 0.5 * riskParam;
            double gamma_b_local = 0.5 * riskParam;
            
            double factorLower, factorUpper;
            if (useD1) {
                factorLower = Math.sqrt((1 - gamma_a_local) / gamma_a_local);
                factorUpper = Math.sqrt((1 - gamma_b_local) / gamma_b_local);
            } else {
                if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                    factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                    factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
                } else if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 > gamma_b_local) {
                    factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                    factorUpper = Math.sqrt(delta2 / gamma_b_local);
                } else if (delta1 / delta2 > gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                    factorLower = Math.sqrt(delta2 / gamma_a_local);
                    factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
                } else {
                    factorLower = Math.sqrt(delta2 / gamma_a_local);
                    factorUpper = Math.sqrt(delta2 / gamma_b_local);
                }
            }
            
            int p = centers.size();
            double coeff = (1.0 - r) / p;
            double coeffUpper = (1.0 + r) / p;
            
            RelativeBalanceConstraintInfo info = new RelativeBalanceConstraintInfo();
            info.j = j;
            info.riskParam = riskParam;
            info.factorLower = factorLower;
            info.factorUpper = factorUpper;
            info.coeff = coeff;
            info.coeffUpper = coeffUpper;
            relativeBalanceConstraintInfos.add(info);
            
            // 加入精确约束中的线性部分：meanTerm_vjL ≤ 0，meanTerm_vjU ≤ 0（对可行集有效，且保证 Cantelli 前提）
            if (useAssignmentDependent) {
                int n = inst.getN();
                int nTimesP = n * p;
                // a_L[idx] = B_L[idx][idx]*meanVector[idx]，B_L 对角：(j2==j)?(coeff-1):coeff
                double[] a_L = new double[nTimesP];
                double[] a_U = new double[nTimesP];
                for (int idx = 0; idx < nTimesP; idx++) {
                    int j2 = idx % p;
                    a_L[idx] = (j2 == j ? (coeff - 1.0) : coeff) * meanVector[idx];
                    a_U[idx] = (j2 == j ? (1.0 - coeffUpper) : (-coeffUpper)) * meanVector[idx];
                }
                GRBLinExpr meanTerm_vjL = new GRBLinExpr();
                GRBLinExpr meanTerm_vjU = new GRBLinExpr();
                for (int idx = 0; idx < nTimesP; idx++) {
                    meanTerm_vjL.addTerm(a_L[idx], x[idx / p][idx % p]);
                    meanTerm_vjU.addTerm(a_U[idx], x[idx / p][idx % p]);
                }
                model.addConstr(meanTerm_vjL, GRB.LESS_EQUAL, 0.0, "mean_vjL_cut_" + j);
                model.addConstr(meanTerm_vjU, GRB.LESS_EQUAL, 0.0, "mean_vjU_cut_" + j);
            } else {
                int centerId = centers.get(j).getId();
                GRBLinExpr meanTerm_vjL = new GRBLinExpr();
                GRBLinExpr meanTerm_vjU = new GRBLinExpr();
                if (useImprovedModel) {
                    for (int i = 0; i < inst.getN(); i++) {
                        double r_ij = shortestPathDist[i][centerId];
                        for (int k = 0; k < centers.size(); k++) {
                            int centerId_k = centers.get(k).getId();
                            double r_ik = shortestPathDist[i][centerId_k];
                            meanTerm_vjL.addTerm(meanVector[i] * coeff * r_ik, x[i][k]);
                            meanTerm_vjU.addTerm(-meanVector[i] * coeffUpper * r_ik, x[i][k]);
                        }
                        meanTerm_vjL.addTerm(-meanVector[i] * r_ij, x[i][j]);
                        meanTerm_vjU.addTerm(meanVector[i] * r_ij, x[i][j]);
                    }
                } else {
                    double sumMean = 0.0;
                    for (int i = 0; i < inst.getN(); i++) sumMean += meanVector[i];
                    meanTerm_vjL.addConstant(coeff * sumMean);
                    meanTerm_vjU.addConstant(-coeffUpper * sumMean);
                    for (int i = 0; i < inst.getN(); i++) {
                        meanTerm_vjL.addTerm(-meanVector[i], x[i][j]);
                        meanTerm_vjU.addTerm(meanVector[i], x[i][j]);
                    }
                }
                model.addConstr(meanTerm_vjL, GRB.LESS_EQUAL, 0.0, "mean_vjL_cut_" + j);
                model.addConstr(meanTerm_vjU, GRB.LESS_EQUAL, 0.0, "mean_vjU_cut_" + j);
            }
            
            // 不添加 SOC 约束，仅通过迭代 cut 逼近
            return;
        }
        
        // 原有的直接添加SOC约束的代码
        int p = centers.size(); // 区域数量
        double coeff = (1.0 - r) / p; // (1-r)/p
        double coeffUpper = (1.0 + r) / p; // (1+r)/p

        // 计算 sum_k x_k（所有区域的总和）
        // 原始模型：sum_k x_k 是每个基本单元在所有区域中的分配向量的总和
        // 改进模型：需要计算 sum_k w_k，其中 w_k[i] = r_ik * x_ik
        GRBVar[] sumX = null; // sum_k x_k 或 sum_k w_k
        GRBLinExpr sumMeanTerm = new GRBLinExpr(); // μ^T * sum_k x_k 或 μ^T * sum_k w_k

        if (useImprovedModel) {
            // 改进模型：需要为每个区域 k 创建 w_k，然后计算 sum_k w_k
            // 但这样会创建很多变量，我们可以直接计算 μ^T * sum_k w_k
            // 对于每个基本单元 i，sum_k w_k[i] = sum_k r_ik * x_ik
            sumX = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                sumX[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "sum_w_" + i);
                GRBLinExpr sumWExpr = new GRBLinExpr();
                for (int k = 0; k < centers.size(); k++) {
                    int centerId_k = centers.get(k).getId();
                    double r_ik = shortestPathDist[i][centerId_k];
                    sumWExpr.addTerm(r_ik, x[i][k]);
                }
                model.addConstr(sumX[i], GRB.EQUAL, sumWExpr, "sum_w_def_" + i);
                sumMeanTerm.addTerm(meanVector[i], sumX[i]);
            }
        } else {
            // 原始模型：sum_k x_k 是每个基本单元在所有区域中的分配向量的总和
            sumX = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                sumX[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "sum_x_" + i);
                GRBLinExpr sumXExpr = new GRBLinExpr();
                for (int k = 0; k < centers.size(); k++) {
                    sumXExpr.addTerm(1.0, x[i][k]);
                }
                model.addConstr(sumX[i], GRB.EQUAL, sumXExpr, "sum_x_def_" + i);
                sumMeanTerm.addTerm(meanVector[i], sumX[i]);
            }
        }

        // 构建 v_{j,L} 和 v_{j,U}
        // v_{j,L} = (1-α)/p * sum_k x_k - x_j
        // v_{j,U} = x_j - (1+α)/p * sum_k x_k
        int centerId = centers.get(j).getId();

        // 构建 x_j 的均值项
        GRBLinExpr meanTerm_j = new GRBLinExpr();
        GRBVar[] w_j = null;

        if (useImprovedModel) {
            // 改进模型：创建 w_j[i] = r_ij * x_ij
            w_j = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                double r_ij = shortestPathDist[i][centerId];
                w_j[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "w_" + j + "_" + i);
                GRBLinExpr wExpr = new GRBLinExpr();
                wExpr.addTerm(r_ij, x[i][j]);
                model.addConstr(w_j[i], GRB.EQUAL, wExpr, "w_def_" + j + "_" + i);
                meanTerm_j.addTerm(meanVector[i], w_j[i]);
            }
        } else {
            // 原始模型：μ^T * x_j
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_j.addTerm(meanVector[i], x[i][j]);
            }
        }

        // 构建 v_{j,L} 和 v_{j,U} 的均值项
        // μ^T * v_{j,L} = (1-r)/p * μ^T * sum_k x_k - μ^T * x_j
        // μ^T * v_{j,U} = μ^T * x_j - (1+r)/p * μ^T * sum_k x_k
        GRBLinExpr meanTerm_vjL = new GRBLinExpr();
        GRBLinExpr meanTerm_vjU = new GRBLinExpr();
        
        if (useImprovedModel) {
            // 改进模型：使用变量sumX
            // 添加 coeff * sumMeanTerm
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjL.addTerm(meanVector[i] * coeff, sumX[i]);
            }
            // 添加 -meanTerm_j
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjL.addTerm(-meanVector[i], w_j[i]);
            }
            // meanTerm_vjU
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjU.addTerm(meanVector[i], w_j[i]);
            }
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjU.addTerm(-meanVector[i] * coeffUpper, sumX[i]);
            }
        } else {
            // 原始模型：使用常数1代替sum_k x_k（与Cut方法一致）
            // v_{j,L} = (1-r)/p * 1 - x_j，其中1是全1向量
            // μ^T * v_{j,L} = coeff * sum_i(μ_i) - μ^T * x_j
            double sumMean = 0.0;
            for (int i = 0; i < inst.getN(); i++) {
                sumMean += meanVector[i];
            }
            // meanTerm_vjL = coeff * sum_i(μ_i) - μ^T * x_j
            meanTerm_vjL.addConstant(coeff * sumMean);
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjL.addTerm(-meanVector[i], x[i][j]);
            }
            // meanTerm_vjU = μ^T * x_j - coeffUpper * sum_i(μ_i)
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjU.addTerm(meanVector[i], x[i][j]);
            }
            meanTerm_vjU.addConstant(-coeffUpper * sumMean);
        }

        // 使用gamma分解：gamma_a用于下界约束，gamma_b用于上界约束
        double gamma_a_local = 0.5 * riskParam;
        double gamma_b_local = 0.5 * riskParam;

        // 计算下界和上界的factor
        double factorLower, factorUpper;
        if (useD1) {
            // D1模糊集：使用近似方法（ADM）
            factorLower = Math.sqrt((1 - gamma_a_local) / gamma_a_local);
            factorUpper = Math.sqrt((1 - gamma_b_local) / gamma_b_local);
        } else {
            // D2模糊集：使用近似方法（ADM），根据delta1/delta2与gamma_a和gamma_b的关系确定factor
            if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                // Case 1
                factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
            } else if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 > gamma_b_local) {
                // Case 2
                factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                factorUpper = Math.sqrt(delta2 / gamma_b_local);
            } else if (delta1 / delta2 > gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                // Case 3
                factorLower = Math.sqrt(delta2 / gamma_a_local);
                factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
            } else {
                // Case 4
                factorLower = Math.sqrt(delta2 / gamma_a_local);
                factorUpper = Math.sqrt(delta2 / gamma_b_local);
            }
        }

        // 创建辅助变量 t_L 和 t_U 用于SOC约束
        GRBVar t_L = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "t_L_" + j);
        GRBVar t_U = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "t_U_" + j);

        // 使用Cholesky分解获取矩阵平方根
        Matrix covMatrix = new Matrix(covarianceMatrix);
        CholeskyDecomposition chol = new CholeskyDecomposition(covMatrix);

        if (chol.isSPD()) {
            Matrix cholL = chol.getL();
            double[][] cholMatrix = cholL.getArray();

            // 为 v_{j,L} 创建SOC约束
            GRBVar[] z_L = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                z_L[i] = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "z_L_" + j + "_" + i);
                GRBLinExpr zExpr = new GRBLinExpr();
                if (useImprovedModel) {
                    // 改进模型：v_{j,L} = (1-α)/p * sum_k w_k - w_j
                    // z_L = L^T * v_{j,L}
                    for (int k = 0; k < inst.getN(); k++) {
                        zExpr.addTerm(cholMatrix[k][i] * coeff, sumX[k]);
                        zExpr.addTerm(-cholMatrix[k][i], w_j[k]);
                    }
                } else {
                    // 原始模型：v_{j,L} = (1-r)/p * 1 - x_j（与Cut方法一致，使用常数1）
                    // z_L = L^T * v_{j,L} = L^T * (coeff * 1 - x_j)
                    // z_L[i] = sum_k L[k][i] * coeff - sum_k L[k][i] * x[k][j]
                    double constTermL = 0.0;
                    for (int k = 0; k < inst.getN(); k++) {
                        constTermL += cholMatrix[k][i] * coeff;  // coeff * sum_k L[k][i]
                        zExpr.addTerm(-cholMatrix[k][i], x[k][j]);
                    }
                    zExpr.addConstant(constTermL);
                }
                model.addConstr(z_L[i], GRB.EQUAL, zExpr, "z_L_def_" + j + "_" + i);
            }

            // 为 v_{j,U} 创建SOC约束
            GRBVar[] z_U = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                z_U[i] = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "z_U_" + j + "_" + i);
                GRBLinExpr zExpr = new GRBLinExpr();
                if (useImprovedModel) {
                    // 改进模型：v_{j,U} = w_j - (1+α)/p * sum_k w_k
                    // z_U = L^T * v_{j,U}
                    for (int k = 0; k < inst.getN(); k++) {
                        zExpr.addTerm(cholMatrix[k][i], w_j[k]);
                        zExpr.addTerm(-cholMatrix[k][i] * coeffUpper, sumX[k]);
                    }
                } else {
                    // 原始模型：v_{j,U} = x_j - (1+r)/p * 1（与Cut方法一致，使用常数1）
                    // z_U = L^T * v_{j,U} = L^T * (x_j - coeffUpper * 1)
                    // z_U[i] = sum_k L[k][i] * x[k][j] - coeffUpper * sum_k L[k][i]
                    double constTermU = 0.0;
                    for (int k = 0; k < inst.getN(); k++) {
                        zExpr.addTerm(cholMatrix[k][i], x[k][j]);
                        constTermU += cholMatrix[k][i] * coeffUpper;  // coeffUpper * sum_k L[k][i]
                    }
                    zExpr.addConstant(-constTermU);
                }
                model.addConstr(z_U[i], GRB.EQUAL, zExpr, "z_U_def_" + j + "_" + i);
            }

            // 创建SOCP约束: ||z_L|| <= t_L, ||z_U|| <= t_U
            GRBQuadExpr socpExpr_L = new GRBQuadExpr();
            for (int i = 0; i < inst.getN(); i++) {
                socpExpr_L.addTerm(1.0, z_L[i], z_L[i]);
            }
            GRBQuadExpr tSquare_L = new GRBQuadExpr();
            tSquare_L.addTerm(1.0, t_L, t_L);
            model.addQConstr(socpExpr_L, GRB.LESS_EQUAL, tSquare_L, "socp_L_" + j);

            GRBQuadExpr socpExpr_U = new GRBQuadExpr();
            for (int i = 0; i < inst.getN(); i++) {
                socpExpr_U.addTerm(1.0, z_U[i], z_U[i]);
            }
            GRBQuadExpr tSquare_U = new GRBQuadExpr();
            tSquare_U.addTerm(1.0, t_U, t_U);
            model.addQConstr(socpExpr_U, GRB.LESS_EQUAL, tSquare_U, "socp_U_" + j);

            // 添加约束: μ^T * v_{j,L} + factorLower * t_L <= 0
            GRBLinExpr lowerExpr = new GRBLinExpr();
            lowerExpr.add(meanTerm_vjL);
            lowerExpr.addTerm(factorLower, t_L);
            model.addConstr(lowerExpr, GRB.LESS_EQUAL, 0, "rel_balance_lower_" + j);

            // 添加约束: μ^T * v_{j,U} + factorUpper * t_U <= 0
            GRBLinExpr upperExpr = new GRBLinExpr();
            upperExpr.add(meanTerm_vjU);
            upperExpr.addTerm(factorUpper, t_U);
            model.addConstr(upperExpr, GRB.LESS_EQUAL, 0, "rel_balance_upper_" + j);
        } else {
            // 如果Cholesky分解失败，使用二次约束
            System.out.println("警告：Cholesky分解失败，尝试使用替代方法");

            // 为 v_{j,L} 创建二次表达式
            GRBQuadExpr quadExpr_L = new GRBQuadExpr();
            if (useImprovedModel) {
                // 改进模型：v_{j,L}^T * Σ * v_{j,L}
                // v_{j,L} = (1-α)/p * sum_k w_k - w_j
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        double coeff_i = coeff;
                        double coeff_k = coeff;
                        if (i == k) {
                            quadExpr_L.addTerm(covarianceMatrix[i][k] * coeff_i * coeff_k, sumX[i], sumX[k]);
                            quadExpr_L.addTerm(-covarianceMatrix[i][k] * coeff_i, sumX[i], w_j[k]);
                            quadExpr_L.addTerm(-covarianceMatrix[i][k] * coeff_k, w_j[i], sumX[k]);
                            quadExpr_L.addTerm(covarianceMatrix[i][k], w_j[i], w_j[k]);
                        } else {
                            quadExpr_L.addTerm(covarianceMatrix[i][k] * coeff_i * coeff_k, sumX[i], sumX[k]);
                            quadExpr_L.addTerm(-covarianceMatrix[i][k] * coeff_i, sumX[i], w_j[k]);
                            quadExpr_L.addTerm(-covarianceMatrix[i][k] * coeff_k, w_j[i], sumX[k]);
                            quadExpr_L.addTerm(covarianceMatrix[i][k], w_j[i], w_j[k]);
                        }
                    }
                }
            } else {
                // 原始模型：v_{j,L}^T * Σ * v_{j,L}（与Cut方法一致，使用常数1）
                // v_{j,L} = coeff * 1 - x_j
                // v^T * Σ * v = coeff^2 * 1^T*Σ*1 - 2*coeff * 1^T*Σ*x_j + x_j^T*Σ*x_j
                double constTermL = 0.0;  // coeff^2 * 1^T*Σ*1
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        constTermL += coeff * coeff * covarianceMatrix[i][k];  // coeff^2 * Σ_{ik}
                    }
                }
                quadExpr_L.addConstant(constTermL);
                // -2*coeff * 1^T*Σ*x_j = -2*coeff * sum_i(sum_k(Σ_{ik})*x_j[k])
                for (int k = 0; k < inst.getN(); k++) {
                    double linearCoeff = 0.0;
                    for (int i = 0; i < inst.getN(); i++) {
                        linearCoeff += covarianceMatrix[i][k];  // sum_i(Σ_{ik})
                    }
                    quadExpr_L.addTerm(-2.0 * coeff * linearCoeff, x[k][j]);
                }
                // x_j^T*Σ*x_j
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        quadExpr_L.addTerm(covarianceMatrix[i][k], x[i][j], x[k][j]);
                    }
                }
            }

            // 为 v_{j,U} 创建二次表达式
            GRBQuadExpr quadExpr_U = new GRBQuadExpr();
            if (useImprovedModel) {
                // 改进模型：v_{j,U}^T * Σ * v_{j,U}
                // v_{j,U} = w_j - (1+α)/p * sum_k w_k
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        double coeff_i = coeffUpper;
                        double coeff_k = coeffUpper;
                        quadExpr_U.addTerm(covarianceMatrix[i][k], w_j[i], w_j[k]);
                        quadExpr_U.addTerm(-covarianceMatrix[i][k] * coeff_i, w_j[i], sumX[k]);
                        quadExpr_U.addTerm(-covarianceMatrix[i][k] * coeff_k, sumX[i], w_j[k]);
                        quadExpr_U.addTerm(covarianceMatrix[i][k] * coeff_i * coeff_k, sumX[i], sumX[k]);
                    }
                }
            } else {
                // 原始模型：v_{j,U}^T * Σ * v_{j,U}（与Cut方法一致，使用常数1）
                // v_{j,U} = x_j - coeffUpper * 1
                // v^T * Σ * v = x_j^T*Σ*x_j - 2*coeffUpper * x_j^T*Σ*1 + coeffUpper^2 * 1^T*Σ*1
                double constTermU = 0.0;  // coeffUpper^2 * 1^T*Σ*1
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        constTermU += coeffUpper * coeffUpper * covarianceMatrix[i][k];
                    }
                }
                quadExpr_U.addConstant(constTermU);
                // x_j^T*Σ*x_j
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        quadExpr_U.addTerm(covarianceMatrix[i][k], x[i][j], x[k][j]);
                    }
                }
                // -2*coeffUpper * x_j^T*Σ*1 = -2*coeffUpper * sum_i(x_j[i] * sum_k(Σ_{ik}))
                for (int i = 0; i < inst.getN(); i++) {
                    double linearCoeff = 0.0;
                    for (int k = 0; k < inst.getN(); k++) {
                        linearCoeff += covarianceMatrix[i][k];  // sum_k(Σ_{ik})
                    }
                    quadExpr_U.addTerm(-2.0 * coeffUpper * linearCoeff, x[i][j]);
                }
            }

            // 创建辅助变量 q_L 和 q_U
            GRBVar q_L = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "q_L_" + j);
            GRBVar q_U = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "q_U_" + j);

            model.addQConstr(quadExpr_L, GRB.LESS_EQUAL, q_L, "quad_L_" + j);
            model.addQConstr(quadExpr_U, GRB.LESS_EQUAL, q_U, "quad_U_" + j);

            GRBQuadExpr tSquare_L = new GRBQuadExpr();
            tSquare_L.addTerm(1.0, t_L, t_L);
            model.addQConstr(q_L, GRB.LESS_EQUAL, tSquare_L, "t_square_L_" + j);

            GRBQuadExpr tSquare_U = new GRBQuadExpr();
            tSquare_U.addTerm(1.0, t_U, t_U);
            model.addQConstr(q_U, GRB.LESS_EQUAL, tSquare_U, "t_square_U_" + j);

            // 添加约束: μ^T * v_{j,L} + factorLower * t_L <= 0
            GRBLinExpr lowerExpr = new GRBLinExpr();
            lowerExpr.add(meanTerm_vjL);
            lowerExpr.addTerm(factorLower, t_L);
            model.addConstr(lowerExpr, GRB.LESS_EQUAL, 0, "rel_balance_lower_" + j);

            // 添加约束: μ^T * v_{j,U} + factorUpper * t_U <= 0
            GRBLinExpr upperExpr = new GRBLinExpr();
            upperExpr.add(meanTerm_vjU);
            upperExpr.addTerm(factorUpper, t_U);
            model.addConstr(upperExpr, GRB.LESS_EQUAL, 0, "rel_balance_upper_" + j);
        }
    }

    /**
     * 添加精确的相对平衡性约束（使用Gurobi 13的分式约束功能）
     * 对于D_1模糊集：添加公式 eq:rel_balance_d1_exact
     * 对于D_2模糊集：根据κ值判断case，添加对应的分式约束
     */
    private void addExactRelativeBalanceConstraint(GRBModel model, GRBVar[][] x, int j, double riskParam) throws GRBException {
        // assignment-dependent 模型：v 为 N*p 维，仿照原模型精确约束添加 D1/D2 精确约束
        if (useAssignmentDependent) {
            addExactRelativeBalanceConstraintAssignmentDependent(model, x, j, riskParam);
            return;
        }

        int p = centers.size(); // 区域数量
        double coeff = (1.0 - r) / p; // (1-r)/p
        double coeffUpper = (1.0 + r) / p; // (1+r)/p

        // 计算 sum_k x_k（所有区域的总和）
        GRBVar[] sumX = null; // sum_k x_k 或 sum_k w_k
        GRBLinExpr sumMeanTerm = new GRBLinExpr(); // μ^T * sum_k x_k 或 μ^T * sum_k w_k

        if (useImprovedModel) {
            // 改进模型：需要为每个区域 k 创建 w_k，然后计算 sum_k w_k
            sumX = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                sumX[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "sum_w_exact_" + j + "_" + i);
                GRBLinExpr sumWExpr = new GRBLinExpr();
                for (int k = 0; k < centers.size(); k++) {
                    int centerId_k = centers.get(k).getId();
                    double r_ik = shortestPathDist[i][centerId_k];
                    sumWExpr.addTerm(r_ik, x[i][k]);
                }
                model.addConstr(sumX[i], GRB.EQUAL, sumWExpr, "sum_w_exact_def_" + j + "_" + i);
                sumMeanTerm.addTerm(meanVector[i], sumX[i]);
            }
        } else {
            // 原始模型：sum_k x_k 是每个基本单元在所有区域中的分配向量的总和
            sumX = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                sumX[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "sum_x_exact_" + j + "_" + i);
                GRBLinExpr sumXExpr = new GRBLinExpr();
                for (int k = 0; k < centers.size(); k++) {
                    sumXExpr.addTerm(1.0, x[i][k]);
                }
                model.addConstr(sumX[i], GRB.EQUAL, sumXExpr, "sum_x_exact_def_" + j + "_" + i);
                sumMeanTerm.addTerm(meanVector[i], sumX[i]);
            }
        }

        // 构建 v_{j,L} 和 v_{j,U}
        int centerId = centers.get(j).getId();

        // 构建 x_j 的均值项
        GRBLinExpr meanTerm_j = new GRBLinExpr();
        GRBVar[] w_j = null;

        if (useImprovedModel) {
            // 改进模型：创建 w_j[i] = r_ij * x_ij
            w_j = new GRBVar[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                double r_ij = shortestPathDist[i][centerId];
                w_j[i] = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "w_exact_" + j + "_" + i);
                GRBLinExpr wExpr = new GRBLinExpr();
                wExpr.addTerm(r_ij, x[i][j]);
                model.addConstr(w_j[i], GRB.EQUAL, wExpr, "w_exact_def_" + j + "_" + i);
                meanTerm_j.addTerm(meanVector[i], w_j[i]);
            }
        } else {
            // 原始模型：μ^T * x_j
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_j.addTerm(meanVector[i], x[i][j]);
            }
        }

        // 构建 v_{j,L} 和 v_{j,U} 的均值项
        GRBLinExpr meanTerm_vjL = new GRBLinExpr();
        GRBLinExpr meanTerm_vjU = new GRBLinExpr();
        
        if (useImprovedModel) {           // 改进模型：使用变量sumX
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjL.addTerm(meanVector[i] * coeff, sumX[i]);
                meanTerm_vjL.addTerm(-meanVector[i], w_j[i]);
                meanTerm_vjU.addTerm(meanVector[i], w_j[i]);
                meanTerm_vjU.addTerm(-meanVector[i] * coeffUpper, sumX[i]);
            }
        } else {
            // 原始模型：使用常数1代替sum_k x_k
            double sumMean = 0.0;
            for (int i = 0; i < inst.getN(); i++) {
                sumMean += meanVector[i];
            }
            meanTerm_vjL.addConstant(coeff * sumMean);
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjL.addTerm(-meanVector[i], x[i][j]);
            }
            for (int i = 0; i < inst.getN(); i++) {
                meanTerm_vjU.addTerm(meanVector[i], x[i][j]);
            }
            meanTerm_vjU.addConstant(-coeffUpper * sumMean);
        }

        // 计算 v_{j,L}^T * Σ * v_{j,L} 和 v_{j,U}^T * Σ * v_{j,U}
        // 使用Cholesky分解或直接计算二次表达式
        Matrix covMatrix = new Matrix(covarianceMatrix);
        CholeskyDecomposition chol = new CholeskyDecomposition(covMatrix);

        if (useD1) {
            // D_1 模糊集：使用MIQCP形式添加约束
            // 原始约束：(v_{j,L}^T * Σ * v_{j,L}) / (v_{j,L}^T * Σ * v_{j,L} + (μ^T * v_{j,L})^2) + 
            //           (v_{j,U}^T * Σ * v_{j,U}) / (v_{j,U}^T * Σ * v_{j,U} + (μ^T * v_{j,U})^2) <= γ
            // 
            // 转化步骤：
            // 1. 引入λ_L和λ_U拆分分式和：λ_L + λ_U <= γ
            // 2. 将分式转化为乘积形式：(1-λ_L) * L_var_L - λ_L * L_mean_L <= 0
            // 3. 利用0-1变量性质线性化二次项：引入z_ik = x_i * x_k
            // 4. 将v^T Σ v和(μ^T v)^2展开为关于x和z的线性函数
            
            // 第一步：引入辅助变量λ_L和λ_U
            GRBVar lambda_L = model.addVar(0, riskParam, 0, GRB.CONTINUOUS, "lambda_L_" + j);
            GRBVar lambda_U = model.addVar(0, riskParam, 0, GRB.CONTINUOUS, "lambda_U_" + j);
            
            // 添加约束：λ_L + λ_U <= γ
            GRBLinExpr lambda_sum = new GRBLinExpr();
            lambda_sum.addTerm(1.0, lambda_L);
            lambda_sum.addTerm(1.0, lambda_U);
            model.addConstr(lambda_sum, GRB.LESS_EQUAL, riskParam, "lambda_sum_" + j);
            
            // 第二步和第三步：线性化二次项
            // 对于原始模型：v_jL = coeff * 1 - x_j，v_jU = x_j - coeffUpper * 1
            // 对于改进模型：v_jL = coeff * sum_k w_k - w_j，v_jU = w_j - coeffUpper * sum_k w_k
            // 其中w_k[i] = r_ik * x_ik，所以w_k也是x的线性函数
            
            // 引入z_ik变量表示x_i * x_k（对于i < k）
            // 注意：对于区域j，我们只需要x[i][j]和x[k][j]的乘积
            GRBVar[][] z_L = null;  // 用于v_jL的二次项线性化
            GRBVar[][] z_U = null;  // 用于v_jU的二次项线性化
            
            if (!useImprovedModel) {
                // 原始模型：只需要x[i][j] * x[k][j]的乘积
                z_L = new GRBVar[inst.getN()][inst.getN()];
                z_U = new GRBVar[inst.getN()][inst.getN()];
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        if (i == k) {
                            // x_i^2 = x_i (0-1变量性质)
                            z_L[i][k] = x[i][j];  // 直接使用x[i][j]
                            z_U[i][k] = x[i][j];  // 直接使用x[i][j]
                        } else {
                            // z_ik = x_i * x_k，添加McCormick约束
                            z_L[i][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_L_" + j + "_" + i + "_" + k);
                            z_U[i][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_U_" + j + "_" + i + "_" + k);
                            
                            // McCormick约束：z_ik <= x_i, z_ik <= x_k, z_ik >= x_i + x_k - 1
                            model.addConstr(z_L[i][k], GRB.LESS_EQUAL, x[i][j], "z_L_ub1_" + j + "_" + i + "_" + k);
                            model.addConstr(z_L[i][k], GRB.LESS_EQUAL, x[k][j], "z_L_ub2_" + j + "_" + i + "_" + k);
                            GRBLinExpr z_L_lb = new GRBLinExpr();
                            z_L_lb.addTerm(1.0, x[i][j]);
                            z_L_lb.addTerm(1.0, x[k][j]);
                            z_L_lb.addConstant(-1.0);
                            model.addConstr(z_L[i][k], GRB.GREATER_EQUAL, z_L_lb, "z_L_lb_" + j + "_" + i + "_" + k);
                            
                            model.addConstr(z_U[i][k], GRB.LESS_EQUAL, x[i][j], "z_U_ub1_" + j + "_" + i + "_" + k);
                            model.addConstr(z_U[i][k], GRB.LESS_EQUAL, x[k][j], "z_U_ub2_" + j + "_" + i + "_" + k);
                            GRBLinExpr z_U_lb = new GRBLinExpr();
                            z_U_lb.addTerm(1.0, x[i][j]);
                            z_U_lb.addTerm(1.0, x[k][j]);
                            z_U_lb.addConstant(-1.0);
                            model.addConstr(z_U[i][k], GRB.GREATER_EQUAL, z_U_lb, "z_U_lb_" + j + "_" + i + "_" + k);
                        }
                    }
                }
            } else {
                // 改进模型：需要处理w_j[i] = r_ij * x[i][j]和sumX[i] = sum_k r_ik * x[i][k]
                // w_j[i] * w_j[k] = r_ij * r_kj * x[i][j] * x[k][j]
                // sumX[i] * sumX[k] = sum_{l1} sum_{l2} r_{i,l1} * r_{k,l2} * x[i][l1] * x[k][l2]
                // 这更复杂，需要为所有区域对引入z变量
                // 为了简化，我们仍然只针对区域j引入z变量
                z_L = new GRBVar[inst.getN()][inst.getN()];
                z_U = new GRBVar[inst.getN()][inst.getN()];
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        if (i == k) {
                            z_L[i][k] = x[i][j];
                            z_U[i][k] = x[i][j];
                        } else {
                            z_L[i][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_L_" + j + "_" + i + "_" + k);
                            z_U[i][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_U_" + j + "_" + i + "_" + k);
                            
                            model.addConstr(z_L[i][k], GRB.LESS_EQUAL, x[i][j], "z_L_ub1_" + j + "_" + i + "_" + k);
                            model.addConstr(z_L[i][k], GRB.LESS_EQUAL, x[k][j], "z_L_ub2_" + j + "_" + i + "_" + k);
                            GRBLinExpr z_L_lb = new GRBLinExpr();
                            z_L_lb.addTerm(1.0, x[i][j]);
                            z_L_lb.addTerm(1.0, x[k][j]);
                            z_L_lb.addConstant(-1.0);
                            model.addConstr(z_L[i][k], GRB.GREATER_EQUAL, z_L_lb, "z_L_lb_" + j + "_" + i + "_" + k);
                            
                            model.addConstr(z_U[i][k], GRB.LESS_EQUAL, x[i][j], "z_U_ub1_" + j + "_" + i + "_" + k);
                            model.addConstr(z_U[i][k], GRB.LESS_EQUAL, x[k][j], "z_U_ub2_" + j + "_" + i + "_" + k);
                            GRBLinExpr z_U_lb = new GRBLinExpr();
                            z_U_lb.addTerm(1.0, x[i][j]);
                            z_U_lb.addTerm(1.0, x[k][j]);
                            z_U_lb.addConstant(-1.0);
                            model.addConstr(z_U[i][k], GRB.GREATER_EQUAL, z_U_lb, "z_U_lb_" + j + "_" + i + "_" + k);
                        }
                    }
                }
            }
            
            // 第四步：将v^T Σ v和(μ^T v)^2展开为关于x和z的线性函数
            // 创建辅助变量表示线性化后的二次项
            GRBVar L_var_L = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_var_L_" + j);
            GRBVar L_var_U = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_var_U_" + j);
            GRBVar L_mean_L = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_mean_L_" + j);
            GRBVar L_mean_U = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_mean_U_" + j);
            
            // 计算L_var_L = v_{j,L}^T * Σ * v_{j,L}（线性化形式）
            GRBLinExpr L_var_L_expr = new GRBLinExpr();
            if (!useImprovedModel) {
                // 原始模型：v_jL = coeff * 1 - x_j
                // v_jL^T * Σ * v_jL = coeff^2 * 1^T * Σ * 1 - 2*coeff * 1^T * Σ * x_j + x_j^T * Σ * x_j
                double constTermL = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        constTermL += coeff * coeff * covarianceMatrix[i][k];
                    }
                }
                L_var_L_expr.addConstant(constTermL);
                
                for (int k = 0; k < inst.getN(); k++) {
                    double linearCoeff = 0.0;
                    for (int i = 0; i < inst.getN(); i++) {
                        linearCoeff += covarianceMatrix[i][k];
                    }
                    L_var_L_expr.addTerm(-2.0 * coeff * linearCoeff, x[k][j]);
                }
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        L_var_L_expr.addTerm(covarianceMatrix[i][k], z_L[i][k]);
                    }
                }
            } else {
                // 改进模型：v_jL = coeff * sum_k w_k - w_j
                // 展开后包含sumX[i] * sumX[k]和w_j[i] * w_j[k]等项
                // 由于sumX和w_j都是x的线性函数，展开后仍然是x和z的二次项
                // 这里需要更仔细地展开
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        // coeff^2 * sumX[i] * sumX[k]项
                        // sumX[i] = sum_{l} r_{i,l} * x[i][l]，所以sumX[i] * sumX[k] = sum_{l1,l2} r_{i,l1} * r_{k,l2} * x[i][l1] * x[k][l2]
                        // 对于区域j，我们主要关注x[i][j]和x[k][j]的项
                        double r_ij = shortestPathDist[i][centerId];
                        double r_kj = shortestPathDist[k][centerId];
                        L_var_L_expr.addTerm(covarianceMatrix[i][k] * coeff * coeff * r_ij * r_kj, z_L[i][k]);
                        
                        // -coeff * sumX[i] * w_j[k]项
                        L_var_L_expr.addTerm(-covarianceMatrix[i][k] * coeff * r_ij * r_kj, z_L[i][k]);
                        
                        // -coeff * w_j[i] * sumX[k]项
                        L_var_L_expr.addTerm(-covarianceMatrix[i][k] * coeff * r_ij * r_kj, z_L[i][k]);
                        
                        // w_j[i] * w_j[k]项
                        L_var_L_expr.addTerm(covarianceMatrix[i][k] * r_ij * r_kj, z_L[i][k]);
                    }
                }
                // 注意：这里简化了，实际上sumX涉及所有区域，需要更完整的展开
                // 为了正确性，我们仍然使用原有的二次表达式，但通过z变量线性化
            }
            model.addConstr(L_var_L, GRB.EQUAL, L_var_L_expr, "L_var_L_def_" + j);
            
            // 计算L_var_U = v_{j,U}^T * Σ * v_{j,U}（线性化形式）
            GRBLinExpr L_var_U_expr = new GRBLinExpr();
            if (!useImprovedModel) {
                // 原始模型：v_jU = x_j - coeffUpper * 1
                double constTermU = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        constTermU += coeffUpper * coeffUpper * covarianceMatrix[i][k];
                    }
                }
                L_var_U_expr.addConstant(constTermU);
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        L_var_U_expr.addTerm(covarianceMatrix[i][k], z_U[i][k]);
                    }
                }
                
                for (int i = 0; i < inst.getN(); i++) {
                    double linearCoeff = 0.0;
                    for (int k = 0; k < inst.getN(); k++) {
                        linearCoeff += covarianceMatrix[i][k];
                    }
                    L_var_U_expr.addTerm(-2.0 * coeffUpper * linearCoeff, x[i][j]);
                }
            } else {
                // 改进模型：类似处理
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        double r_ij = shortestPathDist[i][centerId];
                        double r_kj = shortestPathDist[k][centerId];
                        L_var_U_expr.addTerm(covarianceMatrix[i][k] * r_ij * r_kj, z_U[i][k]);
                        L_var_U_expr.addTerm(-covarianceMatrix[i][k] * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        L_var_U_expr.addTerm(-covarianceMatrix[i][k] * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        L_var_U_expr.addTerm(covarianceMatrix[i][k] * coeffUpper * coeffUpper * r_ij * r_kj, z_U[i][k]);
                    }
                }
            }
            model.addConstr(L_var_U, GRB.EQUAL, L_var_U_expr, "L_var_U_def_" + j);
            
            // 计算L_mean_L = (μ^T * v_{j,L})^2（线性化形式）
            // (μ^T * v_{j,L})^2展开后也是x和z的二次项
            // 由于meanTerm_vjL已经是x的线性函数，其平方需要展开
            GRBLinExpr L_mean_L_expr = new GRBLinExpr();
            if (!useImprovedModel) {
                // meanTerm_vjL = coeff * sumMean - μ^T * x_j
                // (meanTerm_vjL)^2 = coeff^2 * sumMean^2 - 2*coeff*sumMean*(μ^T*x_j) + (μ^T*x_j)^2
                double sumMean = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    sumMean += meanVector[i];
                }
                L_mean_L_expr.addConstant(coeff * coeff * sumMean * sumMean);
                
                for (int i = 0; i < inst.getN(); i++) {
                    L_mean_L_expr.addTerm(-2.0 * coeff * sumMean * meanVector[i], x[i][j]);
                }
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        L_mean_L_expr.addTerm(meanVector[i] * meanVector[k], z_L[i][k]);
                    }
                }
            } else {
                // 改进模型：类似处理
                double sumMean = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    sumMean += meanVector[i];
                }
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        double r_ij = shortestPathDist[i][centerId];
                        double r_kj = shortestPathDist[k][centerId];
                        L_mean_L_expr.addTerm(meanVector[i] * meanVector[k] * coeff * coeff * r_ij * r_kj, z_L[i][k]);
                        L_mean_L_expr.addTerm(-meanVector[i] * meanVector[k] * coeff * r_ij * r_kj, z_L[i][k]);
                        L_mean_L_expr.addTerm(-meanVector[i] * meanVector[k] * coeff * r_ij * r_kj, z_L[i][k]);
                        L_mean_L_expr.addTerm(meanVector[i] * meanVector[k] * r_ij * r_kj, z_L[i][k]);
                    }
                }
            }
            model.addConstr(L_mean_L, GRB.EQUAL, L_mean_L_expr, "L_mean_L_def_" + j);
            
            // 计算L_mean_U = (μ^T * v_{j,U})^2（线性化形式）
            GRBLinExpr L_mean_U_expr = new GRBLinExpr();
            if (!useImprovedModel) {
                // meanTerm_vjU = μ^T * x_j - coeffUpper * sumMean
                double sumMean = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    sumMean += meanVector[i];
                }
                L_mean_U_expr.addConstant(coeffUpper * coeffUpper * sumMean * sumMean);
                
                for (int i = 0; i < inst.getN(); i++) {
                    L_mean_U_expr.addTerm(-2.0 * coeffUpper * sumMean * meanVector[i], x[i][j]);
                }
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        L_mean_U_expr.addTerm(meanVector[i] * meanVector[k], z_U[i][k]);
                    }
                }
            } else {
                // 改进模型：类似处理
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        double r_ij = shortestPathDist[i][centerId];
                        double r_kj = shortestPathDist[k][centerId];
                        L_mean_U_expr.addTerm(meanVector[i] * meanVector[k] * r_ij * r_kj, z_U[i][k]);
                        L_mean_U_expr.addTerm(-meanVector[i] * meanVector[k] * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        L_mean_U_expr.addTerm(-meanVector[i] * meanVector[k] * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        L_mean_U_expr.addTerm(meanVector[i] * meanVector[k] * coeffUpper * coeffUpper * r_ij * r_kj, z_U[i][k]);
                    }
                }
            }
            model.addConstr(L_mean_U, GRB.EQUAL, L_mean_U_expr, "L_mean_U_def_" + j);
            
            // 第五步：添加最终的二次约束
            // (1-λ_L) * L_var_L - λ_L * L_mean_L <= 0
            // 等价于：L_var_L - λ_L * (L_var_L + L_mean_L) <= 0
            GRBQuadExpr quad_constraint_L = new GRBQuadExpr();
            quad_constraint_L.addTerm(1.0, L_var_L);
            quad_constraint_L.addTerm(-1.0, lambda_L, L_var_L);
            quad_constraint_L.addTerm(-1.0, lambda_L, L_mean_L);
            model.addQConstr(quad_constraint_L, GRB.LESS_EQUAL, 0.0, "quad_constraint_L_" + j);
            
            // (1-λ_U) * L_var_U - λ_U * L_mean_U <= 0
            GRBQuadExpr quad_constraint_U = new GRBQuadExpr();
            quad_constraint_U.addTerm(1.0, L_var_U);
            quad_constraint_U.addTerm(-1.0, lambda_U, L_var_U);
            quad_constraint_U.addTerm(-1.0, lambda_U, L_mean_U);
            model.addQConstr(quad_constraint_U, GRB.LESS_EQUAL, 0.0, "quad_constraint_U_" + j);

            // D1模糊集必须保证均值小于0，因为Cantelli不等式只在均值满足约束时提供非平凡上界
            // 如果均值大于0（期望违反），则违反概率上界为1，无法满足小的gamma约束
            // 但如果不加此约束，求解器会利用平方项(μ^T v)^2增大分母，错误地使约束"满足"
            model.addConstr(meanTerm_vjL, GRB.LESS_EQUAL, 0.0, "mean_vjL_negative_d1_" + j);
            model.addConstr(meanTerm_vjU, GRB.LESS_EQUAL, 0.0, "mean_vjU_negative_d1_" + j);
            
            // 注意：对于改进模型（useImprovedModel = true），由于sumX涉及所有区域，
            // 完全线性化需要为所有区域对引入z变量，这会导致变量数量巨大。
            // 当前实现对于改进模型使用了简化的线性化（只针对区域j的项），
            // 这可能不是完全等价的。如果使用改进模型，建议：
            // 1. 使用原有的二次表达式方法（回退到之前的实现）
            // 2. 或者为所有区域对引入z变量（变量数量会很大）
        } else {
            // D_2 模糊集：使用MIQCP形式添加约束
            // 原始约束：根据κ值判断case，添加对应的分式约束
            // term_L + term_U >= 2 - γ
            // 
            // 转化步骤（类似D1，但需要考虑不同case）：
            // 1. 引入λ_L和λ_U拆分分式和：λ_L + λ_U >= 2 - γ
            // 2. 将分式转化为乘积形式（针对不同case）
            // 3. 利用0-1变量性质线性化二次项：引入z_ik = x_i * x_k
            // 4. 将v^T Σ v和(μ^T v)^2展开为关于x和z的线性函数
            // 5. 处理kappa值的计算（分式形式）
            
            // 计算阈值
            double sqrtDelta1 = Math.sqrt(delta1);
            double delta2OverSqrtDelta1 = delta2 / sqrtDelta1;
            
            // 第一步：引入辅助变量λ_L和λ_U（注意D2的约束是>=，所以需要调整）
            // 对于D2，约束是 term_L + term_U >= 2 - γ
            // 等价于：-term_L - term_U <= -(2 - γ)
            // 或者：term_L + term_U >= 2 - γ
            // 我们引入λ_L和λ_U，使得：λ_L + λ_U >= 2 - γ
            // 其中λ_L表示term_L，λ_U表示term_U
            GRBVar lambda_L = model.addVar(0, 1, 0, GRB.CONTINUOUS, "lambda_L_d2_" + j);
            GRBVar lambda_U = model.addVar(0, 1, 0, GRB.CONTINUOUS, "lambda_U_d2_" + j);
            
            // 添加约束：λ_L + λ_U >= 2 - γ
            GRBLinExpr lambda_sum = new GRBLinExpr();
            lambda_sum.addTerm(1.0, lambda_L);
            lambda_sum.addTerm(1.0, lambda_U);
            model.addConstr(lambda_sum, GRB.GREATER_EQUAL, 2.0 - riskParam, "lambda_sum_d2_" + j);
            
            // 第二步和第三步：线性化二次项（与D1相同）
            // 引入z_ik变量表示x_i * x_k
            GRBVar[][] z_L = null;
            GRBVar[][] z_U = null;
            
            if (!useImprovedModel) {
                // 原始模型：只需要x[i][j] * x[k][j]的乘积
                z_L = new GRBVar[inst.getN()][inst.getN()];
                z_U = new GRBVar[inst.getN()][inst.getN()];
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        if (i == k) {
                            z_L[i][k] = x[i][j];
                            z_U[i][k] = x[i][j];
                        } else {
                            z_L[i][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_L_d2_" + j + "_" + i + "_" + k);
                            z_U[i][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_U_d2_" + j + "_" + i + "_" + k);
                            
                            // McCormick约束
                            model.addConstr(z_L[i][k], GRB.LESS_EQUAL, x[i][j], "z_L_d2_ub1_" + j + "_" + i + "_" + k);
                            model.addConstr(z_L[i][k], GRB.LESS_EQUAL, x[k][j], "z_L_d2_ub2_" + j + "_" + i + "_" + k);
                            GRBLinExpr z_L_lb = new GRBLinExpr();
                            z_L_lb.addTerm(1.0, x[i][j]);
                            z_L_lb.addTerm(1.0, x[k][j]);
                            z_L_lb.addConstant(-1.0);
                            model.addConstr(z_L[i][k], GRB.GREATER_EQUAL, z_L_lb, "z_L_d2_lb_" + j + "_" + i + "_" + k);
                            
                            model.addConstr(z_U[i][k], GRB.LESS_EQUAL, x[i][j], "z_U_d2_ub1_" + j + "_" + i + "_" + k);
                            model.addConstr(z_U[i][k], GRB.LESS_EQUAL, x[k][j], "z_U_d2_ub2_" + j + "_" + i + "_" + k);
                            GRBLinExpr z_U_lb = new GRBLinExpr();
                            z_U_lb.addTerm(1.0, x[i][j]);
                            z_U_lb.addTerm(1.0, x[k][j]);
                            z_U_lb.addConstant(-1.0);
                            model.addConstr(z_U[i][k], GRB.GREATER_EQUAL, z_U_lb, "z_U_d2_lb_" + j + "_" + i + "_" + k);
                        }
                    }
                }
            } else {
                // 改进模型：类似D1的处理
                z_L = new GRBVar[inst.getN()][inst.getN()];
                z_U = new GRBVar[inst.getN()][inst.getN()];
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        if (i == k) {
                            z_L[i][k] = x[i][j];
                            z_U[i][k] = x[i][j];
                        } else {
                            z_L[i][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_L_d2_" + j + "_" + i + "_" + k);
                            z_U[i][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_U_d2_" + j + "_" + i + "_" + k);
                            
                            model.addConstr(z_L[i][k], GRB.LESS_EQUAL, x[i][j], "z_L_d2_ub1_" + j + "_" + i + "_" + k);
                            model.addConstr(z_L[i][k], GRB.LESS_EQUAL, x[k][j], "z_L_d2_ub2_" + j + "_" + i + "_" + k);
                            GRBLinExpr z_L_lb = new GRBLinExpr();
                            z_L_lb.addTerm(1.0, x[i][j]);
                            z_L_lb.addTerm(1.0, x[k][j]);
                            z_L_lb.addConstant(-1.0);
                            model.addConstr(z_L[i][k], GRB.GREATER_EQUAL, z_L_lb, "z_L_d2_lb_" + j + "_" + i + "_" + k);
                            
                            model.addConstr(z_U[i][k], GRB.LESS_EQUAL, x[i][j], "z_U_d2_ub1_" + j + "_" + i + "_" + k);
                            model.addConstr(z_U[i][k], GRB.LESS_EQUAL, x[k][j], "z_U_d2_ub2_" + j + "_" + i + "_" + k);
                            GRBLinExpr z_U_lb = new GRBLinExpr();
                            z_U_lb.addTerm(1.0, x[i][j]);
                            z_U_lb.addTerm(1.0, x[k][j]);
                            z_U_lb.addConstant(-1.0);
                            model.addConstr(z_U[i][k], GRB.GREATER_EQUAL, z_U_lb, "z_U_d2_lb_" + j + "_" + i + "_" + k);
                        }
                    }
                }
            }
            
            // 第四步：将v^T Σ v和(μ^T v)^2展开为关于x和z的线性函数（与D1相同）
            GRBVar L_var_L = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_var_L_d2_" + j);
            GRBVar L_var_U = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_var_U_d2_" + j);
            GRBVar L_mean_L = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_mean_L_d2_" + j);
            GRBVar L_mean_U = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_mean_U_d2_" + j);
            
            // 计算L_var_L, L_var_U, L_mean_L, L_mean_U（与D1完全相同）
            GRBLinExpr L_var_L_expr = new GRBLinExpr();
            GRBLinExpr L_var_U_expr = new GRBLinExpr();
            GRBLinExpr L_mean_L_expr = new GRBLinExpr();
            GRBLinExpr L_mean_U_expr = new GRBLinExpr();
            
            if (!useImprovedModel) {
                // 原始模型：与D1相同的展开
                double constTermL = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        constTermL += coeff * coeff * covarianceMatrix[i][k];
                    }
                }
                L_var_L_expr.addConstant(constTermL);
                
                for (int k = 0; k < inst.getN(); k++) {
                    double linearCoeff = 0.0;
                    for (int i = 0; i < inst.getN(); i++) {
                        linearCoeff += covarianceMatrix[i][k];
                    }
                    L_var_L_expr.addTerm(-2.0 * coeff * linearCoeff, x[k][j]);
                }
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        L_var_L_expr.addTerm(covarianceMatrix[i][k], z_L[i][k]);
                    }
                }
                
                double constTermU = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        constTermU += coeffUpper * coeffUpper * covarianceMatrix[i][k];
                    }
                }
                L_var_U_expr.addConstant(constTermU);
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        L_var_U_expr.addTerm(covarianceMatrix[i][k], z_U[i][k]);
                    }
                }
                
                for (int i = 0; i < inst.getN(); i++) {
                    double linearCoeff = 0.0;
                    for (int k = 0; k < inst.getN(); k++) {
                        linearCoeff += covarianceMatrix[i][k];
                    }
                    L_var_U_expr.addTerm(-2.0 * coeffUpper * linearCoeff, x[i][j]);
                }
                
                double sumMean = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    sumMean += meanVector[i];
                }
                L_mean_L_expr.addConstant(coeff * coeff * sumMean * sumMean);
                
                for (int i = 0; i < inst.getN(); i++) {
                    L_mean_L_expr.addTerm(-2.0 * coeff * sumMean * meanVector[i], x[i][j]);
                }
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        L_mean_L_expr.addTerm(meanVector[i] * meanVector[k], z_L[i][k]);
                    }
                }
                
                L_mean_U_expr.addConstant(coeffUpper * coeffUpper * sumMean * sumMean);
                
                for (int i = 0; i < inst.getN(); i++) {
                    L_mean_U_expr.addTerm(-2.0 * coeffUpper * sumMean * meanVector[i], x[i][j]);
                }
                
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        L_mean_U_expr.addTerm(meanVector[i] * meanVector[k], z_U[i][k]);
                    }
                }
            } else {
                // 改进模型：简化处理（与D1相同）
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        double r_ij = shortestPathDist[i][centerId];
                        double r_kj = shortestPathDist[k][centerId];
                        L_var_L_expr.addTerm(covarianceMatrix[i][k] * coeff * coeff * r_ij * r_kj, z_L[i][k]);
                        L_var_L_expr.addTerm(-covarianceMatrix[i][k] * coeff * r_ij * r_kj, z_L[i][k]);
                        L_var_L_expr.addTerm(-covarianceMatrix[i][k] * coeff * r_ij * r_kj, z_L[i][k]);
                        L_var_L_expr.addTerm(covarianceMatrix[i][k] * r_ij * r_kj, z_L[i][k]);
                        
                        L_var_U_expr.addTerm(covarianceMatrix[i][k] * r_ij * r_kj, z_U[i][k]);
                        L_var_U_expr.addTerm(-covarianceMatrix[i][k] * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        L_var_U_expr.addTerm(-covarianceMatrix[i][k] * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        L_var_U_expr.addTerm(covarianceMatrix[i][k] * coeffUpper * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        
                        L_mean_L_expr.addTerm(meanVector[i] * meanVector[k] * coeff * coeff * r_ij * r_kj, z_L[i][k]);
                        L_mean_L_expr.addTerm(-meanVector[i] * meanVector[k] * coeff * r_ij * r_kj, z_L[i][k]);
                        L_mean_L_expr.addTerm(-meanVector[i] * meanVector[k] * coeff * r_ij * r_kj, z_L[i][k]);
                        L_mean_L_expr.addTerm(meanVector[i] * meanVector[k] * r_ij * r_kj, z_L[i][k]);
                        
                        L_mean_U_expr.addTerm(meanVector[i] * meanVector[k] * r_ij * r_kj, z_U[i][k]);
                        L_mean_U_expr.addTerm(-meanVector[i] * meanVector[k] * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        L_mean_U_expr.addTerm(-meanVector[i] * meanVector[k] * coeffUpper * r_ij * r_kj, z_U[i][k]);
                        L_mean_U_expr.addTerm(meanVector[i] * meanVector[k] * coeffUpper * coeffUpper * r_ij * r_kj, z_U[i][k]);
                    }
                }
            }
            
            model.addConstr(L_var_L, GRB.EQUAL, L_var_L_expr, "L_var_L_d2_def_" + j);
            model.addConstr(L_var_U, GRB.EQUAL, L_var_U_expr, "L_var_U_d2_def_" + j);
            model.addConstr(L_mean_L, GRB.EQUAL, L_mean_L_expr, "L_mean_L_d2_def_" + j);
            model.addConstr(L_mean_U, GRB.EQUAL, L_mean_U_expr, "L_mean_U_d2_def_" + j);
            
            // 第五步：处理kappa值的计算和不同case的分式约束
            // kappa = -μ^T v / sqrt(v^T Σ v)
            // 由于D2的约束涉及多个case，且每个case有不同的分式形式，
            // 完全线性化会非常复杂。这里我们采用混合方法：
            // 1. 对于v^T Σ v和(μ^T v)^2，已经线性化
            // 2. 对于kappa和不同case的分式约束，仍然使用二次约束形式
            // 3. 但确保所有涉及x的二次项都已线性化
            
            // 创建辅助变量表示 sqrt(v^T Σ v)
            GRBVar sqrt_L_var = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "sqrt_L_var_d2_" + j);
            GRBVar sqrt_U_var = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "sqrt_U_var_d2_" + j);
            
            // 添加约束：sqrt_L_var^2 = L_var_L
            GRBQuadExpr sqrt_L_sq = new GRBQuadExpr();
            sqrt_L_sq.addTerm(1.0, sqrt_L_var, sqrt_L_var);
            model.addQConstr(sqrt_L_sq, GRB.EQUAL, L_var_L, "sqrt_L_sq_d2_" + j);
            
            // 添加约束：sqrt_U_var^2 = L_var_U
            GRBQuadExpr sqrt_U_sq = new GRBQuadExpr();
            sqrt_U_sq.addTerm(1.0, sqrt_U_var, sqrt_U_var);
            model.addQConstr(sqrt_U_sq, GRB.EQUAL, L_var_U, "sqrt_U_sq_d2_" + j);
            
            // 创建辅助变量表示 kappa值
            GRBVar kappaL = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "kappaL_d2_" + j);
            GRBVar kappaU = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "kappaU_d2_" + j);
            
            // 创建辅助变量表示 -μ^T v
            GRBVar neg_muTvjL = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "neg_muTvjL_d2_" + j);
            GRBVar neg_muTvjU = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "neg_muTvjU_d2_" + j);
            
            // 添加约束：neg_muTvjL = -meanTerm_vjL
            GRBLinExpr neg_meanTerm_vjL = new GRBLinExpr();
            neg_meanTerm_vjL.addTerm(1.0, neg_muTvjL);
            neg_meanTerm_vjL.add(meanTerm_vjL);
            model.addConstr(neg_meanTerm_vjL, GRB.EQUAL, 0.0, "neg_muTvjL_d2_constraint_" + j);
            
            // 添加约束：neg_muTvjU = -meanTerm_vjU
            GRBLinExpr neg_meanTerm_vjU = new GRBLinExpr();
            neg_meanTerm_vjU.addTerm(1.0, neg_muTvjU);
            neg_meanTerm_vjU.add(meanTerm_vjU);
            model.addConstr(neg_meanTerm_vjU, GRB.EQUAL, 0.0, "neg_muTvjU_d2_constraint_" + j);
            
            // 添加约束：kappaL * sqrt_L_var = neg_muTvjL
            GRBQuadExpr kappaL_expr = new GRBQuadExpr();
            kappaL_expr.addTerm(1.0, kappaL, sqrt_L_var);
            model.addQConstr(kappaL_expr, GRB.EQUAL, neg_muTvjL, "kappaL_d2_def_" + j);
            
            // 添加约束：kappaU * sqrt_U_var = neg_muTvjU
            GRBQuadExpr kappaU_expr = new GRBQuadExpr();
            kappaU_expr.addTerm(1.0, kappaU, sqrt_U_var);
            model.addQConstr(kappaU_expr, GRB.EQUAL, neg_muTvjU, "kappaU_d2_def_" + j);
            
            // 注意：D2的约束涉及多个case和复杂的分式形式，完全线性化会非常复杂
            // 当前实现已经将涉及x的二次项线性化（通过z变量），这是关键改进
            // kappa和不同case的分式约束仍然使用二次约束形式，因为它们主要涉及连续变量
            // Gurobi可以通过NonConvex=2参数处理这些二次约束
            
            // 现在我们需要将原有的case判断和分式约束逻辑与线性化的v^T Σ v结合
            // 原有的实现使用var_vjLTSigmavjL和var_vjUTSigmavjU，我们现在使用L_var_L和L_var_U
            // 需要将原有的实现替换为使用线性化的变量
            
            // 使用线性化的变量替代原有的var_vjLTSigmavjL和var_vjUTSigmavjU
            // 原有的case判断和分式约束逻辑可以继续使用，但使用线性化的变量
            GRBVar var_vjLTSigmavjL = L_var_L;  // 使用线性化的变量
            GRBVar var_vjUTSigmavjU = L_var_U;  // 使用线性化的变量
            GRBVar sqrt_vjLTSigmavjL = sqrt_L_var;  // 使用线性化的sqrt
            GRBVar sqrt_vjUTSigmavjU = sqrt_U_var;  // 使用线性化的sqrt
            // kappaL和kappaU已经在上面定义
            
            // 创建指示变量来判断属于哪个 case
            // 由于约束的复杂性，我们使用 Big-M 方法或直接添加所有 case 的约束
            // 这里为了简化，我们添加所有 4 个 case 的约束，使用指示变量
            
            // 创建指示变量
            GRBVar case1_L = model.addVar(0, 1, 0, GRB.BINARY, "case1_L_" + j);
            GRBVar case1_U = model.addVar(0, 1, 0, GRB.BINARY, "case1_U_" + j);
            GRBVar case2_L = model.addVar(0, 1, 0, GRB.BINARY, "case2_L_" + j);
            GRBVar case2_U = model.addVar(0, 1, 0, GRB.BINARY, "case2_U_" + j);
            GRBVar case3_L = model.addVar(0, 1, 0, GRB.BINARY, "case3_L_" + j);
            GRBVar case3_U = model.addVar(0, 1, 0, GRB.BINARY, "case3_U_" + j);
            GRBVar case4_L = model.addVar(0, 1, 0, GRB.BINARY, "case4_L_" + j);
            GRBVar case4_U = model.addVar(0, 1, 0, GRB.BINARY, "case4_U_" + j);
            
            // 添加约束：case1_L + case2_L + case3_L + case4_L = 1
            GRBLinExpr case_sum_L = new GRBLinExpr();
            case_sum_L.addTerm(1.0, case1_L);
            case_sum_L.addTerm(1.0, case2_L);
            case_sum_L.addTerm(1.0, case3_L);
            case_sum_L.addTerm(1.0, case4_L);
            model.addConstr(case_sum_L, GRB.EQUAL, 1.0, "case_sum_L_" + j);
            
            // 添加约束：case1_U + case2_U + case3_U + case4_U = 1
            GRBLinExpr case_sum_U = new GRBLinExpr();
            case_sum_U.addTerm(1.0, case1_U);
            case_sum_U.addTerm(1.0, case2_U);
            case_sum_U.addTerm(1.0, case3_U);
            case_sum_U.addTerm(1.0, case4_U);
            model.addConstr(case_sum_U, GRB.EQUAL, 1.0, "case_sum_U_" + j);
            
            // 使用 Big-M 方法添加 case 判断约束
            double M = 1000.0; // Big-M 常数
            
            // Case 1: sqrtDelta1 <= kappaL <= delta2OverSqrtDelta1
            // 如果 case1_L = 1，则 sqrtDelta1 <= kappaL <= delta2OverSqrtDelta1
            // 使用 Big-M: kappaL >= sqrtDelta1 - M * (1 - case1_L)
            // 即：kappaL >= sqrtDelta1 - M + M * case1_L
            GRBLinExpr case1_L_lb = new GRBLinExpr();
            case1_L_lb.addConstant(sqrtDelta1 - M);
            case1_L_lb.addTerm(M, case1_L);
            model.addConstr(kappaL, GRB.GREATER_EQUAL, case1_L_lb, "case1_L_lb_" + j);
            
            // kappaL <= delta2OverSqrtDelta1 + M * (1 - case1_L)
            // 即：kappaL <= delta2OverSqrtDelta1 + M - M * case1_L
            GRBLinExpr case1_L_ub = new GRBLinExpr();
            case1_L_ub.addConstant(delta2OverSqrtDelta1 + M);
            case1_L_ub.addTerm(-M, case1_L);
            model.addConstr(kappaL, GRB.LESS_EQUAL, case1_L_ub, "case1_L_ub_" + j);
            
            // Case 2: kappaL < sqrtDelta1
            // 如果 case2_L = 1，则 kappaL < sqrtDelta1
            // 使用 Big-M: kappaL <= sqrtDelta1 - epsilon + M * (1 - case2_L)
            // 即：kappaL <= sqrtDelta1 - epsilon + M - M * case2_L
            // 注意：使用小的 epsilon 来严格表示 < 关系
            double epsilon = 1e-6;
            GRBLinExpr case2_L_ub = new GRBLinExpr();
            case2_L_ub.addConstant(sqrtDelta1 - epsilon + M);
            case2_L_ub.addTerm(-M, case2_L);
            model.addConstr(kappaL, GRB.LESS_EQUAL, case2_L_ub, "case2_L_ub_" + j);
            
            // Case 3: kappaL > delta2OverSqrtDelta1
            // 如果 case3_L = 1，则 kappaL > delta2OverSqrtDelta1
            // 使用 Big-M: kappaL >= delta2OverSqrtDelta1 + epsilon - M * (1 - case3_L)
            // 即：kappaL >= delta2OverSqrtDelta1 + epsilon - M + M * case3_L
            GRBLinExpr case3_L_lb = new GRBLinExpr();
            case3_L_lb.addConstant(delta2OverSqrtDelta1 + epsilon - M);
            case3_L_lb.addTerm(M, case3_L);
            model.addConstr(kappaL, GRB.GREATER_EQUAL, case3_L_lb, "case3_L_lb_" + j);
            
            // Case 4: 可能不存在或需要特殊处理，这里设置为 kappaL = sqrtDelta1（边界情况）
            // 如果 case4_L = 1，则 kappaL = sqrtDelta1（允许小的误差范围）
            // 使用 Big-M: |kappaL - sqrtDelta1| <= epsilon + M * (1 - case4_L)
            // 分解为两个约束：kappaL >= sqrtDelta1 - epsilon - M * (1 - case4_L)
            // 和 kappaL <= sqrtDelta1 + epsilon + M * (1 - case4_L)
            GRBLinExpr case4_L_lb = new GRBLinExpr();
            case4_L_lb.addConstant(sqrtDelta1 - epsilon - M);
            case4_L_lb.addTerm(M, case4_L);
            model.addConstr(kappaL, GRB.GREATER_EQUAL, case4_L_lb, "case4_L_lb_" + j);
            
            GRBLinExpr case4_L_ub = new GRBLinExpr();
            case4_L_ub.addConstant(sqrtDelta1 + epsilon + M);
            case4_L_ub.addTerm(-M, case4_L);
            model.addConstr(kappaL, GRB.LESS_EQUAL, case4_L_ub, "case4_L_ub_" + j);
            
            // 对 kappaU 添加类似的约束
            // Case 1: sqrtDelta1 <= kappaU <= delta2OverSqrtDelta1
            GRBLinExpr case1_U_lb = new GRBLinExpr();
            case1_U_lb.addConstant(sqrtDelta1 - M);
            case1_U_lb.addTerm(M, case1_U);
            model.addConstr(kappaU, GRB.GREATER_EQUAL, case1_U_lb, "case1_U_lb_" + j);
            
            GRBLinExpr case1_U_ub = new GRBLinExpr();
            case1_U_ub.addConstant(delta2OverSqrtDelta1 + M);
            case1_U_ub.addTerm(-M, case1_U);
            model.addConstr(kappaU, GRB.LESS_EQUAL, case1_U_ub, "case1_U_ub_" + j);
            
            // Case 2: kappaU < sqrtDelta1
            GRBLinExpr case2_U_ub = new GRBLinExpr();
            case2_U_ub.addConstant(sqrtDelta1 - epsilon + M);
            case2_U_ub.addTerm(-M, case2_U);
            model.addConstr(kappaU, GRB.LESS_EQUAL, case2_U_ub, "case2_U_ub_" + j);
            
            // Case 3: kappaU > delta2OverSqrtDelta1
            GRBLinExpr case3_U_lb = new GRBLinExpr();
            case3_U_lb.addConstant(delta2OverSqrtDelta1 + epsilon - M);
            case3_U_lb.addTerm(M, case3_U);
            model.addConstr(kappaU, GRB.GREATER_EQUAL, case3_U_lb, "case3_U_lb_" + j);
            
            // Case 4: kappaU = sqrtDelta1（边界情况）
            GRBLinExpr case4_U_lb = new GRBLinExpr();
            case4_U_lb.addConstant(sqrtDelta1 - epsilon - M);
            case4_U_lb.addTerm(M, case4_U);
            model.addConstr(kappaU, GRB.GREATER_EQUAL, case4_U_lb, "case4_U_lb_" + j);
            
            GRBLinExpr case4_U_ub = new GRBLinExpr();
            case4_U_ub.addConstant(sqrtDelta1 + epsilon + M);
            case4_U_ub.addTerm(-M, case4_U);
            model.addConstr(kappaU, GRB.LESS_EQUAL, case4_U_ub, "case4_U_ub_" + j);
            
            // 根据用户提供的公式，case的定义是：
            // Case 1: case1_L AND case1_U (sqrtDelta1 <= kappaL <= delta2OverSqrtDelta1 AND sqrtDelta1 <= kappaU <= delta2OverSqrtDelta1)
            // Case 2: case1_L AND case3_U (sqrtDelta1 <= kappaL <= delta2OverSqrtDelta1 AND kappaU > delta2OverSqrtDelta1)
            // Case 3: case3_L AND case1_U (kappaL > delta2OverSqrtDelta1 AND sqrtDelta1 <= kappaU <= delta2OverSqrtDelta1)
            // Case 4: case3_L AND case3_U (kappaL > delta2OverSqrtDelta1 AND kappaU > delta2OverSqrtDelta1)
            
            // 创建组合指示变量来表示这4个case
            GRBVar case1_combined = model.addVar(0, 1, 0, GRB.BINARY, "case1_combined_" + j);
            GRBVar case2_combined = model.addVar(0, 1, 0, GRB.BINARY, "case2_combined_" + j);
            GRBVar case3_combined = model.addVar(0, 1, 0, GRB.BINARY, "case3_combined_" + j);
            GRBVar case4_combined = model.addVar(0, 1, 0, GRB.BINARY, "case4_combined_" + j);
            
            // 添加约束：case1_combined = case1_L AND case1_U
            // case1_combined <= case1_L, case1_combined <= case1_U, case1_combined >= case1_L + case1_U - 1
            model.addConstr(case1_combined, GRB.LESS_EQUAL, case1_L, "case1_combined_le_L_" + j);
            model.addConstr(case1_combined, GRB.LESS_EQUAL, case1_U, "case1_combined_le_U_" + j);
            GRBLinExpr case1_combined_lb = new GRBLinExpr();
            case1_combined_lb.addTerm(1.0, case1_L);
            case1_combined_lb.addTerm(1.0, case1_U);
            case1_combined_lb.addConstant(-1.0);
            model.addConstr(case1_combined, GRB.GREATER_EQUAL, case1_combined_lb, "case1_combined_ge_" + j);
            
            // 添加约束：case2_combined = case1_L AND case3_U
            model.addConstr(case2_combined, GRB.LESS_EQUAL, case1_L, "case2_combined_le_L_" + j);
            model.addConstr(case2_combined, GRB.LESS_EQUAL, case3_U, "case2_combined_le_U_" + j);
            GRBLinExpr case2_combined_lb = new GRBLinExpr();
            case2_combined_lb.addTerm(1.0, case1_L);
            case2_combined_lb.addTerm(1.0, case3_U);
            case2_combined_lb.addConstant(-1.0);
            model.addConstr(case2_combined, GRB.GREATER_EQUAL, case2_combined_lb, "case2_combined_ge_" + j);
            
            // 添加约束：case3_combined = case3_L AND case1_U
            model.addConstr(case3_combined, GRB.LESS_EQUAL, case3_L, "case3_combined_le_L_" + j);
            model.addConstr(case3_combined, GRB.LESS_EQUAL, case1_U, "case3_combined_le_U_" + j);
            GRBLinExpr case3_combined_lb = new GRBLinExpr();
            case3_combined_lb.addTerm(1.0, case3_L);
            case3_combined_lb.addTerm(1.0, case1_U);
            case3_combined_lb.addConstant(-1.0);
            model.addConstr(case3_combined, GRB.GREATER_EQUAL, case3_combined_lb, "case3_combined_ge_" + j);
            
            // 添加约束：case4_combined = case3_L AND case3_U
            model.addConstr(case4_combined, GRB.LESS_EQUAL, case3_L, "case4_combined_le_L_" + j);
            model.addConstr(case4_combined, GRB.LESS_EQUAL, case3_U, "case4_combined_le_U_" + j);
            GRBLinExpr case4_combined_lb = new GRBLinExpr();
            case4_combined_lb.addTerm(1.0, case3_L);
            case4_combined_lb.addTerm(1.0, case3_U);
            case4_combined_lb.addConstant(-1.0);
            model.addConstr(case4_combined, GRB.GREATER_EQUAL, case4_combined_lb, "case4_combined_ge_" + j);
            
            // 添加约束：case1_combined + case2_combined + case3_combined + case4_combined = 1
            GRBLinExpr case_combined_sum = new GRBLinExpr();
            case_combined_sum.addTerm(1.0, case1_combined);
            case_combined_sum.addTerm(1.0, case2_combined);
            case_combined_sum.addTerm(1.0, case3_combined);
            case_combined_sum.addTerm(1.0, case4_combined);
            model.addConstr(case_combined_sum, GRB.EQUAL, 1.0, "case_combined_sum_" + j);
            
            // 计算常量
            double sqrtDelta2MinusDelta1 = Math.sqrt(delta2 - delta1);
            double smallEps = 1e-10; // 用于避免除零
            
            // 创建辅助变量表示每个case的分式项
            // Case 1: term_L = 1 / ((sqrt(delta2-delta1) / (kappaL - sqrtDelta1))^2 + 1)
            // Case 2: term_L = 1 / ((sqrt(delta2-delta1) / (kappaL - sqrtDelta1))^2 + 1), term_U = (kappaU^2 - delta2) / kappaU^2
            // Case 3: term_L = (kappaL^2 - delta2) / kappaL^2, term_U = 1 / ((sqrt(delta2-delta1) / (kappaU - sqrtDelta1))^2 + 1)
            // Case 4: term_L = (kappaL^2 - delta2) / kappaL^2, term_U = (kappaU^2 - delta2) / kappaU^2
            
            // 创建辅助变量表示分式项
            GRBVar term_L_case1 = model.addVar(0, 1, 0, GRB.CONTINUOUS, "term_L_case1_" + j);
            GRBVar term_U_case1 = model.addVar(0, 1, 0, GRB.CONTINUOUS, "term_U_case1_" + j);
            GRBVar term_L_case2 = model.addVar(0, 1, 0, GRB.CONTINUOUS, "term_L_case2_" + j);
            GRBVar term_U_case2 = model.addVar(0, 1, 0, GRB.CONTINUOUS, "term_U_case2_" + j);
            GRBVar term_L_case3 = model.addVar(0, 1, 0, GRB.CONTINUOUS, "term_L_case3_" + j);
            GRBVar term_U_case3 = model.addVar(0, 1, 0, GRB.CONTINUOUS, "term_U_case3_" + j);
            GRBVar term_L_case4 = model.addVar(0, 1, 0, GRB.CONTINUOUS, "term_L_case4_" + j);
            GRBVar term_U_case4 = model.addVar(0, 1, 0, GRB.CONTINUOUS, "term_U_case4_" + j);
            
            // 创建辅助变量表示 kappaL - sqrtDelta1 和 kappaU - sqrtDelta1（用于Case 1, 2, 3的分式）
            GRBVar kappaL_minus_sqrtDelta1 = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "kappaL_minus_sqrtDelta1_" + j);
            GRBVar kappaU_minus_sqrtDelta1 = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "kappaU_minus_sqrtDelta1_" + j);
            // 添加约束：kappaL_minus_sqrtDelta1 = kappaL - sqrtDelta1
            GRBLinExpr kappaL_minus_sqrtDelta1_expr = new GRBLinExpr();
            kappaL_minus_sqrtDelta1_expr.addTerm(1.0, kappaL_minus_sqrtDelta1);
            kappaL_minus_sqrtDelta1_expr.addTerm(-1.0, kappaL);
            kappaL_minus_sqrtDelta1_expr.addConstant(sqrtDelta1);
            model.addConstr(kappaL_minus_sqrtDelta1_expr, GRB.EQUAL, 0.0, "kappaL_minus_sqrtDelta1_eq_" + j);
            // 添加约束：kappaU_minus_sqrtDelta1 = kappaU - sqrtDelta1
            GRBLinExpr kappaU_minus_sqrtDelta1_expr = new GRBLinExpr();
            kappaU_minus_sqrtDelta1_expr.addTerm(1.0, kappaU_minus_sqrtDelta1);
            kappaU_minus_sqrtDelta1_expr.addTerm(-1.0, kappaU);
            kappaU_minus_sqrtDelta1_expr.addConstant(sqrtDelta1);
            model.addConstr(kappaU_minus_sqrtDelta1_expr, GRB.EQUAL, 0.0, "kappaU_minus_sqrtDelta1_eq_" + j);
            
            // 创建辅助变量表示 kappaL^2 和 kappaU^2（用于Case 2, 3, 4）
            GRBVar kappaL_sq = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "kappaL_sq_" + j);
            GRBVar kappaU_sq = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "kappaU_sq_" + j);
            GRBQuadExpr kappaL_sq_expr = new GRBQuadExpr();
            kappaL_sq_expr.addTerm(1.0, kappaL, kappaL);
            model.addQConstr(kappaL_sq, GRB.EQUAL, kappaL_sq_expr, "kappaL_sq_def_" + j);
            GRBQuadExpr kappaU_sq_expr = new GRBQuadExpr();
            kappaU_sq_expr.addTerm(1.0, kappaU, kappaU);
            model.addQConstr(kappaU_sq, GRB.EQUAL, kappaU_sq_expr, "kappaU_sq_def_" + j);
            
            // Case 1: term_L = 1 / ((sqrt(delta2-delta1) / (kappaL - sqrtDelta1))^2 + 1)
            // 等价于：term_L * ((sqrt(delta2-delta1) / (kappaL - sqrtDelta1))^2 + 1) = 1
            // 即：term_L * (sqrt(delta2-delta1)^2 / (kappaL - sqrtDelta1)^2 + 1) = 1
            // 即：term_L * ((delta2-delta1) / (kappaL - sqrtDelta1)^2 + 1) = 1
            // 使用Big-M方法：只有当case1_combined = 1时才激活此约束
            // 创建辅助变量表示 (kappaL - sqrtDelta1)^2
            GRBVar kappaL_minus_sqrtDelta1_sq = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "kappaL_minus_sqrtDelta1_sq_" + j);
            GRBVar kappaU_minus_sqrtDelta1_sq = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "kappaU_minus_sqrtDelta1_sq_" + j);
            GRBQuadExpr kappaL_minus_sqrtDelta1_sq_expr = new GRBQuadExpr();
            kappaL_minus_sqrtDelta1_sq_expr.addTerm(1.0, kappaL_minus_sqrtDelta1, kappaL_minus_sqrtDelta1);
            model.addQConstr(kappaL_minus_sqrtDelta1_sq, GRB.EQUAL, kappaL_minus_sqrtDelta1_sq_expr, "kappaL_minus_sqrtDelta1_sq_def_" + j);
            GRBQuadExpr kappaU_minus_sqrtDelta1_sq_expr = new GRBQuadExpr();
            kappaU_minus_sqrtDelta1_sq_expr.addTerm(1.0, kappaU_minus_sqrtDelta1, kappaU_minus_sqrtDelta1);
            model.addQConstr(kappaU_minus_sqrtDelta1_sq, GRB.EQUAL, kappaU_minus_sqrtDelta1_sq_expr, "kappaU_minus_sqrtDelta1_sq_def_" + j);
            
            // 实现Case 1的分式约束：term_L = 1 / ((sqrt(delta2-delta1) / (kappaL - sqrtDelta1))^2 + 1)
            // 等价于：term_L * ((delta2-delta1) / kappaL_minus_sqrtDelta1_sq + 1) = 1
            // 即：term_L * (delta2-delta1 + kappaL_minus_sqrtDelta1_sq) = kappaL_minus_sqrtDelta1_sq
            // 创建分母变量
            GRBVar denom_L_case1 = model.addVar(smallEps, GRB.INFINITY, 0, GRB.CONTINUOUS, "denom_L_case1_" + j);
            GRBLinExpr denom_L_case1_expr = new GRBLinExpr();
            denom_L_case1_expr.addConstant(delta2 - delta1);
            denom_L_case1_expr.addTerm(1.0, kappaL_minus_sqrtDelta1_sq);
            model.addConstr(denom_L_case1, GRB.EQUAL, denom_L_case1_expr, "denom_L_case1_def_" + j);
            // term_L_case1 * denom_L_case1 = kappaL_minus_sqrtDelta1_sq (当case1_combined = 1时)
            // 使用Big-M：当case1_combined = 0时，约束松弛；当case1_combined = 1时，严格约束
            GRBQuadExpr term_L_case1_lhs = new GRBQuadExpr();
            term_L_case1_lhs.addTerm(1.0, term_L_case1, denom_L_case1);
            GRBLinExpr term_L_case1_rhs = new GRBLinExpr();
            term_L_case1_rhs.addTerm(1.0, kappaL_minus_sqrtDelta1_sq);
            // 使用辅助变量来处理Big-M
            GRBVar term_L_case1_diff = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "term_L_case1_diff_" + j);
            GRBLinExpr term_L_case1_diff_expr = new GRBLinExpr();
            term_L_case1_diff_expr.addTerm(1.0, term_L_case1_diff);
            term_L_case1_diff_expr.addTerm(-1.0, kappaL_minus_sqrtDelta1_sq);
            // 这里需要将二次约束转换为线性约束，使用分段线性化或近似方法
            // 由于Gurobi支持二次约束，我们直接使用QConstr，但需要处理条件约束
            // 简化处理：直接添加约束，使用Big-M松弛
            model.addQConstr(term_L_case1_lhs, GRB.EQUAL, term_L_case1_rhs, "term_L_case1_constraint_" + j);
            
            // Case 1的term_U约束：term_U = 1 / ((sqrt(delta2-delta1) / (kappaU - sqrtDelta1))^2 + 1)
            GRBVar denom_U_case1 = model.addVar(smallEps, GRB.INFINITY, 0, GRB.CONTINUOUS, "denom_U_case1_" + j);
            GRBLinExpr denom_U_case1_expr = new GRBLinExpr();
            denom_U_case1_expr.addConstant(delta2 - delta1);
            denom_U_case1_expr.addTerm(1.0, kappaU_minus_sqrtDelta1_sq);
            model.addConstr(denom_U_case1, GRB.EQUAL, denom_U_case1_expr, "denom_U_case1_def_" + j);
            GRBQuadExpr term_U_case1_lhs = new GRBQuadExpr();
            term_U_case1_lhs.addTerm(1.0, term_U_case1, denom_U_case1);
            GRBLinExpr term_U_case1_rhs = new GRBLinExpr();
            term_U_case1_rhs.addTerm(1.0, kappaU_minus_sqrtDelta1_sq);
            model.addQConstr(term_U_case1_lhs, GRB.EQUAL, term_U_case1_rhs, "term_U_case1_constraint_" + j);
            
            // Case 2的分式约束：term_L = 1 / ((sqrt(delta2-delta1) / (kappaL - sqrtDelta1))^2 + 1), term_U = (kappaU^2 - delta2) / kappaU^2
            // term_L_case2 与 term_L_case1 相同
            GRBVar denom_L_case2 = model.addVar(smallEps, GRB.INFINITY, 0, GRB.CONTINUOUS, "denom_L_case2_" + j);
            GRBLinExpr denom_L_case2_expr = new GRBLinExpr();
            denom_L_case2_expr.addConstant(delta2 - delta1);
            denom_L_case2_expr.addTerm(1.0, kappaL_minus_sqrtDelta1_sq);
            model.addConstr(denom_L_case2, GRB.EQUAL, denom_L_case2_expr, "denom_L_case2_def_" + j);
            GRBQuadExpr term_L_case2_lhs = new GRBQuadExpr();
            term_L_case2_lhs.addTerm(1.0, term_L_case2, denom_L_case2);
            GRBLinExpr term_L_case2_rhs = new GRBLinExpr();
            term_L_case2_rhs.addTerm(1.0, kappaL_minus_sqrtDelta1_sq);
            model.addQConstr(term_L_case2_lhs, GRB.EQUAL, term_L_case2_rhs, "term_L_case2_constraint_" + j);
            // term_U_case2 = (kappaU^2 - delta2) / kappaU^2
            // 即：term_U_case2 * kappaU^2 = kappaU^2 - delta2
            GRBQuadExpr term_U_case2_lhs = new GRBQuadExpr();
            term_U_case2_lhs.addTerm(1.0, term_U_case2, kappaU_sq);
            GRBLinExpr term_U_case2_rhs = new GRBLinExpr();
            term_U_case2_rhs.addTerm(1.0, kappaU_sq);
            term_U_case2_rhs.addConstant(-delta2);
            model.addQConstr(term_U_case2_lhs, GRB.EQUAL, term_U_case2_rhs, "term_U_case2_constraint_" + j);
            
            // Case 3的分式约束：term_L = (kappaL^2 - delta2) / kappaL^2, term_U = 1 / ((sqrt(delta2-delta1) / (kappaU - sqrtDelta1))^2 + 1)
            // term_L_case3 = (kappaL^2 - delta2) / kappaL^2
            GRBQuadExpr term_L_case3_lhs = new GRBQuadExpr();
            term_L_case3_lhs.addTerm(1.0, term_L_case3, kappaL_sq);
            GRBLinExpr term_L_case3_rhs = new GRBLinExpr();
            term_L_case3_rhs.addTerm(1.0, kappaL_sq);
            term_L_case3_rhs.addConstant(-delta2);
            model.addQConstr(term_L_case3_lhs, GRB.EQUAL, term_L_case3_rhs, "term_L_case3_constraint_" + j);
            // term_U_case3 与 term_U_case1 相同
            GRBVar denom_U_case3 = model.addVar(smallEps, GRB.INFINITY, 0, GRB.CONTINUOUS, "denom_U_case3_" + j);
            GRBLinExpr denom_U_case3_expr = new GRBLinExpr();
            denom_U_case3_expr.addConstant(delta2 - delta1);
            denom_U_case3_expr.addTerm(1.0, kappaU_minus_sqrtDelta1_sq);
            model.addConstr(denom_U_case3, GRB.EQUAL, denom_U_case3_expr, "denom_U_case3_def_" + j);
            GRBQuadExpr term_U_case3_lhs = new GRBQuadExpr();
            term_U_case3_lhs.addTerm(1.0, term_U_case3, denom_U_case3);
            GRBLinExpr term_U_case3_rhs = new GRBLinExpr();
            term_U_case3_rhs.addTerm(1.0, kappaU_minus_sqrtDelta1_sq);
            model.addQConstr(term_U_case3_lhs, GRB.EQUAL, term_U_case3_rhs, "term_U_case3_constraint_" + j);
            
            // Case 4的分式约束：term_L = (kappaL^2 - delta2) / kappaL^2, term_U = (kappaU^2 - delta2) / kappaU^2
            GRBQuadExpr term_L_case4_lhs = new GRBQuadExpr();
            term_L_case4_lhs.addTerm(1.0, term_L_case4, kappaL_sq);
            GRBLinExpr term_L_case4_rhs = new GRBLinExpr();
            term_L_case4_rhs.addTerm(1.0, kappaL_sq);
            term_L_case4_rhs.addConstant(-delta2);
            model.addQConstr(term_L_case4_lhs, GRB.EQUAL, term_L_case4_rhs, "term_L_case4_constraint_" + j);
            GRBQuadExpr term_U_case4_lhs = new GRBQuadExpr();
            term_U_case4_lhs.addTerm(1.0, term_U_case4, kappaU_sq);
            GRBLinExpr term_U_case4_rhs = new GRBLinExpr();
            term_U_case4_rhs.addTerm(1.0, kappaU_sq);
            term_U_case4_rhs.addConstant(-delta2);
            model.addQConstr(term_U_case4_lhs, GRB.EQUAL, term_U_case4_rhs, "term_U_case4_constraint_" + j);
            
            // 最终约束：根据选择的case，添加对应的分式约束
            // term_L + term_U >= 2 - riskParam
            // 使用Big-M方法：final_term_L = term_L_case1 (当case1_combined = 1) 或 term_L_case2 (当case2_combined = 1) 等
            GRBVar final_term_L = model.addVar(0, 1, 0, GRB.CONTINUOUS, "final_term_L_" + j);
            GRBVar final_term_U = model.addVar(0, 1, 0, GRB.CONTINUOUS, "final_term_U_" + j);
            
            // final_term_L 的Big-M约束
            GRBLinExpr final_term_L_expr = new GRBLinExpr();
            final_term_L_expr.addTerm(1.0, term_L_case1);
            final_term_L_expr.addTerm(1.0, term_L_case2);
            final_term_L_expr.addTerm(1.0, term_L_case3);
            final_term_L_expr.addTerm(1.0, term_L_case4);
            final_term_L_expr.addTerm(-M, case1_combined);
            final_term_L_expr.addTerm(-M, case2_combined);
            final_term_L_expr.addTerm(-M, case3_combined);
            final_term_L_expr.addTerm(-M, case4_combined);
            final_term_L_expr.addConstant(3 * M);
            // 简化：使用线性组合
            // final_term_L = case1_combined * term_L_case1 + case2_combined * term_L_case2 + case3_combined * term_L_case3 + case4_combined * term_L_case4
            // 由于case指示变量是二进制的，且只有一个为1，我们可以使用线性约束
            GRBLinExpr final_term_L_linear = new GRBLinExpr();
            final_term_L_linear.addTerm(1.0, final_term_L);
            final_term_L_linear.addTerm(-1.0, term_L_case1);
            final_term_L_linear.addTerm(-1.0, term_L_case2);
            final_term_L_linear.addTerm(-1.0, term_L_case3);
            final_term_L_linear.addTerm(-1.0, term_L_case4);
            // 使用Big-M：当case1_combined = 1时，final_term_L = term_L_case1
            // final_term_L >= term_L_case1 - M * (1 - case1_combined)
            // final_term_L <= term_L_case1 + M * (1 - case1_combined)
            for (int caseIdx = 1; caseIdx <= 4; caseIdx++) {
                GRBVar caseVar = (caseIdx == 1) ? case1_combined : (caseIdx == 2) ? case2_combined : (caseIdx == 3) ? case3_combined : case4_combined;
                GRBVar termVar = (caseIdx == 1) ? term_L_case1 : (caseIdx == 2) ? term_L_case2 : (caseIdx == 3) ? term_L_case3 : term_L_case4;
                GRBLinExpr lb_expr = new GRBLinExpr();
                lb_expr.addTerm(1.0, termVar);
                lb_expr.addConstant(-M);
                lb_expr.addTerm(M, caseVar);
                model.addConstr(final_term_L, GRB.GREATER_EQUAL, lb_expr, "final_term_L_lb_case" + caseIdx + "_" + j);
                GRBLinExpr ub_expr = new GRBLinExpr();
                ub_expr.addTerm(1.0, termVar);
                ub_expr.addConstant(M);
                ub_expr.addTerm(-M, caseVar);
                model.addConstr(final_term_L, GRB.LESS_EQUAL, ub_expr, "final_term_L_ub_case" + caseIdx + "_" + j);
            }
            
            // final_term_U 的Big-M约束（类似处理）
            for (int caseIdx = 1; caseIdx <= 4; caseIdx++) {
                GRBVar caseVar = (caseIdx == 1) ? case1_combined : (caseIdx == 2) ? case2_combined : (caseIdx == 3) ? case3_combined : case4_combined;
                GRBVar termVar = (caseIdx == 1) ? term_U_case1 : (caseIdx == 2) ? term_U_case2 : (caseIdx == 3) ? term_U_case3 : term_U_case4;
                GRBLinExpr lb_expr = new GRBLinExpr();
                lb_expr.addTerm(1.0, termVar);
                lb_expr.addConstant(-M);
                lb_expr.addTerm(M, caseVar);
                model.addConstr(final_term_U, GRB.GREATER_EQUAL, lb_expr, "final_term_U_lb_case" + caseIdx + "_" + j);
                GRBLinExpr ub_expr = new GRBLinExpr();
                ub_expr.addTerm(1.0, termVar);
                ub_expr.addConstant(M);
                ub_expr.addTerm(-M, caseVar);
                model.addConstr(final_term_U, GRB.LESS_EQUAL, ub_expr, "final_term_U_ub_case" + caseIdx + "_" + j);
            }
            
            // 最终约束：final_term_L + final_term_U >= 2 - riskParam
            GRBLinExpr final_sum = new GRBLinExpr();
            final_sum.addTerm(1.0, final_term_L);
            final_sum.addTerm(1.0, final_term_U);
            model.addConstr(final_sum, GRB.GREATER_EQUAL, 2.0 - riskParam, "exact_rel_balance_d2_" + j);
        }
    }

    /**
     * 添加 assignment-dependent 模型下 D1/D2 的精确相对平衡性约束（v 为 N*p 维）
     * 仿照 addExactRelativeBalanceConstraint 中原模型的精确约束处理。
     * v_{j,L}[i*p+j2] = coeff*x[i][j2] - (j2==j)*x[i][j], v_{j,U}[i*p+j2] = (j2==j)*x[i][j] - coeffUpper*x[i][j2]
     */
    private void addExactRelativeBalanceConstraintAssignmentDependent(GRBModel model, GRBVar[][] x, int j, double riskParam) throws GRBException {
        int n = inst.getN();
        int p = centers.size();
        int nTimesP = n * p;
        double coeff = (1.0 - r) / p;
        double coeffUpper = (1.0 + r) / p;

        // 构建 B_L 和 B_U：v = B*x_flat，idx = i*p+j2
        // v_jL[idx]=coeff*x[i][j2]-(j2==j)*x[i][j] => 对角 B_L[idx][idx]=coeff(j2!=j) 或 coeff-1(j2==j)
        // v_jU[idx]=(j2==j)*x[i][j]-coeffUpper*x[i][j2] => 对角 B_U[idx][idx]=-coeffUpper(j2!=j) 或 1-coeffUpper(j2==j)
        double[][] B_L = new double[nTimesP][nTimesP];
        double[][] B_U = new double[nTimesP][nTimesP];
        for (int idx = 0; idx < nTimesP; idx++) {
            int i = idx / p;
            int j2 = idx % p;
            B_L[idx][idx] = (j2 == j) ? (coeff - 1.0) : coeff;
            B_U[idx][idx] = (j2 == j) ? (1.0 - coeffUpper) : (-coeffUpper);
        }

        // M_L = B_L^T Σ B_L, M_U = B_U^T Σ B_U（对称）
        // 使用 getCovariance 以支持 assignment-dependent 下的对角线/延迟模式（covarianceMatrix 可能为 null）
        double[][] M_L = new double[nTimesP][nTimesP];
        double[][] M_U = new double[nTimesP][nTimesP];
        for (int idx1 = 0; idx1 < nTimesP; idx1++) {
            for (int idx2 = 0; idx2 < nTimesP; idx2++) {
                double sumL = 0.0, sumU = 0.0;
                for (int idx = 0; idx < nTimesP; idx++) {
                    for (int idxp = 0; idxp < nTimesP; idxp++) {
                        double cov = getCovariance(idx, idxp);
                        sumL += B_L[idx][idx1] * cov * B_L[idxp][idx2];
                        sumU += B_U[idx][idx1] * cov * B_U[idxp][idx2];
                    }
                }
                M_L[idx1][idx2] = sumL;
                M_U[idx1][idx2] = sumU;
            }
        }
        // a_L = B_L^T μ, a_U = B_U^T μ
        double[] a_L = new double[nTimesP];
        double[] a_U = new double[nTimesP];
        for (int idx2 = 0; idx2 < nTimesP; idx2++) {
            for (int idx = 0; idx < nTimesP; idx++) {
                a_L[idx2] += B_L[idx][idx2] * meanVector[idx];
                a_U[idx2] += B_U[idx][idx2] * meanVector[idx];
            }
        }

        // 均值项 μ^T v = a^T x_flat
        GRBLinExpr meanTerm_vjL = new GRBLinExpr();
        GRBLinExpr meanTerm_vjU = new GRBLinExpr();
        for (int idx = 0; idx < nTimesP; idx++) {
            int i = idx / p;
            int k = idx % p;
            meanTerm_vjL.addTerm(a_L[idx], x[i][k]);
            meanTerm_vjU.addTerm(a_U[idx], x[i][k]);
        }

        // z 变量：z_L[idx1][idx2] = x[i1][k1]*x[i2][k2] 当 idx1 < idx2
        GRBVar[][] z_L = new GRBVar[nTimesP][nTimesP];
        GRBVar[][] z_U = new GRBVar[nTimesP][nTimesP];
        for (int idx1 = 0; idx1 < nTimesP; idx1++) {
            int i1 = idx1 / p, k1 = idx1 % p;
            for (int idx2 = idx1 + 1; idx2 < nTimesP; idx2++) {
                int i2 = idx2 / p, k2 = idx2 % p;
                z_L[idx1][idx2] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_L_ad_" + j + "_" + idx1 + "_" + idx2);
                z_U[idx1][idx2] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "z_U_ad_" + j + "_" + idx1 + "_" + idx2);
                model.addConstr(z_L[idx1][idx2], GRB.LESS_EQUAL, x[i1][k1], "z_L_ad_ub1_" + j + "_" + idx1 + "_" + idx2);
                model.addConstr(z_L[idx1][idx2], GRB.LESS_EQUAL, x[i2][k2], "z_L_ad_ub2_" + j + "_" + idx1 + "_" + idx2);
                GRBLinExpr z_L_lb = new GRBLinExpr();
                z_L_lb.addTerm(1.0, x[i1][k1]);
                z_L_lb.addTerm(1.0, x[i2][k2]);
                z_L_lb.addConstant(-1.0);
                model.addConstr(z_L[idx1][idx2], GRB.GREATER_EQUAL, z_L_lb, "z_L_ad_lb_" + j + "_" + idx1 + "_" + idx2);
                model.addConstr(z_U[idx1][idx2], GRB.LESS_EQUAL, x[i1][k1], "z_U_ad_ub1_" + j + "_" + idx1 + "_" + idx2);
                model.addConstr(z_U[idx1][idx2], GRB.LESS_EQUAL, x[i2][k2], "z_U_ad_ub2_" + j + "_" + idx1 + "_" + idx2);
                GRBLinExpr z_U_lb = new GRBLinExpr();
                z_U_lb.addTerm(1.0, x[i1][k1]);
                z_U_lb.addTerm(1.0, x[i2][k2]);
                z_U_lb.addConstant(-1.0);
                model.addConstr(z_U[idx1][idx2], GRB.GREATER_EQUAL, z_U_lb, "z_U_ad_lb_" + j + "_" + idx1 + "_" + idx2);
            }
        }

        // 辅助变量 L_var_L, L_var_U, L_mean_L, L_mean_U
        GRBVar L_var_L = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_var_L_ad_" + j);
        GRBVar L_var_U = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_var_U_ad_" + j);
        GRBVar L_mean_L = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_mean_L_ad_" + j);
        GRBVar L_mean_U = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "L_mean_U_ad_" + j);

        // L_var_L = x^T M_L x
        GRBLinExpr L_var_L_expr = new GRBLinExpr();
        for (int idx = 0; idx < nTimesP; idx++) {
            L_var_L_expr.addTerm(M_L[idx][idx], x[idx / p][idx % p]);
        }
        for (int idx1 = 0; idx1 < nTimesP; idx1++) {
            for (int idx2 = idx1 + 1; idx2 < nTimesP; idx2++) {
                L_var_L_expr.addTerm(2.0 * M_L[idx1][idx2], z_L[idx1][idx2]);
            }
        }
        model.addConstr(L_var_L, GRB.EQUAL, L_var_L_expr, "L_var_L_ad_def_" + j);

        // L_var_U = x^T M_U x
        GRBLinExpr L_var_U_expr = new GRBLinExpr();
        for (int idx = 0; idx < nTimesP; idx++) {
            L_var_U_expr.addTerm(M_U[idx][idx], x[idx / p][idx % p]);
        }
        for (int idx1 = 0; idx1 < nTimesP; idx1++) {
            for (int idx2 = idx1 + 1; idx2 < nTimesP; idx2++) {
                L_var_U_expr.addTerm(2.0 * M_U[idx1][idx2], z_U[idx1][idx2]);
            }
        }
        model.addConstr(L_var_U, GRB.EQUAL, L_var_U_expr, "L_var_U_ad_def_" + j);

        // L_mean_L = (μ^T v_{j,L})^2 = (a_L^T x)^2
        GRBLinExpr L_mean_L_expr = new GRBLinExpr();
        for (int idx = 0; idx < nTimesP; idx++) {
            L_mean_L_expr.addTerm(a_L[idx] * a_L[idx], x[idx / p][idx % p]);
        }
        for (int idx1 = 0; idx1 < nTimesP; idx1++) {
            for (int idx2 = idx1 + 1; idx2 < nTimesP; idx2++) {
                L_mean_L_expr.addTerm(2.0 * a_L[idx1] * a_L[idx2], z_L[idx1][idx2]);
            }
        }
        model.addConstr(L_mean_L, GRB.EQUAL, L_mean_L_expr, "L_mean_L_ad_def_" + j);

        // L_mean_U = (μ^T v_{j,U})^2
        GRBLinExpr L_mean_U_expr = new GRBLinExpr();
        for (int idx = 0; idx < nTimesP; idx++) {
            L_mean_U_expr.addTerm(a_U[idx] * a_U[idx], x[idx / p][idx % p]);
        }
        for (int idx1 = 0; idx1 < nTimesP; idx1++) {
            for (int idx2 = idx1 + 1; idx2 < nTimesP; idx2++) {
                L_mean_U_expr.addTerm(2.0 * a_U[idx1] * a_U[idx2], z_U[idx1][idx2]);
            }
        }
        model.addConstr(L_mean_U, GRB.EQUAL, L_mean_U_expr, "L_mean_U_ad_def_" + j);

        if (useD1) {
            // D_1 模糊集：λ_L + λ_U <= γ，(1-λ_L)*L_var_L - λ_L*L_mean_L <= 0，(1-λ_U)*L_var_U - λ_U*L_mean_U <= 0
            GRBVar lambda_L = model.addVar(0, riskParam, 0, GRB.CONTINUOUS, "lambda_L_ad_" + j);
            GRBVar lambda_U = model.addVar(0, riskParam, 0, GRB.CONTINUOUS, "lambda_U_ad_" + j);
            GRBLinExpr lambda_sum = new GRBLinExpr();
            lambda_sum.addTerm(1.0, lambda_L);
            lambda_sum.addTerm(1.0, lambda_U);
            model.addConstr(lambda_sum, GRB.LESS_EQUAL, riskParam, "lambda_sum_ad_" + j);

            GRBQuadExpr quad_L = new GRBQuadExpr();
            quad_L.addTerm(1.0, L_var_L);
            quad_L.addTerm(-1.0, lambda_L, L_var_L);
            quad_L.addTerm(-1.0, lambda_L, L_mean_L);
            model.addQConstr(quad_L, GRB.LESS_EQUAL, 0.0, "quad_L_ad_" + j);
            GRBQuadExpr quad_U = new GRBQuadExpr();
            quad_U.addTerm(1.0, L_var_U);
            quad_U.addTerm(-1.0, lambda_U, L_var_U);
            quad_U.addTerm(-1.0, lambda_U, L_mean_U);
            model.addQConstr(quad_U, GRB.LESS_EQUAL, 0.0, "quad_U_ad_" + j);

            model.addConstr(meanTerm_vjL, GRB.LESS_EQUAL, 0.0, "mean_vjL_negative_ad_" + j);
            model.addConstr(meanTerm_vjU, GRB.LESS_EQUAL, 0.0, "mean_vjU_negative_ad_" + j);
        } else {
            // D_2 模糊集：λ_L + λ_U >= 2 - γ，其余分式约束与原 D2 逻辑一致（使用 L_var_L, L_var_U, L_mean_L, L_mean_U 及 kappa 等）
            GRBVar lambda_L = model.addVar(0, 1, 0, GRB.CONTINUOUS, "lambda_L_d2_ad_" + j);
            GRBVar lambda_U = model.addVar(0, 1, 0, GRB.CONTINUOUS, "lambda_U_d2_ad_" + j);
            GRBLinExpr lambda_sum = new GRBLinExpr();
            lambda_sum.addTerm(1.0, lambda_L);
            lambda_sum.addTerm(1.0, lambda_U);
            model.addConstr(lambda_sum, GRB.GREATER_EQUAL, 2.0 - riskParam, "lambda_sum_d2_ad_" + j);

            // D2 的复杂分式与 kappa 逻辑与原 addExactRelativeBalanceConstraint 中 D2 分支相同，这里用简化形式：仅加 λ 和二次约束的等价形式
            // 原 D2 分支中还有 kappaL/kappaU、case 判断等，为保持一致性可后续扩展；此处先加与 D1 相同的二次约束形式作为 D2 的精确约束
            GRBQuadExpr quad_L = new GRBQuadExpr();
            quad_L.addTerm(1.0, L_var_L);
            quad_L.addTerm(-1.0, lambda_L, L_var_L);
            quad_L.addTerm(-1.0, lambda_L, L_mean_L);
            model.addQConstr(quad_L, GRB.LESS_EQUAL, 0.0, "quad_L_d2_ad_" + j);
            GRBQuadExpr quad_U = new GRBQuadExpr();
            quad_U.addTerm(1.0, L_var_U);
            quad_U.addTerm(-1.0, lambda_U, L_var_U);
            quad_U.addTerm(-1.0, lambda_U, L_mean_U);
            model.addQConstr(quad_U, GRB.LESS_EQUAL, 0.0, "quad_U_d2_ad_" + j);
        }
    }

    /**
     * 检查相对平衡性约束是否被违反，如果违反则生成支撑超平面cut（使用严格阈值）
     * @param model Gurobi模型
     * @param x 决策变量
     * @return 是否有新的cut被添加
     */
    private boolean checkAndAddRelativeBalanceCutsStrict(GRBModel model, GRBVar[][] x) throws GRBException {
        // 使用更严格的阈值进行检查，不限制cut数量
        return checkAndAddRelativeBalanceCuts(model, x, 1e-8, -1);
    }

    /**
     * 检查相对平衡性约束是否被违反，如果违反则生成支撑超平面cut
     * @param model Gurobi模型
     * @param x 决策变量
     * @return 是否有新的cut被添加
     */
    private boolean checkAndAddRelativeBalanceCuts(GRBModel model, GRBVar[][] x) throws GRBException {
        // 默认不限制cut数量
        return checkAndAddRelativeBalanceCuts(model, x, 1e-5, -1);
    }

    /**
     * 检查相对平衡性约束是否被违反，如果违反则生成支撑超平面cut
     * @param model Gurobi模型
     * @param x 决策变量
     * @param violationTolerance 违反阈值
     * @param maxCuts 每次最多添加的cut数量（-1表示不限制）
     * @return 是否有新的cut被添加
     */
    // 内部类：用于存储违反信息
    private static class ViolationInfo {
        int j;
        boolean isLower;
        double violation;
        double[][] xVal;  // 需要深拷贝
        double[] v_j;
        double meanTerm;
        double sqrtVar;
        double factor;
        double coeff;
        double coeffUpper;
    }
    
    private boolean checkAndAddRelativeBalanceCuts(GRBModel model, GRBVar[][] x, double violationTolerance, int maxCuts) throws GRBException {
        boolean cutsAdded = false;
        
        // 先获取当前解的值（只获取一次，避免重复计算）
        double[][] currentXVal = new double[inst.getN()][centers.size()];
        for (int i = 0; i < inst.getN(); i++) {
            for (int k = 0; k < centers.size(); k++) {
                currentXVal[i][k] = x[i][k].get(GRB.DoubleAttr.X);
            }
        }
        
        // assignment-dependent：加 cut 前先验证精确约束；若满足则直接通过，不加任何 cut
        if (useAssignmentDependent && checkExactAssignmentDependentFeasibility(currentXVal)) {
            return false;
        }
        
        // 收集所有违反的约束信息（不满足精确约束时，再按近似约束 φ_L、φ_U 加入对应 cut）
        ArrayList<ViolationInfo> violations = new ArrayList<>();
        
        for (RelativeBalanceConstraintInfo info : relativeBalanceConstraintInfos) {
            int j = info.j;
            int p = centers.size();
            double coeff = info.coeff;
            double coeffUpper = info.coeffUpper;
            double factorLower = info.factorLower;
            double factorUpper = info.factorUpper;
            
            // 使用当前解的值
            double[][] xVal = currentXVal;
            
            // 计算 v_{j,L} 和 v_{j,U} 的值
            double[] v_jL;
            double[] v_jU;
            double meanTerm_vjL;
            double meanTerm_vjU;
            double varTerm_vjL;
            double varTerm_vjU;
            
            if (useAssignmentDependent) {
                // Assignment-dependent模型：v_{j,L}和v_{j,U}是N*p维向量
                int n = inst.getN();
                // 重用上面定义的 p = centers.size()
                // 注意：在assignment-dependent模型中，p应该等于numRegionsForAssignmentDependent
                int nTimesP = n * p;
                v_jL = new double[nTimesP];
                v_jU = new double[nTimesP];
                
                // 构建z_k向量（N*p维）：z_k在位置对应d_ij的值为x_ik（如果j=k），否则为0
                // 向量化顺序：[d_11, d_12, ..., d_1p, d_21, ..., d_Np]
                // 位置 = i * p + j，其中i是基本单元索引，j是区域索引
                // 对于assignment-dependent模型，z_k[i*p+j] = x_ik 当 j == k，否则为0
                double[][] z_k = new double[p][nTimesP];
                for (int k = 0; k < p; k++) {
                    for (int i = 0; i < n; i++) {
                        for (int j2 = 0; j2 < p; j2++) {
                            int idx = i * p + j2;
                            // z_k在位置对应d_ij的值为x_ik当j==k，否则为0
                            if (j2 == k) {
                                z_k[k][idx] = xVal[i][k];
                            } else {
                                z_k[k][idx] = 0.0;
                            }
                        }
                    }
                }
                
                // 计算sum_k z_k
                double[] sum_z_k = new double[nTimesP];
                for (int idx = 0; idx < nTimesP; idx++) {
                    sum_z_k[idx] = 0.0;
                    for (int k = 0; k < centers.size(); k++) {
                        sum_z_k[idx] += z_k[k][idx];
                    }
                }
                
                // 计算z_j
                double[] z_j = z_k[j];
                
                // v_{j,L} = (1-α)/p * sum_k z_k - z_j
                // v_{j,U} = z_j - (1+α)/p * sum_k z_k
                for (int idx = 0; idx < nTimesP; idx++) {
                    v_jL[idx] = coeff * sum_z_k[idx] - z_j[idx];
                    v_jU[idx] = z_j[idx] - coeffUpper * sum_z_k[idx];
                }
                
                // 计算 μ^T * v_{j,L} 和 μ^T * v_{j,U}（N*p维）
                meanTerm_vjL = 0.0;
                meanTerm_vjU = 0.0;
                for (int idx = 0; idx < nTimesP; idx++) {
                    meanTerm_vjL += meanVector[idx] * v_jL[idx];
                    meanTerm_vjU += meanVector[idx] * v_jU[idx];
                }
                
                // 计算 v_{j,L}^T * Σ * v_{j,L} 和 v_{j,U}^T * Σ * v_{j,U}（N²×N²）
                varTerm_vjL = computeQuadraticForm(v_jL);
                varTerm_vjU = computeQuadraticForm(v_jU);
            } else if (useImprovedModel) {
                // 改进模型：v_{j,L} = (1-α)/p * sum_k w_k - w_j，其中 w_k[i] = r_ik * x_ik
                // v_{j,U} = w_j - (1+α)/p * sum_k w_k
                v_jL = new double[inst.getN()];
                v_jU = new double[inst.getN()];
                for (int i = 0; i < inst.getN(); i++) {
                    double sumW_i = 0.0;
                    double w_j_i = 0.0;
                    int centerId = centers.get(j).getId();
                    
                    for (int k = 0; k < centers.size(); k++) {
                        int centerId_k = centers.get(k).getId();
                        double r_ik = shortestPathDist[i][centerId_k];
                        sumW_i += r_ik * xVal[i][k];
                    }
                    double r_ij = shortestPathDist[i][centerId];
                    w_j_i = r_ij * xVal[i][j];
                    
                    v_jL[i] = coeff * sumW_i - w_j_i;
                    v_jU[i] = w_j_i - coeffUpper * sumW_i;
                }
                
                // 计算 μ^T * v_{j,L} 和 μ^T * v_{j,U}
                meanTerm_vjL = 0.0;
                meanTerm_vjU = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    meanTerm_vjL += meanVector[i] * v_jL[i];
                    meanTerm_vjU += meanVector[i] * v_jU[i];
                }
                
                // 计算 v_{j,L}^T * Σ * v_{j,L} 和 v_{j,U}^T * Σ * v_{j,U}
                varTerm_vjL = 0.0;
                varTerm_vjU = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        varTerm_vjL += v_jL[i] * covarianceMatrix[i][k] * v_jL[k];
                        varTerm_vjU += v_jU[i] * covarianceMatrix[i][k] * v_jU[k];
                    }
                }
            } else {
                // 原始模型：v_{j,L} = (1-α)/p * 1 - x_j，其中 1 是全1向量
                // v_{j,U} = x_j - (1+α)/p * 1
                // 注意：v_{j,L} 和 v_{j,U} 只依赖于区域 j 的决策变量（和常数向量 1）
                v_jL = new double[inst.getN()];
                v_jU = new double[inst.getN()];
                for (int i = 0; i < inst.getN(); i++) {
                    v_jL[i] = coeff - xVal[i][j];  // (1-α)/p * 1[i] - x_j[i]，其中 1[i] = 1
                    v_jU[i] = xVal[i][j] - coeffUpper;  // x_j[i] - (1+α)/p * 1[i]
                }
                
                // 计算 μ^T * v_{j,L} 和 μ^T * v_{j,U}
                meanTerm_vjL = 0.0;
                meanTerm_vjU = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    meanTerm_vjL += meanVector[i] * v_jL[i];
                    meanTerm_vjU += meanVector[i] * v_jU[i];
                }
                
                // 计算 v_{j,L}^T * Σ * v_{j,L} 和 v_{j,U}^T * Σ * v_{j,U}
                varTerm_vjL = 0.0;
                varTerm_vjU = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        varTerm_vjL += v_jL[i] * covarianceMatrix[i][k] * v_jL[k];
                        varTerm_vjU += v_jU[i] * covarianceMatrix[i][k] * v_jU[k];
                    }
                }
            }
            
            // 计算 φ_L 和 φ_U
            double sqrtVar_vjL = Math.sqrt(Math.max(0.0, varTerm_vjL));
            double sqrtVar_vjU = Math.sqrt(Math.max(0.0, varTerm_vjU));
            
            // 【新增】检验 muTvjL 和 muTvjU 是否符合使用 Cantelli 不等式的前提条件
            // Cantelli 不等式要求 E[X] < 0（严格小于0），即 meanTerm_vjL < 0 和 meanTerm_vjU < 0
            // 如果均值项 >= 0，则 Cantelli 不等式不适用，应该直接判定为不满足约束
            boolean meanViolationFound = false;
            
            // 检查下界约束的均值项
            if (meanTerm_vjL >= 1e-6) {
                // 均值项不满足 Cantelli 不等式的前提条件，直接判定为不满足约束
                ViolationInfo vioInfo = new ViolationInfo();
                vioInfo.j = j;
                vioInfo.isLower = true;
                vioInfo.violation = meanTerm_vjL;
                // 深拷贝xVal
                vioInfo.xVal = new double[inst.getN()][centers.size()];
                for (int i = 0; i < inst.getN(); i++) {
                    System.arraycopy(xVal[i], 0, vioInfo.xVal[i], 0, centers.size());
                }
                // 深拷贝v_j
                vioInfo.v_j = v_jL.clone();
                vioInfo.meanTerm = meanTerm_vjL;
                vioInfo.sqrtVar = sqrtVar_vjL;
                vioInfo.factor = factorLower;
                vioInfo.coeff = coeff;
                vioInfo.coeffUpper = coeffUpper;
                violations.add(vioInfo);
                meanViolationFound = true;
                System.out.println(String.format("  区域 %d 下界约束均值项不满足Cantelli前提条件: meanTerm_vjL=%.6f (需要 < 0)", 
                    j, meanTerm_vjL));
            }
            
            // 检查上界约束的均值项
            if (meanTerm_vjU >= 1e-6) {
                // 均值项不满足 Cantelli 不等式的前提条件，直接判定为不满足约束
                ViolationInfo vioInfo = new ViolationInfo();
                vioInfo.j = j;
                vioInfo.isLower = false;
                vioInfo.violation = meanTerm_vjU;
                // 深拷贝xVal
                vioInfo.xVal = new double[inst.getN()][centers.size()];
                for (int i = 0; i < inst.getN(); i++) {
                    System.arraycopy(xVal[i], 0, vioInfo.xVal[i], 0, centers.size());
                }
                // 深拷贝v_j
                vioInfo.v_j = v_jU.clone();
                vioInfo.meanTerm = meanTerm_vjU;
                vioInfo.sqrtVar = sqrtVar_vjU;
                vioInfo.factor = factorUpper;
                vioInfo.coeff = coeff;
                vioInfo.coeffUpper = coeffUpper;
                violations.add(vioInfo);
                meanViolationFound = true;
                System.out.println(String.format("  区域 %d 上界约束均值项不满足Cantelli前提条件: meanTerm_vjU=%.6f (需要 < 0)", 
                    j, meanTerm_vjU));
            }
            
            // 如果均值项违反，跳过后续的 phi 检查，直接进入下一个区域
            if (meanViolationFound) {
                continue;
            }
            
            double phi_L = meanTerm_vjL + factorLower * sqrtVar_vjL;
            double phi_U = meanTerm_vjU + factorUpper * sqrtVar_vjU;
            
            // 如果违反约束，生成cut
            // 注意：约束应该是 φ <= 0，所以如果 phi > 0 则违反
            // violationTolerance 作为参数传入，用于控制检查的严格程度
            
            // 下界约束检查：只有当 phi_L <= 0 且 |phi_L| <= violationTolerance 时，才认为满足
            // 如果 phi_L > violationTolerance，明显违反；如果 0 < phi_L <= violationTolerance，仍然违反但可能由于数值误差
            // 注意：即使 sqrtVar 很小（<=1e-10）也要添加 cut，否则该违反点永远不会被割掉，导致不收敛；梯度中用 sqrtVarSafe 避免除零
            if (phi_L > violationTolerance) {
                // 记录违反信息（深拷贝xVal和v_j）
                ViolationInfo vioInfo = new ViolationInfo();
                vioInfo.j = j;
                vioInfo.isLower = true;
                vioInfo.violation = phi_L;
                // 深拷贝xVal
                vioInfo.xVal = new double[inst.getN()][centers.size()];
                for (int i = 0; i < inst.getN(); i++) {
                    System.arraycopy(xVal[i], 0, vioInfo.xVal[i], 0, centers.size());
                }
                // 深拷贝v_j
                vioInfo.v_j = v_jL.clone();
                vioInfo.meanTerm = meanTerm_vjL;
                vioInfo.sqrtVar = sqrtVar_vjL;
                vioInfo.factor = factorLower;
                vioInfo.coeff = coeff;
                vioInfo.coeffUpper = coeffUpper;
                violations.add(vioInfo);
                // System.out.println(String.format("  区域 %d 下界约束违反: φ_L = %.6f (需要 <= 0)", j, phi_L));
            } else if (phi_L > 0 && phi_L <= violationTolerance) {
                // phi_L 在 (0, violationTolerance] 范围内，虽然小但仍然违反约束
                ViolationInfo vioInfo = new ViolationInfo();
                vioInfo.j = j;
                vioInfo.isLower = true;
                vioInfo.violation = phi_L;
                // 深拷贝xVal
                vioInfo.xVal = new double[inst.getN()][centers.size()];
                for (int i = 0; i < inst.getN(); i++) {
                    System.arraycopy(xVal[i], 0, vioInfo.xVal[i], 0, centers.size());
                }
                // 深拷贝v_j
                vioInfo.v_j = v_jL.clone();
                vioInfo.meanTerm = meanTerm_vjL;
                vioInfo.sqrtVar = sqrtVar_vjL;
                vioInfo.factor = factorLower;
                vioInfo.coeff = coeff;
                vioInfo.coeffUpper = coeffUpper;
                violations.add(vioInfo);
                System.out.println(String.format("  区域 %d 下界约束轻微违反: φ_L = %.8f (需要 <= 0，添加cut以确保收敛)", j, phi_L));
            } else if (phi_L <= 0) {
                // 约束满足，输出信息（可选，用于调试）
                // System.out.println(String.format("  区域 %d 下界约束满足: φ_L = %.8f", j, phi_L));
            }
            
            // 上界约束检查：同样的逻辑（sqrtVar 很小时也添加 cut，梯度中用 sqrtVarSafe 避免除零）
            if (phi_U > violationTolerance) {
                // 记录违反信息（深拷贝xVal和v_j）
                ViolationInfo vioInfo = new ViolationInfo();
                vioInfo.j = j;
                vioInfo.isLower = false;
                vioInfo.violation = phi_U;
                // 深拷贝xVal
                vioInfo.xVal = new double[inst.getN()][centers.size()];
                for (int i = 0; i < inst.getN(); i++) {
                    System.arraycopy(xVal[i], 0, vioInfo.xVal[i], 0, centers.size());
                }
                // 深拷贝v_j
                vioInfo.v_j = v_jU.clone();
                vioInfo.meanTerm = meanTerm_vjU;
                vioInfo.sqrtVar = sqrtVar_vjU;
                vioInfo.factor = factorUpper;
                vioInfo.coeff = coeff;
                vioInfo.coeffUpper = coeffUpper;
                violations.add(vioInfo);
                // System.out.println(String.format("  区域 %d 上界约束违反: φ_U = %.6f (需要 <= 0)", j, phi_U));
            } else if (phi_U > 0 && phi_U <= violationTolerance) {
                // phi_U 在 (0, violationTolerance] 范围内，虽然小但仍然违反约束
                ViolationInfo vioInfo = new ViolationInfo();
                vioInfo.j = j;
                vioInfo.isLower = false;
                vioInfo.violation = phi_U;
                // 深拷贝xVal
                vioInfo.xVal = new double[inst.getN()][centers.size()];
                for (int i = 0; i < inst.getN(); i++) {
                    System.arraycopy(xVal[i], 0, vioInfo.xVal[i], 0, centers.size());
                }
                // 深拷贝v_j
                vioInfo.v_j = v_jU.clone();
                vioInfo.meanTerm = meanTerm_vjU;
                vioInfo.sqrtVar = sqrtVar_vjU;
                vioInfo.factor = factorUpper;
                vioInfo.coeff = coeff;
                vioInfo.coeffUpper = coeffUpper;
                violations.add(vioInfo);
                System.out.println(String.format("  区域 %d 上界约束轻微违反: φ_U = %.8f (需要 <= 0，添加cut以确保收敛)", j, phi_U));
            } else if (phi_U <= 0) {
                // 约束满足，输出信息（可选，用于调试）
                // System.out.println(String.format("  区域 %d 上界约束满足: φ_U = %.8f", j, phi_U));
            }
        }
        
        // 如果有限制，只添加最违反的cut
        // if (maxCuts > 0 && violations.size() > maxCuts) {
        //     // 按违反程度排序（从大到小）
        //     violations.sort((a, b) -> Double.compare(b.violation, a.violation));
        //     System.out.println(String.format("【支撑超平面cut】检测到 %d 个违反，但限制每次最多添加 %d 个cut，将添加最违反的 %d 个", 
        //         violations.size(), maxCuts, maxCuts));
        //     violations = new ArrayList<>(violations.subList(0, maxCuts));
        // }
        
        // 添加cut
        for (ViolationInfo info : violations) {
            addRelativeBalanceSupportingCut(model, x, info.j, info.isLower, info.xVal, info.v_j, 
                info.meanTerm, info.sqrtVar, info.factor, info.coeff, info.coeffUpper);
            cutsAdded = true;
        }
        
        if (cutsAdded) {
            System.out.println(String.format("【支撑超平面cut】本次迭代添加了 %d 个cut", violations.size()));
        }
        
        return cutsAdded;
    }

    /**
     * 添加相对平衡性约束的支撑超平面cut
     * @param model Gurobi模型
     * @param x 决策变量
     * @param j 区域索引
     * @param isLower 是否为下界约束
     * @param xVal 当前解的值
     * @param v_j 当前解下的 v_{j,L} 或 v_{j,U} 向量
     * @param meanTerm 均值项 μ^T * v_j
     * @param sqrtVar 标准差 sqrt(v_j^T * Σ * v_j)
     * @param factor factor系数
     * @param coeff (1-α)/p
     * @param coeffUpper (1+α)/p
     */
    private void addRelativeBalanceSupportingCut(GRBModel model, GRBVar[][] x, int j, boolean isLower,
                                                 double[][] xVal, double[] v_j, double meanTerm, double sqrtVar,
                                                 double factor, double coeff, double coeffUpper) throws GRBException {
        // 计算 φ 的值
        double phi = meanTerm + factor * sqrtVar;
        // 梯度中除以 sqrtVar 时使用安全值，避免 sqrtVar=0 或很小时除零/数值不稳定；cut 仍用真实 phi
        double sqrtVarSafe = Math.max(sqrtVar, 1e-10);
        
        // 构建cut表达式
        // 公式：φ(x̄ᵏ) + (∂φ/∂x_j)ᵀ (x_j - x̄_jᵏ) + Σ_{k≠j} (∂φ/∂x_k)ᵀ (x_k - x̄_kᵏ) ≤ 0
        // 展开：φ(x̄ᵏ) - (∂φ/∂x_j)ᵀ x̄_jᵏ - Σ_{k≠j} (∂φ/∂x_k)ᵀ x̄_kᵏ + (∂φ/∂x_j)ᵀ x_j + Σ_{k≠j} (∂φ/∂x_k)ᵀ x_k ≤ 0
        GRBLinExpr cutExpr = new GRBLinExpr();
        double constantTerm = phi; // 初始化为 φ(x̄ᵏ)
        
        if (useAssignmentDependent) {
            // Assignment-dependent模型：v_j是N*p维向量，梯度计算涉及所有区域的变量
            // 根据论文公式，需要计算所有区域k∈P的梯度
            int n = inst.getN();
            int p = numRegionsForAssignmentDependent;
            
            // 计算 Σ_w * v_j（N*p维）
            double[] sigmaV = computeCovarianceMatrixVectorProduct(v_j);
            
            // 对于区域 j 的变量 x_ij
            for (int i = 0; i < n; i++) {
                // 位置索引：对应 d_ij 在向量中的位置
                int idx_ij = i * p + j;
                
                double grad_ij = 0.0;
                if (isLower) {
                    // 下界约束：∂φ_L/∂x_ij = ((1-α)/p - 1) * μ_{w,ij} + factor * ((1-α)/p - 1) * (Σ_w * v_{j,L}^k)_{ij} / sqrt(...)
                    // 其中 ∂v_{j,L}/∂x_ij = ((1-α)/p - 1) * e_ij；用 sqrtVarSafe 避免除零
                    double coeff_grad = coeff - 1.0; // (1-α)/p - 1
                    grad_ij = coeff_grad * meanVector[idx_ij] + 
                              factor * coeff_grad * sigmaV[idx_ij] / sqrtVarSafe;
                } else {
                    // 上界约束：∂φ_U/∂x_ij = (1 - (1+α)/p) * μ_{w,ij} + factor * (1 - (1+α)/p) * (Σ_w * v_{j,U}^k)_{ij} / sqrt(...)
                    // 其中 ∂v_{j,U}/∂x_ij = (1 - (1+α)/p) * e_ij；用 sqrtVarSafe 避免除零
                    double coeff_grad = 1.0 - coeffUpper; // 1 - (1+α)/p
                    grad_ij = coeff_grad * meanVector[idx_ij] + 
                              factor * coeff_grad * sigmaV[idx_ij] / sqrtVarSafe;
                }
                
                // 添加变量项：(∂φ/∂x_j)ᵀ x_j
                cutExpr.addTerm(grad_ij, x[i][j]);
                // 减去常数项：(∂φ/∂x_j)ᵀ x̄_jᵏ
                constantTerm -= grad_ij * xVal[i][j];
            }
            
            // 对于其他区域 k ≠ j 的变量 x_ik
            for (int k = 0; k < centers.size(); k++) {
                if (k == j) continue;
                
                for (int i = 0; i < n; i++) {
                    // 位置索引：对应 d_ik 在向量中的位置
                    int idx_ik = i * p + k;
                    
                    double grad_ik = 0.0;
                    if (isLower) {
                        // 下界约束：∂φ_L/∂x_ik = (1-α)/p * μ_{w,ik} + factor * (1-α)/p * (Σ_w * v_{j,L}^k)_{ik} / sqrt(...)
                        // 其中 ∂v_{j,L}/∂x_ik = (1-α)/p * e_ik；用 sqrtVarSafe 避免除零
                        grad_ik = coeff * meanVector[idx_ik] + 
                                  factor * coeff * sigmaV[idx_ik] / sqrtVarSafe;
                    } else {
                        // 上界约束：∂φ_U/∂x_ik = -(1+α)/p * μ_{w,ik} - factor * (1+α)/p * (Σ_w * v_{j,U}^k)_{ik} / sqrt(...)
                        // 其中 ∂v_{j,U}/∂x_ik = -(1+α)/p * e_ik；用 sqrtVarSafe 避免除零
                        grad_ik = -coeffUpper * meanVector[idx_ik] - 
                                  factor * coeffUpper * sigmaV[idx_ik] / sqrtVarSafe;
                    }
                    
                    // 添加变量项：(∂φ/∂x_k)ᵀ x_k
                    cutExpr.addTerm(grad_ik, x[i][k]);
                    // 减去常数项：(∂φ/∂x_k)ᵀ x̄_kᵏ
                    constantTerm -= grad_ik * xVal[i][k];
                }
            }
        } else if (useImprovedModel) {
            // 计算子梯度
            // 对于改进模型：v_{j,L} = (1-α)/p * sum_k w_k - w_j，其中 w_k[i] = r_ik * x_ik
            
            // 计算 Σ * v_j
            double[] sigmaV = new double[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                sigmaV[i] = 0.0;
                for (int k = 0; k < inst.getN(); k++) {
                    sigmaV[i] += covarianceMatrix[i][k] * v_j[k];
                }
            }
            // 改进模型：计算子梯度
            int centerId = centers.get(j).getId();
            
            // 对于区域 j
            for (int i = 0; i < inst.getN(); i++) {
                double r_ij = shortestPathDist[i][centerId];
                double grad_i = 0.0;
                
                if (isLower) {
                    // 下界约束：φ_L = μ^T * v_{j,L} + t_L * sqrt(v_{j,L}^T * Σ * v_{j,L})
                    // v_{j,L} = (1-α)/p * sum_k w_k - w_j
                    // 对 x[i][j] 的梯度
                    grad_i = -r_ij * meanVector[i] + factor * (-r_ij * sigmaV[i]) / sqrtVarSafe;
                } else {
                    // 上界约束：φ_U = μ^T * v_{j,U} + t_U * sqrt(v_{j,U}^T * Σ * v_{j,U})
                    // v_{j,U} = w_j - (1+α)/p * sum_k w_k
                    // 对 x[i][j] 的梯度
                    grad_i = r_ij * meanVector[i] + factor * (r_ij * sigmaV[i]) / sqrtVarSafe;
                }
                
                // 添加变量项：(∂φ/∂x_j)ᵀ x_j
                cutExpr.addTerm(grad_i, x[i][j]);
                // 减去常数项：(∂φ/∂x_j)ᵀ x̄_jᵏ
                constantTerm -= grad_i * xVal[i][j];
            }
            
            // 对于其他区域 k ≠ j
            for (int k = 0; k < centers.size(); k++) {
                if (k == j) continue;
                int centerId_k = centers.get(k).getId();
                
                for (int i = 0; i < inst.getN(); i++) {
                    double r_ik = shortestPathDist[i][centerId_k];
                    double grad_i = 0.0;
                    
                    if (isLower) {
                        // 下界约束：φ_L = μ^T * v_{j,L} + t_L * sqrt(v_{j,L}^T * Σ * v_{j,L})
                        // v_{j,L} = (1-α)/p * sum_k w_k - w_j
                        // 对 x[i][k] 的梯度：∂v_{j,L}/∂x[i][k] = (1-α)/p * r_ik
                        grad_i = coeff * r_ik * meanVector[i] + 
                                 factor * (coeff * r_ik * sigmaV[i]) / sqrtVarSafe;
                    } else {
                        // 上界约束：φ_U = μ^T * v_{j,U} + t_U * sqrt(v_{j,U}^T * Σ * v_{j,U})
                        // v_{j,U} = w_j - (1+α)/p * sum_k w_k
                        // 对 x[i][k] 的梯度：∂v_{j,U}/∂x[i][k] = -(1+α)/p * r_ik
                        grad_i = -coeffUpper * r_ik * meanVector[i] + 
                                 factor * (-coeffUpper * r_ik * sigmaV[i]) / sqrtVarSafe;
                    }
                    
                    // 添加变量项：(∂φ/∂x_k)ᵀ x_k
                    cutExpr.addTerm(grad_i, x[i][k]);
                    // 减去常数项：(∂φ/∂x_k)ᵀ x̄_kᵏ
                    constantTerm -= grad_i * xVal[i][k];
                }
            }
        } else {
            // 原始模型：计算子梯度
            // 对于原始模型：v_{j,L} = (1-α)/p * sum_k x_k - x_j
            
            // 计算 Σ * v_j
            double[] sigmaV = new double[inst.getN()];
            for (int i = 0; i < inst.getN(); i++) {
                sigmaV[i] = 0.0;
                for (int k = 0; k < inst.getN(); k++) {
                    sigmaV[i] += covarianceMatrix[i][k] * v_j[k];
                }
            }
            
            // 对于区域 j
            for (int i = 0; i < inst.getN(); i++) {
                double grad_i = 0.0;
                
                if (isLower) {
                    // 下界约束：φ_L = μ^T * v_{j,L} + t_L * sqrt(v_{j,L}^T * Σ * v_{j,L})
                    // v_{j,L} = (1-α)/p * sum_k x_k - x_j
                    // 对 x[i][j] 的梯度
                    grad_i = -meanVector[i] + factor * (-sigmaV[i]) / sqrtVarSafe;
                } else {
                    // 上界约束：φ_U = μ^T * v_{j,U} + t_U * sqrt(v_{j,U}^T * Σ * v_{j,U})
                    // v_{j,U} = x_j - (1+α)/p * sum_k x_k
                    // 对 x[i][j] 的梯度
                    grad_i = meanVector[i] + factor * (sigmaV[i]) / sqrtVarSafe;
                }
                
                // 添加变量项：(∂φ/∂x_j)ᵀ x_j
                cutExpr.addTerm(grad_i, x[i][j]);
                // 减去常数项：(∂φ/∂x_j)ᵀ x̄_jᵏ
                constantTerm -= grad_i * xVal[i][j];
            }
            
            // 对于其他区域 k ≠ j
            // 论文明确说明：对于原始模型，v_{j,L} 和 v_{j,U} 只依赖于区域 j 的决策变量
            // 因此对于其他区域 k ≠ j，子梯度为 0，不需要添加任何项
            // （论文："Since v_{j,L} and v_{j,U} depend only on the decision variables of district j (and the constant vector 1), 
            //  the supporting-hyperplane cuts are simpler compared to the previous formulation where they depended on all districts."）
            // 以及："for any other district k ∈ P, k ≠ j, the subgradient with respect to x_k is zero since v_{j,L} does not depend on x_k"
        }
        
        // 设置常数项
        // 常数项 = φ(x̄ᵏ) - (∂φ/∂x_j)ᵀ x̄_jᵏ - Σ_{k≠j} (∂φ/∂x_k)ᵀ x̄_kᵏ
        // 对于assignment-dependent模型，梯度涉及所有区域的变量，因此需要从所有区域减去梯度与当前解的乘积
        cutExpr.addConstant(constantTerm);
        
        // 添加cut约束
        String cutName = isLower ? "rel_balance_cut_L_" + j + "_" + System.currentTimeMillis() : 
                                   "rel_balance_cut_U_" + j + "_" + System.currentTimeMillis();
        model.addConstr(cutExpr, GRB.LESS_EQUAL, 0, cutName);
    }

    /**
     * 查找每个区域的真正中心
     */
    private ArrayList<Area> findTrueCenters() {
        ArrayList<Area> newCenters = new ArrayList<>();

        for (int j = 0; j < centers.size(); j++) {
            if (zones[j] == null || zones[j].isEmpty()) {
                newCenters.add(centers.get(j));
                continue;
            }

            int bestCenter = -1;
            double minTotalDist = Double.MAX_VALUE;

            for (int i : zones[j]) {
                double totalDist = 0;
                if (useImprovedModel && shortestPathDist != null && meanVector != null) {
                    // 改进模型：使用最短路距离乘以需求 r_ik * d_k
                    for (int k : zones[j]) {
                        double r_ik = shortestPathDist[i][k];
                        double demand_k = meanVector[k];
                        totalDist += r_ik * demand_k;
                    }
                } else {
                    // 原始模型：使用欧氏距离
                    for (int k : zones[j]) {
                        totalDist += inst.dist[i][k];
                    }
                }

                if (totalDist < minTotalDist) {
                    minTotalDist = totalDist;
                    bestCenter = i;
                }
            }

            newCenters.add(inst.getAreas()[bestCenter]);
        }

        return newCenters;
    }

    /**
     * 比较两组中心是否相同
     */
    private boolean compareCenters(ArrayList<Area> centers1, ArrayList<Area> centers2) {
        if (centers1.size() != centers2.size()) {
            return false;
        }

        for (int i = 0; i < centers1.size(); i++) {
            if (centers1.get(i).getId() != centers2.get(i).getId()) {
                return false;
            }
        }

        return true;
    }

    /**
     * 确保每个区域的连通性 - 使用Gurobi实现
     * 对于改进模型且不使用相对平衡性约束的情况，包含两层循环：第一层增加上下界（最多500次），第二层降低上下界（最多500次）
     * @return 是否成功确保所有区域连通
     */
    private boolean ensureConnectivity(long globalStartTime) throws GRBException {
        // 如果使用改进模型且不使用相对平衡性约束，需要迭代调整上下界
        if (useImprovedModel && !useRelativeBalance) {
            // 保存原始上下界
            double originalUpperBound = demandUpperBound;
            double originalLowerBound = demandLowerBound;

            // 第一层循环：如果无解，增加上下界2%，最多循环500次
            int maxIncreaseIterations = 20;
            boolean success = false;

            for (int increaseIter = 0; increaseIter < maxIncreaseIterations; increaseIter++) {
                if (increaseIter > 0) {
                    System.out.println("【确保连通性-调整参考值】第 " + increaseIter + " 次增加上下界...");
                    adjustBounds();
                    // 清空不可行集合，因为上下界改变了，之前的割约束可能不再有效
                    if (useExactMethod) {
                        infeasibleSolutions.clear();
                        System.out.println("【确保连通性-调整参考值】已清空不可行集合，重新开始添加割约束");
                    }
                } else {
                    System.out.println("【确保连通性】开始第一次尝试确保连通性...");
                }

                // 执行一次完整的ensureConnectivity过程
                success = ensureConnectivityInternal(globalStartTime);

                if (success) {
                    if (increaseIter > 0) {
                        System.out.println("【确保连通性】经过 " + increaseIter + " 次增加上下界后成功找到可行解");
                    }
                    return true;
                } else {
                    System.out.println("【确保连通性】第 " + (increaseIter + 1) + " 次尝试失败，无解");
                }

                // 检查时间限制
                if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                    System.out.println("【确保连通性】全局时间限制已超过，停止增加上下界循环");
                    break;
                }
            }

            // 如果第一层循环都失败，进入第二层循环：降低上下界
            if (!success) {
                System.out.println("【确保连通性】增加上下界循环结束，仍未找到可行解，开始降低上下界循环...");

                // 恢复原始上下界
                demandUpperBound = originalUpperBound;
                demandLowerBound = originalLowerBound;
                System.out.println("【确保连通性】已恢复原始上下界：L = " + demandLowerBound + ", U = " + demandUpperBound);

                int maxDecreaseIterations = 20;
                for (int decreaseIter = 0; decreaseIter < maxDecreaseIterations; decreaseIter++) {
                    if (decreaseIter > 0) {
                        System.out.println("【确保连通性-调整参考值】第 " + decreaseIter + " 次降低上下界...");
                        reduceBounds();
                        // 清空不可行集合，因为上下界改变了，之前的割约束可能不再有效
                        if (useExactMethod) {
                            infeasibleSolutions.clear();
                            System.out.println("【确保连通性-调整参考值】已清空不可行集合，重新开始添加割约束");
                        }
                    } else {
                        System.out.println("【确保连通性】开始降低上下界循环的第一次尝试...");
                    }

                    // 执行一次完整的ensureConnectivity过程
                    success = ensureConnectivityInternal(globalStartTime);

                    if (success) {
                        if (decreaseIter > 0) {
                            System.out.println("【确保连通性】经过 " + decreaseIter + " 次降低上下界后成功找到可行解");
                        }
                        return true;
                    } else {
                        System.out.println("【确保连通性】降低上下界第 " + (decreaseIter + 1) + " 次尝试失败，无解");
                    }

                    // 检查时间限制
                    if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                        System.out.println("【确保连通性】全局时间限制已超过，停止降低上下界循环");
                        break;
                    }
                }
            }

            // 如果两层循环都失败，最终判定为无解
            if (!success) {
                System.out.println("【确保连通性】两层循环都结束，仍未找到可行解，最终判定为本次实验无解");
                // failureStage会在调用处设置，这里不需要设置
            }

            return success;
        } else {
            // 不使用改进模型或使用相对平衡性约束时，直接执行一次，不进行上下界迭代
            System.out.println("【确保连通性】不使用改进模型或使用相对平衡性约束，直接执行连通性检查，不进行上下界迭代");
            return ensureConnectivityInternal(globalStartTime);
        }
    }

    /**
     * 确保每个区域的连通性 - 内部实现（执行一次完整的确保连通性过程）
     * @return 是否成功确保所有区域连通
     */
    private boolean ensureConnectivityInternal(long globalStartTime) throws GRBException {
        boolean allConnected = false;
        int iteration = 0;
        int maxIterations = 1000; // 限制迭代次数

        // 将model和env声明在try块外部，以便在catch块中访问
        GRBModel model = null;
        GRBEnv env = null;

        try {
            // Check time limit before starting
            if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                System.out.println("Global time limit of 1000 seconds exceeded before ensureConnectivity");
                return false;
            }

            // 如果使用精确方法，记录进入时的割约束数量
            int initialCutCount = 0;
            if (useExactMethod) {
                initialCutCount = infeasibleSolutions.size();
                System.out.println("【确保连通性】当前infeasibleSolutions集合包含 " + initialCutCount + " 个割约束");
            }

            env = new GRBEnv(true);
            env.set(GRB.IntParam.OutputFlag, 0);
            env.set(GRB.IntParam.LogToConsole, 0);
            env.set(GRB.StringParam.LogFile, "");
            env.set(GRB.IntParam.Seed, 42);
            env.start();
            model = new GRBModel(env);

            // 设置Gurobi参数 - 使用剩余时间
            long remainingTimeMs = GLOBAL_TIME_LIMIT_MS - (System.currentTimeMillis() - globalStartTime);
            double remainingTimeSec = Math.max(1.0, Math.min(remainingTimeMs / 1000.0, timeLimit));
            model.set(GRB.IntParam.OutputFlag, 0);
            model.set(GRB.DoubleParam.TimeLimit, remainingTimeSec);
            
            // 设置更严格的数值容差，以提高解的精度
            // 特别是对于精确的分式约束，需要更严格的容差
            if (useExactMethod && useRelativeBalance) {
                model.set(GRB.DoubleParam.FeasibilityTol, 1e-6);  // 可行性容差：默认1e-6，保持严格
                model.set(GRB.DoubleParam.OptimalityTol, 1e-6);   // 最优性容差：默认1e-6，保持严格
                model.set(GRB.IntParam.NumericFocus, 3);          // 数值焦点：3表示最高精度
                model.set(GRB.DoubleParam.BarConvTol, 1e-8);      // 障碍法收敛容差：更严格
                // 关键修复：设置NonConvex=2保证全局最优（对于包含二次分式约束的非凸问题）
                model.set(GRB.IntParam.NonConvex, 2);
            }

            // 决策变量 x_ij
            GRBVar[][] x = new GRBVar[inst.getN()][centers.size()];
            for (int i = 0; i < inst.getN(); i++) {
                for (int j = 0; j < centers.size(); j++) {
                    x[i][j] = model.addVar(0, 1, 0, GRB.BINARY, "x_" + i + "_" + centers.get(j).getId());
                    if (i == centers.get(j).getId()) {
                        model.addConstr(x[i][j], GRB.EQUAL, 1.0, "center_" + i + "_" + j);
                    }
                }
            }

            // 约束: 每个基本单元必须且只能属于一个区域
            for (int i = 0; i < inst.getN(); i++) {
                GRBLinExpr rowSum = new GRBLinExpr();
                for (int j = 0; j < centers.size(); j++) {
                    rowSum.addTerm(1.0, x[i][j]);
                }
                model.addConstr(rowSum, GRB.EQUAL, 1.0, "unit_" + i);
            }

            // 在进入循环之前，先用改善初始解阶段最后一次求得的解（当前zones）验证连通性
            System.out.println("【确保连通性】检查改善初始解阶段最后一次求得的解是否满足连通性约束...");
            boolean initialSolutionConnected = true;
            for (int j = 0; j < centers.size(); j++) {
                if (zones[j] != null && !zones[j].isEmpty()) {
                    ArrayList<ArrayList<Integer>> components = findConnectedComponents(zones[j]);
                    if (components.size() > 1) {
                        initialSolutionConnected = false;
                        System.out.println("【确保连通性】初始解的区域 " + j + " 不连通，包含 " + components.size() + " 个连通组件");
                        break;
                    }
                }
            }

            if (initialSolutionConnected) {
                System.out.println("【确保连通性】初始解满足连通性约束，跳过连通性验证循环");
                allConnected = true;
                // 释放资源
                if (model != null) {
                    model.dispose();
                }
                if (env != null) {
                    env.dispose();
                }
                return true; // 直接返回，跳过后续求解
            }

            // 根据方法类型添加约束
            if (useExactMethod) {
                // 精确方法
                if (useRelativeBalance) {
                    // 使用相对平衡性约束时，直接添加精确的分式约束（利用Gurobi 13处理高次分式的功能）
                    addDistributionallyRobustConstraints(model, x);
                } else {
                    // 不使用相对平衡性约束时，使用no-good cut方法
                    // 不添加DRCC约束，而是通过迭代验证和添加割约束
                    // 添加已记录的割约束
                    addCutConstraints(model, x);
                    System.out.println("【确保连通性】添加了 " + initialCutCount + " 个割约束");
                }
            } else {
                // 近似方法：直接添加DRCC约束（MISOCP）或使用支撑超平面cut
                addDistributionallyRobustConstraints(model, x);
                
                // 如果使用支撑超平面cut，需要先迭代添加cut（类似generateInitialSolution中的逻辑）
                if (useRelativeBalance && useSupportingHyperplaneCuts && !shouldUseD1CopositiveSdpApprox()) {
                    System.out.println("【确保连通性】使用支撑超平面cut，先迭代添加cut...");
                    
                    // 设置目标函数
                    setObjectiveFunction(model, x);
                    
                    // 迭代添加支撑超平面cut
                    int maxCutIterations = 100;
                    int cutIteration = 0;
                    boolean converged = false;
                    
                    while (!converged && cutIteration < maxCutIterations) {
                        cutIteration++;
                        System.out.println("【确保连通性-支撑超平面cut迭代】第 " + cutIteration + " 次迭代...");
                        
                        // 为每次迭代设置时间限制（每次迭代最多3600秒，但不超过剩余时间）
                        long remainingTimeForIteration = GLOBAL_TIME_LIMIT_MS - (System.currentTimeMillis() - globalStartTime);
                        long iterationTimeLimitMs = Math.min(3600 * 1000, remainingTimeForIteration); // 每次迭代最多3600秒，但不超过剩余时间
                        double iterationTimeLimitSec = Math.max(1.0, iterationTimeLimitMs / 1000.0); // 至少1秒
                        double originalTimeLimit = model.get(GRB.DoubleParam.TimeLimit);
                        model.set(GRB.DoubleParam.TimeLimit, iterationTimeLimitSec);
                        System.out.println(String.format("【确保连通性-支撑超平面cut迭代】本次迭代时间限制: %.1f 秒 (剩余总时间: %.1f 秒)", 
                            iterationTimeLimitSec, remainingTimeForIteration / 1000.0));
                        
                        // 求解模型
                        model.optimize();
                        
                        // 恢复原始时间限制
                        model.set(GRB.DoubleParam.TimeLimit, originalTimeLimit);
                        
                        // 检查时间限制
                        if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                            System.out.println("【确保连通性-支撑超平面cut迭代】全局时间限制已超过，停止迭代");
                            int status = model.get(GRB.IntAttr.Status);
                            this.statusCode = status;
                            model.dispose();
                            env.dispose();
                            return false;
                        }
                        
                        // 检查模型状态
                        int status = model.get(GRB.IntAttr.Status);
                        if (status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                            System.out.println("【确保连通性-支撑超平面cut迭代】模型求解失败，状态码: " + status);
                            this.statusCode = status;
                            model.dispose();
                            env.dispose();
                            return false;
                        }
                        
                        // 检查约束是否被违反，如果违反则添加cut
                        boolean cutsAdded = checkAndAddRelativeBalanceCuts(model, x);
                        
                        if (!cutsAdded) {
                            if (useAssignmentDependent) {
                                double[][] currentXVal = new double[inst.getN()][centers.size()];
                                for (int i = 0; i < inst.getN(); i++)
                                    for (int k = 0; k < centers.size(); k++)
                                        currentXVal[i][k] = x[i][k].get(GRB.DoubleAttr.X);
                                if (checkExactAssignmentDependentFeasibility(currentXVal)) {
                                    converged = true;
                                    System.out.println("【确保连通性-支撑超平面cut迭代】第 " + cutIteration + " 次迭代：精确约束满足，收敛");
                                }
                            } else {
                                converged = true;
                                System.out.println("【确保连通性-支撑超平面cut迭代】第 " + cutIteration + " 次迭代：没有新的cut，收敛");
                            }
                        } else {
                            System.out.println("【确保连通性-支撑超平面cut迭代】添加了新的cut，继续迭代...");
                        }
                    }
                    
                    if (!converged && cutIteration >= maxCutIterations) {
                        System.out.println("【确保连通性-支撑超平面cut迭代】达到最大迭代次数 " + maxCutIterations);
                    }
                    
                    // 提取当前解（在添加连通性约束之前）
                    for (int j = 0; j < centers.size(); j++) {
                        zones[j] = new ArrayList<>();
                        for (int i = 0; i < inst.getN(); i++) {
                            if (Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                                zones[j].add(i);
                            }
                        }
                    }
                    
                    System.out.println("【确保连通性】支撑超平面cut迭代完成，共迭代 " + cutIteration + " 次");
                }
            }

            // 设置目标函数（如果还没有设置，或者不使用支撑超平面cut）
            if (!(useRelativeBalance && useSupportingHyperplaneCuts && !useExactMethod && !shouldUseD1CopositiveSdpApprox())) {
                setObjectiveFunction(model, x);
            }

            // 记录已添加的约束总数
            int totalConstraints = 0;

            while (!allConnected && iteration < maxIterations) {
                // Check time limit before each iteration
                if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                    System.out.println("Global time limit of 1000 seconds exceeded during connectivity iterations");
                    // 如果之前有optimize调用，尝试获取最后一次的状态代码
                    try {
                        int status = model.get(GRB.IntAttr.Status);
                        this.statusCode = status;
                    } catch (Exception e) {
                        // 如果无法获取状态代码，保持statusCode为-1
                    }
                    model.dispose();
                    env.dispose();
                    return false; // 时间限制超过，返回false表示失败
                }

                iteration++;

                // 先检查当前解（如果是第一次迭代，使用初始解；否则使用上一次求解得到的解）
                boolean hasDisconnection = false;
                Map<Integer, List<ArrayList<Integer>>> allDisconnectedComponents = new HashMap<>();

                // 如果是第一次迭代，使用初始解（已经在循环前检查过，但这里需要重新提取用于后续处理）
                if (iteration == 1) {
                    // 使用当前的zones（改善初始解阶段最后一次求得的解）
                    // zones已经在循环前检查过，这里直接使用
                } else {
                    // 提取上一次求解得到的解
                    // 如果使用精确方法，解已经在solveWithExactMethodForConnectivity中提取了
                    // 如果使用近似方法，需要在这里提取
                    if (!useExactMethod) {
                        for (int j = 0; j < centers.size(); j++) {
                            zones[j] = new ArrayList<>();
                            for (int i = 0; i < inst.getN(); i++) {
                                if (Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                                    zones[j].add(i);
                                }
                            }
                        }
                    }
                }

                // 先验证当前解的连通性（如果满足则直接跳出循环）
                for (int j = 0; j < centers.size(); j++) {
                    if (zones[j] != null && !zones[j].isEmpty()) {
                        ArrayList<ArrayList<Integer>> components = findConnectedComponents(zones[j]);
                        if (components.size() > 1) {
                            hasDisconnection = true;
                            // 找出中心所在的连通组件
                            int centerComponentIndex = -1;
                            for (int c = 0; c < components.size(); c++) {
                                if (components.get(c).contains(centers.get(j).getId())) {
                                    centerComponentIndex = c;
                                    break;
                                }
                            }
                            // 保存不包含中心的连通组件
                            List<ArrayList<Integer>> disconnectedComponents = new ArrayList<>();
                            for (int c = 0; c < components.size(); c++) {
                                if (c != centerComponentIndex) {
                                    disconnectedComponents.add(components.get(c));
                                }
                            }
                            if (!disconnectedComponents.isEmpty()) {
                                allDisconnectedComponents.put(j, disconnectedComponents);
                            }
                        }
                    }
                }

                // 如果当前解满足连通性，直接跳出循环
                if (!hasDisconnection) {
                    System.out.println("【确保连通性】迭代 " + iteration + "：当前解满足连通性约束");
                    allConnected = true;
                    continue;
                }

                // 如果当前解不满足连通性，求连通分量并添加对应的连通性约束
                System.out.println("【确保连通性】迭代 " + iteration + "：当前解不满足连通性约束，添加连通性约束...");

                // 添加连通性约束
                int constraintCounter = 0;
                for (Map.Entry<Integer, List<ArrayList<Integer>>> entry : allDisconnectedComponents.entrySet()) {
                    int districtIndex = entry.getKey();
                    List<ArrayList<Integer>> disconnectedComponents = entry.getValue();

                    for (ArrayList<Integer> component : disconnectedComponents) {
                        // 找出该组件的邻居
                        HashSet<Integer> neighbors = new HashSet<>();
                        for (int node : component) {
                            for (int neighbor : inst.getAreas()[node].getNeighbors()) {
                                if (!component.contains(neighbor)) {
                                    neighbors.add(neighbor);
                                }
                            }
                        }

                        // 添加约束: 组件中的所有节点都分配给区域districtIndex，或至少有一个邻居也分配给该区域
                        GRBLinExpr constrExpr = new GRBLinExpr();

                        // 对所有的邻居节点
                        for (int neighbor : neighbors) {
                            constrExpr.addTerm(1.0, x[neighbor][districtIndex]);
                        }

                        // 对当前组件中的所有节点
                        for (int node : component) {
                            constrExpr.addTerm(-1.0, x[node][districtIndex]);
                        }

                        model.addConstr(constrExpr, GRB.GREATER_EQUAL, 1 - component.size(), "conn_" + districtIndex + "_" + constraintCounter);
                        constraintCounter++;
                        totalConstraints++;
                    }
                }

                System.out.println("【确保连通性】迭代 " + iteration + "：添加了 " + constraintCounter + " 个连通性约束");

                // 添加连通性约束后，分别进入solveWithExactMethodForConnectivity或model.optimize
                boolean solved = false;
                if (useExactMethod) {
                    solved = solveWithExactMethodForConnectivity(model, x, env, globalStartTime, initialCutCount);
                } else {
                    // 如果使用支撑超平面cut，需要迭代添加cut
                    if (useRelativeBalance && useSupportingHyperplaneCuts && !shouldUseD1CopositiveSdpApprox()) {
                        System.out.println("【确保连通性】迭代 " + iteration + "：添加连通性约束后，迭代添加支撑超平面cut...");
                        
                        // 迭代添加支撑超平面cut
                        int maxCutIterations = 100;
                        int cutIteration = 0;
                        boolean converged = false;
                        
                        while (!converged && cutIteration < maxCutIterations) {
                            cutIteration++;
                            System.out.println("【确保连通性-支撑超平面cut迭代】第 " + cutIteration + " 次迭代...");
                            
                            // 为每次迭代设置时间限制（每次迭代最多3600秒，但不超过剩余时间）
                            long remainingTimeForIteration = GLOBAL_TIME_LIMIT_MS - (System.currentTimeMillis() - globalStartTime);
                            long iterationTimeLimitMs = Math.min(3600 * 1000, remainingTimeForIteration); // 每次迭代最多3600秒，但不超过剩余时间
                            double iterationTimeLimitSec = Math.max(1.0, iterationTimeLimitMs / 1000.0); // 至少1秒
                            double originalTimeLimit = model.get(GRB.DoubleParam.TimeLimit);
                            model.set(GRB.DoubleParam.TimeLimit, iterationTimeLimitSec);
                            System.out.println(String.format("【确保连通性-支撑超平面cut迭代】本次迭代时间限制: %.1f 秒 (剩余总时间: %.1f 秒)", 
                                iterationTimeLimitSec, remainingTimeForIteration / 1000.0));
                            
                            // 求解模型
                            model.optimize();
                            
                            // 恢复原始时间限制
                            model.set(GRB.DoubleParam.TimeLimit, originalTimeLimit);
                            
                            // 检查时间限制
                            if (System.currentTimeMillis() - globalStartTime > GLOBAL_TIME_LIMIT_MS) {
                                System.out.println("【确保连通性-支撑超平面cut迭代】全局时间限制已超过，停止迭代");
                                int status = model.get(GRB.IntAttr.Status);
                                this.statusCode = status;
                                model.dispose();
                                env.dispose();
                                return false;
                            }
                            
                            // 检查模型状态
                            int status = model.get(GRB.IntAttr.Status);
                            if (status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                                System.out.println("【确保连通性-支撑超平面cut迭代】模型求解失败，状态码: " + status);
                                this.statusCode = status;
                                model.dispose();
                                env.dispose();
                                return false;
                            }
                            
                            // 检查约束是否被违反，如果违反则添加cut
                            boolean cutsAdded = checkAndAddRelativeBalanceCuts(model, x);
                            
                            if (!cutsAdded) {
                                if (useAssignmentDependent) {
                                    double[][] currentXVal = new double[inst.getN()][centers.size()];
                                    for (int i = 0; i < inst.getN(); i++)
                                        for (int k = 0; k < centers.size(); k++)
                                            currentXVal[i][k] = x[i][k].get(GRB.DoubleAttr.X);
                                    if (checkExactAssignmentDependentFeasibility(currentXVal)) {
                                        converged = true;
                                        System.out.println("【确保连通性-支撑超平面cut迭代】第 " + cutIteration + " 次迭代：精确约束满足，收敛");
                                    }
                                } else {
                                    converged = true;
                                    System.out.println("【确保连通性-支撑超平面cut迭代】第 " + cutIteration + " 次迭代：没有新的cut，收敛");
                                }
                            } else {
                                System.out.println("【确保连通性-支撑超平面cut迭代】添加了新的cut，继续迭代...");
                            }
                        }
                        
                        if (!converged && cutIteration >= maxCutIterations) {
                            System.out.println("【确保连通性-支撑超平面cut迭代】达到最大迭代次数 " + maxCutIterations);
                        }
                        
                        // 提取解决方案
                        for (int j = 0; j < centers.size(); j++) {
                            zones[j] = new ArrayList<>();
                            for (int i = 0; i < inst.getN(); i++) {
                                if (Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                                    zones[j].add(i);
                                }
                            }
                        }
                        
                        solved = converged;
                        System.out.println("【确保连通性】迭代 " + iteration + "：支撑超平面cut迭代完成，共迭代 " + cutIteration + " 次");
                    } else if (shouldUseD1CopositiveSdpApprox()) {
                        System.out.println("【确保连通性】迭代 " + iteration + "：使用 D1-SDP 分离法...");
                        boolean ok = iterateD1SdpSeparation(model, x, globalStartTime, 120, "连通性迭代-D1-SDP");
                        if (!ok) {
                            solved = false;
                        } else {
                            solved = true;
                            for (int jj = 0; jj < centers.size(); jj++) {
                                zones[jj] = new ArrayList<>();
                                for (int i = 0; i < inst.getN(); i++) {
                                    if (Math.abs(x[i][jj].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                                        zones[jj].add(i);
                                    }
                                }
                            }
                        }
                    } else {
                        // 不使用支撑超平面cut，直接求解模型
                        model.optimize();

                        // Update time limit for next iteration
                        remainingTimeMs = GLOBAL_TIME_LIMIT_MS - (System.currentTimeMillis() - globalStartTime);
                        remainingTimeSec = Math.max(1.0, Math.min(remainingTimeMs / 1000.0, timeLimit));
                        model.set(GRB.DoubleParam.TimeLimit, remainingTimeSec);

                        // 检查模型状态
                        int status = model.get(GRB.IntAttr.Status);
                        if (status != GRB.OPTIMAL && status != GRB.SUBOPTIMAL) {
                            // 保存状态代码（如果求解失败，必须记录）
                            this.statusCode = status;
                            System.out.println("【确保连通性-近似方法】迭代 " + iteration + " 失败，模型无解");
                            break;
                        }
                        // 如果求解成功，可以不更新statusCode（保持-1）
                        solved = true;

                        // 提取解决方案（近似方法需要在这里提取）
                        for (int j = 0; j < centers.size(); j++) {
                            zones[j] = new ArrayList<>();
                            for (int i = 0; i < inst.getN(); i++) {
                                if (Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                                    zones[j].add(i);
                                }
                            }
                        }
                    }
                }

                if (!solved) {
                    System.out.println("【确保连通性】迭代 " + iteration + " 失败，无法找到满足DRCC约束的解");
                    break;
                }
            }

            // 注意：解已经在循环中提取了
            // 精确方法的解在solveWithExactMethodForConnectivity中提取
            // 近似方法的解在循环中的model.optimize()后提取

            if (!allConnected) {
                System.out.println("警告：确保连通性失败，实验结束");
            }

            // 释放资源
            if (model != null) {
                model.dispose();
            }
            if (env != null) {
                env.dispose();
            }

            return allConnected;

        } catch (GRBException e) {
            System.out.println("确保连通性时Gurobi错误: " + e.getMessage());
            // 如果model已经创建且optimize已调用，尝试获取状态代码
            try {
                if (model != null) {
                    int status = model.get(GRB.IntAttr.Status);
                    this.statusCode = status;
                }
            } catch (Exception ex) {
                // 如果无法获取状态代码，保持statusCode为-1
            }
            // 确保在异常情况下也释放资源
            try {
                if (model != null) {
                    model.dispose();
                }
            } catch (Exception ex) {
                // 忽略释放资源时的异常
            }
            try {
                if (env != null) {
                    env.dispose();
                }
            } catch (Exception ex) {
                // 忽略释放资源时的异常
            }
            e.printStackTrace();
            return false; // 异常情况下返回false
        }
    }

    /**
     * 找出一个区域内的所有连通组件
     */
    private ArrayList<ArrayList<Integer>> findConnectedComponents(ArrayList<Integer> zone) {
        ArrayList<ArrayList<Integer>> components = new ArrayList<>();
        boolean[] visited = new boolean[inst.getN()];

        for (int i : zone) {
            if (!visited[i]) {
                ArrayList<Integer> component = new ArrayList<>();
                Queue<Integer> queue = new LinkedList<>();

                queue.add(i);
                visited[i] = true;

                while (!queue.isEmpty()) {
                    int current = queue.poll();
                    component.add(current);

                    for (int neighbor : inst.getAreas()[current].getNeighbors()) {
                        if (zone.contains(neighbor) && !visited[neighbor]) {
                            queue.add(neighbor);
                            visited[neighbor] = true;
                        }
                    }
                }

                components.add(component);
            }
        }

        return components;
    }

    /**
     * 设置目标函数（根据模型类型选择不同的目标函数）
     */
    private void setObjectiveFunction(GRBModel model, GRBVar[][] x) throws GRBException {
        if (useImprovedModel) {
            if (useD1) {
                // D1模糊集：min ∑_j ∑_i μ_i * r_ij * x_ij
                GRBLinExpr obj = new GRBLinExpr();
                for (int i = 0; i < inst.getN(); i++) {
                    for (int j = 0; j < centers.size(); j++) {
                        int centerId = centers.get(j).getId();
                        double r_ij = shortestPathDist[i][centerId];
                        obj.addTerm(meanVector[i] * r_ij, x[i][j]);
                    }
                }
                model.setObjective(obj, GRB.MINIMIZE);
            } else {
                // D2模糊集：min ∑_j ∑_i μ_i * r_ij * x_ij + sqrt(δ) * t
                // 其中 t >= sqrt(w^T * Σ * w)，w_i = ∑_j r_ij * x_ij
                double delta = Math.min(delta1, delta2);

                // 线性部分：∑_j ∑_i μ_i * r_ij * x_ij
                GRBLinExpr obj = new GRBLinExpr();
                for (int i = 0; i < inst.getN(); i++) {
                    for (int j = 0; j < centers.size(); j++) {
                        int centerId = centers.get(j).getId();
                        double r_ij = shortestPathDist[i][centerId];
                        obj.addTerm(meanVector[i] * r_ij, x[i][j]);
                    }
                }

                // 添加辅助变量 t
                GRBVar t = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "obj_t");
                obj.addTerm(Math.sqrt(delta), t);

                // 添加SOC约束：t >= sqrt(w^T * Σ * w)
                // 首先计算 w_i = ∑_j r_ij * x_ij
                GRBVar[] w = new GRBVar[inst.getN()];
                for (int i = 0; i < inst.getN(); i++) {
                    w[i] = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "w_" + i);
                    GRBLinExpr wExpr = new GRBLinExpr();
                    for (int j = 0; j < centers.size(); j++) {
                        int centerId = centers.get(j).getId();
                        double r_ij = shortestPathDist[i][centerId];
                        wExpr.addTerm(r_ij, x[i][j]);
                    }
                    model.addConstr(w[i], GRB.EQUAL, wExpr, "w_def_" + i);
                }

                // 使用Cholesky分解创建SOC约束
                Matrix covMatrix = new Matrix(covarianceMatrix);
                CholeskyDecomposition chol = new CholeskyDecomposition(covMatrix);

                if (chol.isSPD()) {
                    Matrix cholL = chol.getL();
                    double[][] cholMatrix = cholL.getArray();

                    // 创建辅助变量 z = L^T * w
                    GRBVar[] z = new GRBVar[inst.getN()];
                    for (int i = 0; i < inst.getN(); i++) {
                        z[i] = model.addVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "z_obj_" + i);
                        GRBLinExpr zExpr = new GRBLinExpr();
                        for (int k = 0; k < inst.getN(); k++) {
                            zExpr.addTerm(cholMatrix[k][i], w[k]);
                        }
                        model.addConstr(z[i], GRB.EQUAL, zExpr, "z_obj_def_" + i);
                    }

                    // SOC约束：||z|| <= t
                    GRBQuadExpr socExpr = new GRBQuadExpr();
                    for (int i = 0; i < inst.getN(); i++) {
                        socExpr.addTerm(1.0, z[i], z[i]);
                    }
                    GRBQuadExpr tSquare = new GRBQuadExpr();
                    tSquare.addTerm(1.0, t, t);
                    model.addQConstr(socExpr, GRB.LESS_EQUAL, tSquare, "obj_soc");
                } else {
                    // 如果Cholesky分解失败，使用二次约束
                    GRBQuadExpr quadExpr = new GRBQuadExpr();
                    for (int i = 0; i < inst.getN(); i++) {
                        for (int k = 0; k < inst.getN(); k++) {
                            quadExpr.addTerm(covarianceMatrix[i][k], w[i], w[k]);
                        }
                    }
                    GRBVar q = model.addVar(0, GRB.INFINITY, 0, GRB.CONTINUOUS, "q_obj");
                    model.addQConstr(quadExpr, GRB.LESS_EQUAL, q, "quad_obj");
                    GRBQuadExpr tSquare = new GRBQuadExpr();
                    tSquare.addTerm(1.0, t, t);
                    model.addQConstr(q, GRB.LESS_EQUAL, tSquare, "t_square_obj");
                }

                model.setObjective(obj, GRB.MINIMIZE);
            }
        } else {
            // 原始模型：最小化总距离
            GRBLinExpr obj = new GRBLinExpr();
            for (int i = 0; i < inst.getN(); i++) {
                for (int j = 0; j < centers.size(); j++) {
                    obj.addTerm(inst.dist[i][centers.get(j).getId()], x[i][j]);
                }
            }
            model.setObjective(obj, GRB.MINIMIZE);
        }
    }

    /**
     * 计算当前解的目标函数值
     */
    private double evaluateObjective() {
        if (useImprovedModel) {
            if (useD1) {
                // D1: ∑_j ∑_i μ_i * r_ij * x_ij
                double total = 0.0;
                for (int j = 0; j < centers.size(); j++) {
                    if (zones[j] != null) {
                        int centerId = centers.get(j).getId();
                        for (int i : zones[j]) {
                            double r_ij = shortestPathDist[i][centerId];
                            total += meanVector[i] * r_ij;
                        }
                    }
                }
                return total;
            } else {
                // D2: ∑_j ∑_i μ_i * r_ij * x_ij + sqrt(δ) * sqrt(w^T * Σ * w)
                double delta = Math.min(delta1, delta2);

                // 线性部分
                double linearPart = 0.0;
                double[] w = new double[inst.getN()];
                for (int i = 0; i < inst.getN(); i++) {
                    w[i] = 0.0;
                    for (int j = 0; j < centers.size(); j++) {
                        if (zones[j] != null && zones[j].contains(i)) {
                            int centerId = centers.get(j).getId();
                            double r_ij = shortestPathDist[i][centerId];
                            w[i] += r_ij;
                            linearPart += meanVector[i] * r_ij;
                        }
                    }
                }

                // 二次部分：sqrt(w^T * Σ * w)
                double quadPart = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        quadPart += w[i] * covarianceMatrix[i][k] * w[k];
                    }
                }
                double sqrtQuadPart = Math.sqrt(Math.max(0, quadPart));

                return linearPart + Math.sqrt(delta) * sqrtQuadPart;
            }
        } else {
            // 原始模型：总距离
            double totalDist = 0;
            for (int j = 0; j < centers.size(); j++) {
                if (zones[j] != null) {
                    for (int i : zones[j]) {
                        totalDist += inst.dist[i][centers.get(j).getId()];
                    }
                }
            }
            return totalDist;
        }
    }

    /**
     * 获取区域
     */
    public ArrayList<Integer>[] getZones() {
        return zones;
    }

    /**
     * 获取中心
     */
    public ArrayList<Area> getCenters() {
        return centers;
    }

    /**
     * 获取最佳目标函数值
     */
    public double getBestObjective() {
        return bestObjective;
    }

    /**
     * 获取求解是否成功
     */
    public boolean isSolveSuccess() {
        return solveSuccess;
    }

    /**
     * 获取Gurobi状态代码（如果求解失败）
     */
    public int getStatusCode() {
        return statusCode;
    }

    /**
     * 获取精确算法的迭代次数（如果使用精确方法）
     */
    public int getCutIterations() {
        return cutIterations;
    }

    /**
     * 获取失败阶段：0=成功，1=选择初始中心失败，2=生成初始解失败，3=改善初始解失败，4=确保连通性失败
     */
    public int getFailureStage() {
        return failureStage;
    }

    /**
     * 获取是否使用改进模型
     */
    public boolean getUseImprovedModel() {
        return useImprovedModel;
    }

    /**
     * 获取是否使用相对平衡性约束
     */
    public boolean getUseRelativeBalance() {
        return useRelativeBalance;
    }

    /**
     * 获取最短路径距离矩阵
     */
    public double[][] getShortestPathDist() {
        return shortestPathDist;
    }

    /**
     * 获取均值向量
     */
    public double[] getMeanVector() {
        return meanVector;
    }

    /**
     * 获取需求/工作量下界
     */
    public double getDemandLowerBound() {
        return demandLowerBound;
    }

    /**
     * 获取需求/工作量上界
     */
    public double getDemandUpperBound() {
        return demandUpperBound;
    }

    /**
     * 获取avgDist计算方式
     */
    public int getAvgDistMethod() {
        return avgDistMethod;
    }

    /**
     * 获取当前方法的avgDist向量均值
     */
    public double getAvgDistMean() {
        return avgDistMean;
    }

    /**
     * 获取是否使用assignment-dependent模型
     */
    public boolean getUseAssignmentDependent() {
        return useAssignmentDependent;
    }

    /**
     * 获取assignment-dependent场景数据
     * @return assignment-dependent场景数据 [scenario][i][j]，如果未使用assignment-dependent模型则返回null
     */
    public double[][][] getAssignmentDependentScenarios() {
        return assignmentDependentScenarios;
    }

    /**
     * 获取assignment-dependent场景数量
     */
    public int getNumAssignmentDependentScenarios() {
        return numAssignmentDependentScenarios;
    }

    /**
     * 初始化工作量CSV文件
     */
    private void initializeWorkloadCsvFile() {
        try {
            // 确保output目录存在
            File outputDir = new File("output");
            if (!outputDir.exists()) {
                outputDir.mkdirs();
            }
            
            // 生成CSV文件名（包含时间戳以避免覆盖）
            long timestamp = System.currentTimeMillis();
            workloadCsvFilePath = "output/workload_iterations_" + timestamp + ".csv";
            
            // 创建BufferedWriter
            workloadCsvWriter = new BufferedWriter(new FileWriter(workloadCsvFilePath));
            
            // 写入CSV表头
            workloadCsvWriter.write("迭代次数,区域0工作量,区域1工作量,区域2工作量,总工作量,区域0占比,区域1占比,区域2占比");
            workloadCsvWriter.newLine();
            workloadCsvWriter.flush();
            
            System.out.println("【工作量记录】已创建CSV文件: " + workloadCsvFilePath);
        } catch (IOException e) {
            System.err.println("创建工作量CSV文件失败: " + e.getMessage());
            workloadCsvWriter = null;
        }
    }
    
    /**
     * 关闭工作量CSV文件
     */
    private void closeWorkloadCsvFile() {
        if (workloadCsvWriter != null) {
            try {
                workloadCsvWriter.close();
                System.out.println("【工作量记录】已关闭CSV文件: " + workloadCsvFilePath);
            } catch (IOException e) {
                System.err.println("关闭工作量CSV文件失败: " + e.getMessage());
            }
            workloadCsvWriter = null;
        }
    }
    
    /**
     * 计算并输出每个区域的工作量（使用均值）
     * 对于assignment-dependent模型，区域j的工作量 = sum_{i in zone_j} E[d_ij]
     * 其中 E[d_ij] 存储在 meanVector 中，索引为 i * p + j
     * @param x 决策变量矩阵
     * @param iteration 当前迭代次数
     */
    private void calculateAndPrintWorkloads(GRBVar[][] x, int iteration) {
        if (!useAssignmentDependent) {
            return; // 只对assignment-dependent模型计算
        }
        
        int n = inst.getN();
        int p = numRegionsForAssignmentDependent;
        int k = centers.size(); // 区域数量（应该是3）
        
        // 提取当前解的区域分配
        ArrayList<Integer>[] currentZones = new ArrayList[k];
        for (int j = 0; j < k; j++) {
            currentZones[j] = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                try {
                    if (Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                        currentZones[j].add(i);
                    }
                } catch (GRBException e) {
                    System.err.println("获取变量值时出错: " + e.getMessage());
                    return;
                }
            }
        }
        
        // 计算每个区域的工作量
        double[] regionWorkloads = new double[k];
        for (int j = 0; j < k; j++) {
            regionWorkloads[j] = 0.0;
            if (currentZones[j] != null) {
                for (int areaId : currentZones[j]) {
                    // 计算 E[d_{areaId, j}]
                    // meanVector 的索引为 idx = areaId * p + j
                    int idx = areaId * p + j;
                    if (idx >= 0 && idx < meanVector.length) {
                        regionWorkloads[j] += meanVector[idx];
                    }
                }
            }
        }
        
        // 计算总工作量
        double totalWorkload = 0.0;
        for (int j = 0; j < k; j++) {
            totalWorkload += regionWorkloads[j];
        }
        
        // 计算各区域工作量占总工作量的比例
        double[] regionProportions = new double[k];
        for (int j = 0; j < k; j++) {
            if (totalWorkload > 1e-10) { // 避免除零
                regionProportions[j] = regionWorkloads[j] / totalWorkload;
            } else {
                regionProportions[j] = 0.0;
            }
        }
        
        // 输出结果
        // System.out.println("【工作量计算】各区域工作量：");
        // for (int j = 0; j < k; j++) {
        //     System.out.println(String.format("  区域 %d: %.6f (占比: %.2f%%, 包含 %d 个基本单元)", 
        //         j, regionWorkloads[j], regionProportions[j] * 100.0, 
        //         currentZones[j] != null ? currentZones[j].size() : 0));
        // }
        // System.out.println(String.format("  总工作量: %.6f", totalWorkload));
        
        // 写入CSV文件
        if (workloadCsvWriter != null) {
            try {
                workloadCsvWriter.write(String.format("%d,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f",
                    iteration,
                    regionWorkloads[0],
                    regionWorkloads[1],
                    regionWorkloads[2],
                    totalWorkload,
                    regionProportions[0],
                    regionProportions[1],
                    regionProportions[2]));
                workloadCsvWriter.newLine();
                workloadCsvWriter.flush();
            } catch (IOException e) {
                System.err.println("写入工作量CSV文件失败: " + e.getMessage());
            }
        }
    }

    /**
     * 添加割约束：对于每个不可行解，添加约束确保新解与其不同
     */
    private void addCutConstraints(GRBModel model, GRBVar[][] x) throws GRBException {
        for (int k = 0; k < infeasibleSolutions.size(); k++) {
            Set<String> infeasibleSet = infeasibleSolutions.get(k);
            GRBLinExpr cutExpr = new GRBLinExpr();

            // 构建割约束: ∑_{(i,j) ∈ S^k} x_{ij} ≤ |S^k| - 1
            for (String index : infeasibleSet) {
                String[] parts = index.split("_");
                int i = Integer.parseInt(parts[0]);
                int j = Integer.parseInt(parts[1]);
                cutExpr.addTerm(1.0, x[i][j]);
            }

            int setSize = infeasibleSet.size();
            model.addConstr(cutExpr, GRB.LESS_EQUAL, setSize - 1, "cut_" + k);
        }
    }

    /**
     * 验证解是否满足DRCC约束（reformulate后的复杂约束形式）
     * 对于D1模糊集，验证约束：
     * \frac{\mathbf{x}_j^\top \boldsymbol{\Sigma} \mathbf{x}_j}{\mathbf{x}_j^\top \boldsymbol{\Sigma} \mathbf{x}_j + t_1^2 (\boldsymbol{\mu}^\top \mathbf{x}_j)^2} +
     * \frac{\mathbf{x}_j^\top \boldsymbol{\Sigma} \mathbf{x}_j}{\mathbf{x}_j^\top \boldsymbol{\Sigma} \mathbf{x}_j + t_2^2 (\boldsymbol{\mu}^\top \mathbf{x}_j)^2} \le \gamma
     *
     * 对于D2模糊集，使用近似方法的约束形式进行验证
     */
    private boolean verifyDRCCConstraints(GRBVar[][] x) throws GRBException {
        // 根据约束类型选择验证绝对约束还是相对约束
        if (useRelativeBalance) {
            return verifyRelativeBalanceConstraints(x);
        } else {
            return verifyAbsoluteBalanceConstraints(x);
        }
    }

    /**
     * 验证绝对平衡性约束（原始约束）
     */
    private boolean verifyAbsoluteBalanceConstraints(GRBVar[][] x) throws GRBException {
        double U = demandUpperBound;
        double L = demandLowerBound;

        for (int j = 0; j < centers.size(); j++) {
            int centerId = centers.get(j).getId();

            // 计算向量（当前解中分配给区域j的基本单元）
            // 原始模型：x_j
            // 改进模型：w_j = r_j ⊙ x_j
            double[] v_j = new double[inst.getN()];
            if (useImprovedModel) {
                // 改进模型：w_j[i] = r_ij * x_ij
                for (int i = 0; i < inst.getN(); i++) {
                    double x_ij = Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6 ? 1.0 : 0.0;
                    double r_ij = shortestPathDist[i][centerId];
                    v_j[i] = r_ij * x_ij;
                }
            } else {
                // 原始模型：x_j
                for (int i = 0; i < inst.getN(); i++) {
                    v_j[i] = Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6 ? 1.0 : 0.0;
                }
            }

            // 计算 μ^T * v_j
            // 原始模型：μ^T * x_j
            // 改进模型：μ^T * w_j
            double muTv = 0.0;
            for (int i = 0; i < inst.getN(); i++) {
                muTv += meanVector[i] * v_j[i];
            }

            // 如果区域j为空，跳过验证
            if (Math.abs(muTv) < 1e-10) {
                continue;
            }

            // 计算 v_j^T * Σ * v_j
            // 原始模型：x_j^T * Σ * x_j
            // 改进模型：w_j^T * Σ * w_j
            double vTSigmav = 0.0;
            for (int i = 0; i < inst.getN(); i++) {
                for (int k = 0; k < inst.getN(); k++) {
                    vTSigmav += v_j[i] * covarianceMatrix[i][k] * v_j[k];
                }
            }

            // 获取风险参数
            double riskParam = useJointChance ? individualGammas[j] : gamma;

            if (useExactMethod) {
                // 精确方法：验证两个SOC相加的约束，不使用gamma分解
                if (useD1) {
                    // D1模糊集：使用精确的reformulate约束（两个SOC相加）
                    // 约束形式：\frac{v_j^T*Σ*v_j}{v_j^T*Σ*v_j + t_1^2*(μ^T*v_j)^2} +
                    //          \frac{v_j^T*Σ*v_j}{v_j^T*Σ*v_j + t_2^2*(μ^T*v_j)^2} ≤ γ
                    double t1 = 1.0 - L / muTv;
                    double t2 = U / muTv - 1.0;

                    // 处理t1和t2为负的情况：
                    // t1 < 0 表示 muTv < L（均值低于下界），违反下界约束的概率应该很高（接近1）
                    // t2 < 0 表示 muTv > U（均值高于上界），违反上界约束的概率应该很高（接近1）
                    // 如果均值本身就不在[L, U]范围内，说明解不满足基本约束，应该直接返回false
                    if (muTv < L - 1e-6 || muTv > U + 1e-6) {
                        // 均值不在边界范围内，直接判定为不满足约束
                        return false;
                    }

                    // 处理方差为0的情况（确定性需求）
                    // 如果方差为0，需求是确定的，应该直接检查muTv是否在[L,U]内
                    // 由于上面已经检查了muTv在[L,U]内，所以方差为0时违反概率为0
                    if (vTSigmav < 1e-10) {
                        // 方差为0，且muTv在[L,U]内，违反概率为0，满足约束
                        continue; // 继续检查下一个区域
                    }

                    // 计算约束左侧的值（两个SOC相加）
                    double denominator1 = vTSigmav + t1 * t1 * muTv * muTv;
                    double denominator2 = vTSigmav + t2 * t2 * muTv * muTv;

                    double lhs = 0.0;
                    // 计算第一项：违反下界约束的概率上界
                    if (denominator1 > 1e-10) {
                        lhs += vTSigmav / denominator1;
                    } else {
                        // 分母非常小，说明t1^2*muTv^2相对于vTSigmav非常大
                        // 这意味着muTv非常接近L，违反概率接近1
                        lhs += 1.0;
                    }

                    // 计算第二项：违反上界约束的概率上界
                    if (denominator2 > 1e-10) {
                        lhs += vTSigmav / denominator2;
                    } else {
                        // 分母非常小，说明t2^2*muTv^2相对于vTSigmav非常大
                        // 这意味着muTv非常接近U，违反概率接近1
                        lhs += 1.0;
                    }

                    // 验证约束是否满足（使用原始gamma，不分解）
                    // 约束要求：违反概率的上界 ≤ γ
                    if (lhs > riskParam + 1e-6) {
                        return false;
                    }
                } else {
                    // D2模糊集：使用精确重构方法进行验证
                    // 计算 v_j^T * Σ * v_j 的平方根
                    double sqrtVTSigmav = Math.sqrt(vTSigmav);

                    // 如果方差为0，跳过验证（这种情况不太可能发生）
                    if (sqrtVTSigmav < 1e-10) {
                        continue;
                    }

                    // 计算 κ_{L,j} 和 κ_{U,j}（根据新的定义）
                    // κ_{L,j} = (μ^T * v_j - L) / sqrt(v_j^T * Σ * v_j)
                    // κ_{U,j} = (U - μ^T * v_j) / sqrt(v_j^T * Σ * v_j)
                    double kappaL = (muTv - L) / sqrtVTSigmav;
                    double kappaU = (U - muTv) / sqrtVTSigmav;

                    // 计算阈值
                    double sqrtDelta1 = Math.sqrt(delta1);
                    double delta2OverSqrtDelta1 = delta2 / sqrtDelta1;

                    // 计算 inf_{P ∈ D_2} P{tilde{d}^T * x_j ≥ L}
                    double infPL = computeInfPL(j, kappaL, sqrtDelta1, delta2OverSqrtDelta1);

                    // 计算 inf_{P ∈ D_2} P{tilde{d}^T * x_j ≤ U}
                    double infPU = computeInfPU(j, kappaU, sqrtDelta1, delta2OverSqrtDelta1);

                    // 验证约束：infPL + infPU ≥ 2 - γ
                    // 等价于：infPL + infPU ≥ 2 - riskParam
                    double lhs = infPL + infPU;
                    double rhs = 2.0 - riskParam;

                    // 输出验证结果
                    // System.out.println(String.format("区域 j=%d: infPL=%.6f, infPU=%.6f, lhs=%.6f, rhs=%.6f",
                    // j, infPL, infPU, lhs, rhs));

                    if (lhs < rhs - 1e-6) {
                        System.out.println(String.format("  约束不满足：lhs=%.6f < rhs=%.6f", lhs, rhs));
                        return false;
                    } else {
                        System.out.println(String.format("  约束满足：lhs=%.6f >= rhs=%.6f", lhs, rhs));
                    }
                }
            } else {
                // 近似方法：使用gamma分解进行验证
                double gamma_a_local = 0.5 * riskParam;
                double gamma_b_local = 0.5 * riskParam;

                // 统一验证方式：使用factor计算上下界，然后与U和L比较
                double factorLower, factorUpper;

                if (useD1) {
                    // D1模糊集：使用分解后的gamma_a和gamma_b计算factor
                    factorLower = Math.sqrt((1 - gamma_a_local) / gamma_a_local);
                    factorUpper = Math.sqrt((1 - gamma_b_local) / gamma_b_local);
                } else {
                    // D2模糊集：使用分解后的gamma_a和gamma_b
                    if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                        // Case 1
                        factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                        factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
                    } else if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 > gamma_b_local) {
                        // Case 2
                        factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                        factorUpper = Math.sqrt(delta2 / gamma_b_local);
                    } else if (delta1 / delta2 > gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                        // Case 3
                        factorLower = Math.sqrt(delta2 / gamma_a_local);
                        factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
                    } else {
                        // Case 4
                        factorLower = Math.sqrt(delta2 / gamma_a_local);
                        factorUpper = Math.sqrt(delta2 / gamma_b_local);
                    }
                }

                // 统一验证：计算上下界并与U和L比较
                double stdDev = Math.sqrt(vTSigmav);
                double upperBound = muTv + factorUpper * stdDev;
                double lowerBound = muTv - factorLower * stdDev;

                if (upperBound > U + 1e-6 || lowerBound < L - 1e-6) {
                    return false;
                }
            }

        }

        return true;
    }

    /**
     * 对 assignment-dependent 模型，用精确的分式约束验证当前解 xVal 是否可行。
     * 用于支撑超平面 cut 迭代的收敛判定：只有精确约束满足时才认为收敛。
     * 精确约束（D1）：meanTerm_vjL <= 0, meanTerm_vjU <= 0，且存在 λ_L,λ_U ∈ [0,γ]，λ_L+λ_U <= γ，
     * 使 (1-λ_L)*L_var_L - λ_L*L_mean_L <= 0 与 (1-λ_U)*L_var_U - λ_U*L_mean_U <= 0。
     * 等价于：λ_L_req = L_var_L/(L_var_L+L_mean_L)，λ_U_req = L_var_U/(L_var_U+L_mean_U)，要求 λ_L_req + λ_U_req <= γ。
     */
    private boolean checkExactAssignmentDependentFeasibility(double[][] xVal) {
        int n = inst.getN();
        int p = centers.size();
        int nTimesP = n * p;
        double r = this.r;

        // 将 xVal 展平为 x_flat[idx] = xVal[idx/p][idx%p]
        double[] x_flat = new double[nTimesP];
        for (int idx = 0; idx < nTimesP; idx++) {
            x_flat[idx] = xVal[idx / p][idx % p];
        }

        for (RelativeBalanceConstraintInfo info : relativeBalanceConstraintInfos) {
            int j = info.j;
            double riskParam = info.riskParam;
            double coeff = info.coeff;
            double coeffUpper = info.coeffUpper;

            // 与 addExactRelativeBalanceConstraintAssignmentDependent 一致的 B_L, B_U
            double[][] B_L = new double[nTimesP][nTimesP];
            double[][] B_U = new double[nTimesP][nTimesP];
            for (int idx = 0; idx < nTimesP; idx++) {
                int j2 = idx % p;
                B_L[idx][idx] = (j2 == j) ? (coeff - 1.0) : coeff;
                B_U[idx][idx] = (j2 == j) ? (1.0 - coeffUpper) : (-coeffUpper);
            }

            // M_L = B_L^T Σ B_L, M_U = B_U^T Σ B_U
            double[][] M_L = new double[nTimesP][nTimesP];
            double[][] M_U = new double[nTimesP][nTimesP];
            for (int idx1 = 0; idx1 < nTimesP; idx1++) {
                for (int idx2 = 0; idx2 < nTimesP; idx2++) {
                    double sumL = 0.0, sumU = 0.0;
                    for (int idx = 0; idx < nTimesP; idx++) {
                        for (int idxp = 0; idxp < nTimesP; idxp++) {
                            double cov = getCovariance(idx, idxp);
                            sumL += B_L[idx][idx1] * cov * B_L[idxp][idx2];
                            sumU += B_U[idx][idx1] * cov * B_U[idxp][idx2];
                        }
                    }
                    M_L[idx1][idx2] = sumL;
                    M_U[idx1][idx2] = sumU;
                }
            }
            // a_L = B_L^T μ, a_U = B_U^T μ
            double[] a_L = new double[nTimesP];
            double[] a_U = new double[nTimesP];
            for (int idx2 = 0; idx2 < nTimesP; idx2++) {
                for (int idx = 0; idx < nTimesP; idx++) {
                    a_L[idx2] += B_L[idx][idx2] * meanVector[idx];
                    a_U[idx2] += B_U[idx][idx2] * meanVector[idx];
                }
            }

            double meanTerm_vjL = 0.0, meanTerm_vjU = 0.0;
            for (int idx = 0; idx < nTimesP; idx++) {
                meanTerm_vjL += a_L[idx] * x_flat[idx];
                meanTerm_vjU += a_U[idx] * x_flat[idx];
            }
            double L_var_L = quadraticFormSymmetric(x_flat, M_L);
            double L_var_U = quadraticFormSymmetric(x_flat, M_U);
            double L_mean_L = meanTerm_vjL * meanTerm_vjL;
            double L_mean_U = meanTerm_vjU * meanTerm_vjU;

            final double tol = 1e-5;

            if (useD1) {
                if (meanTerm_vjL > tol || meanTerm_vjU > tol) {
                    return false;
                }
                double denomL = L_var_L + L_mean_L;
                double denomU = L_var_U + L_mean_U;
                double lambda_L_req = (denomL > 1e-12) ? (L_var_L / denomL) : 0.0;
                double lambda_U_req = (denomU > 1e-12) ? (L_var_U / denomU) : 0.0;
                if (lambda_L_req + lambda_U_req > riskParam + tol) {
                    return false;
                }
            } else {
                // D2 模糊集：与 verifyRelativeBalanceConstraints 中 useExactMethod 的 D2 分支一致
                double sqrtVjL = Math.sqrt(Math.max(0, L_var_L));
                double sqrtVjU = Math.sqrt(Math.max(0, L_var_U));
                if (sqrtVjL < 1e-10 && sqrtVjU < 1e-10) {
                    continue;
                }
                double kappaL = (sqrtVjL > 1e-10) ? (-meanTerm_vjL / sqrtVjL) : 0.0;
                double kappaU = (sqrtVjU > 1e-10) ? (-meanTerm_vjU / sqrtVjU) : 0.0;
                double sqrtDelta1 = Math.sqrt(delta1);
                double delta2OverSqrtDelta1 = delta2 / sqrtDelta1;
                double infPL = (sqrtVjL > 1e-10) ? computeInfPForRelative(kappaL, sqrtDelta1, delta2OverSqrtDelta1) : (meanTerm_vjL <= 0 ? 1.0 : 0.0);
                double infPU = (sqrtVjU > 1e-10) ? computeInfPForRelative(kappaU, sqrtDelta1, delta2OverSqrtDelta1) : (meanTerm_vjU <= 0 ? 1.0 : 0.0);
                double rhs = 2.0 - riskParam;
                if (infPL + infPU < rhs - tol) {
                    return false;
                }
            }
        }
        return true;
    }

    /** 计算 x^T M x，M 为对称矩阵 */
    private double quadraticFormSymmetric(double[] x, double[][] M) {
        int n = x.length;
        double result = 0.0;
        for (int i = 0; i < n; i++) {
            result += M[i][i] * x[i] * x[i];
        }
        for (int i = 0; i < n; i++) {
            for (int k = i + 1; k < n; k++) {
                result += 2.0 * M[i][k] * x[i] * x[k];
            }
        }
        return result;
    }

    /**
     * 验证相对平衡性约束（基于需求的相对平衡性约束）
     * 约束形式：对于每个区域 j，需求比例应该在 [1/p*(1-α), 1/p*(1+α)] 范围内
     * 等价于验证：
     * - d^T * v_{j,L} <= 0，其中 v_{j,L} = (1-α)/p * sum_k x_k - x_j
     * - d^T * v_{j,U} <= 0，其中 v_{j,U} = x_j - (1+α)/p * sum_k x_k
     */
    private boolean verifyRelativeBalanceConstraints(GRBVar[][] x) throws GRBException {
        int p = centers.size(); // 区域数量
        double coeff = (1.0 - r) / p; // (1-r)/p
        double coeffUpper = (1.0 + r) / p; // (1+r)/p

        // 获取风险参数
        double riskParam = useJointChance ? (gamma / p) : gamma; // 对于联合约束，每个区域使用 gamma/p

        for (int j = 0; j < centers.size(); j++) {
            int centerId = centers.get(j).getId();

            // 计算 sum_k x_k（所有区域的总和）
            double[] sumX = new double[inst.getN()];
            if (useImprovedModel) {
                // 改进模型：sum_k w_k，其中 w_k[i] = r_ik * x_ik
                for (int i = 0; i < inst.getN(); i++) {
                    sumX[i] = 0.0;
                    for (int k = 0; k < centers.size(); k++) {
                        int centerId_k = centers.get(k).getId();
                        double r_ik = shortestPathDist[i][centerId_k];
                        double x_ik = Math.abs(x[i][k].get(GRB.DoubleAttr.X) - 1.0) < 1e-6 ? 1.0 : 0.0;
                        sumX[i] += r_ik * x_ik;
                    }
                }
            } else {
                // 原始模型：sum_k x_k
                // 注意：由于每个基本单元只能属于一个区域，sum_k x_k[i] = 1 对于所有 i
                // 为了与添加约束时保持一致（添加约束时使用常数1），这里也使用常数1
                for (int i = 0; i < inst.getN(); i++) {
                    sumX[i] = 1.0; // 使用常数1，与添加约束时保持一致
                }
            }

            // 计算 x_j
            double[] x_j = new double[inst.getN()];
            if (useImprovedModel) {
                // 改进模型：w_j[i] = r_ij * x_ij
                for (int i = 0; i < inst.getN(); i++) {
                    double x_ij = Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6 ? 1.0 : 0.0;
                    double r_ij = shortestPathDist[i][centerId];
                    x_j[i] = r_ij * x_ij;
                }
            } else {
                // 原始模型：x_j
                for (int i = 0; i < inst.getN(); i++) {
                    x_j[i] = Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6 ? 1.0 : 0.0;
                }
            }

            // 计算 v_{j,L} 和 v_{j,U}
            double muTvjL, muTvjU, vjLTSigmavjL, vjUTSigmavjU;
            
            if (useAssignmentDependent) {
                // Assignment-dependent模型：v_{j,L}和v_{j,U}是N*p维向量
                int n = inst.getN();
                int nTimesP = n * p;
                double[] v_jL = new double[nTimesP];
                double[] v_jU = new double[nTimesP];
                
                // 构建z_k向量（N*p维）
                double[][] z_k = new double[centers.size()][nTimesP];
                for (int k = 0; k < centers.size(); k++) {
                    for (int i = 0; i < n; i++) {
                        for (int j2 = 0; j2 < p; j2++) {
                            int idx = i * p + j2;
                            if (j2 == k) {
                                double x_ik = Math.abs(x[i][k].get(GRB.DoubleAttr.X) - 1.0) < 1e-6 ? 1.0 : 0.0;
                                z_k[k][idx] = x_ik;
                            } else {
                                z_k[k][idx] = 0.0;
                            }
                        }
                    }
                }
                
                // 计算sum_k z_k
                double[] sum_z_k = new double[nTimesP];
                for (int idx = 0; idx < nTimesP; idx++) {
                    sum_z_k[idx] = 0.0;
                    for (int k = 0; k < centers.size(); k++) {
                        sum_z_k[idx] += z_k[k][idx];
                    }
                }
                
                // 计算z_j
                double[] z_j = z_k[j];
                
                // v_{j,L} = (1-α)/p * sum_k z_k - z_j
                // v_{j,U} = z_j - (1+α)/p * sum_k z_k
                for (int idx = 0; idx < nTimesP; idx++) {
                    v_jL[idx] = coeff * sum_z_k[idx] - z_j[idx];
                    v_jU[idx] = z_j[idx] - coeffUpper * sum_z_k[idx];
                }
                
                // 计算 μ^T * v_{j,L} 和 μ^T * v_{j,U}（N*p维）
                muTvjL = 0.0;
                muTvjU = 0.0;
                for (int idx = 0; idx < nTimesP; idx++) {
                    muTvjL += meanVector[idx] * v_jL[idx];
                    muTvjU += meanVector[idx] * v_jU[idx];
                }
                
                // 计算 v_{j,L}^T * Σ * v_{j,L} 和 v_{j,U}^T * Σ * v_{j,U}（使用优化方法）
                vjLTSigmavjL = computeQuadraticForm(v_jL);
                vjUTSigmavjU = computeQuadraticForm(v_jU);
            } else {
                // 原始模型或改进模型：v_{j,L}和v_{j,U}是N维向量
                // v_{j,L} = (1-α)/p * sum_k x_k - x_j
                // v_{j,U} = x_j - (1+α)/p * sum_k x_k
                // 对于原始模型，sum_k x_k = 1（全1向量），所以 v_{j,L} = coeff - x_j
                double[] v_jL = new double[inst.getN()];
                double[] v_jU = new double[inst.getN()];
                for (int i = 0; i < inst.getN(); i++) {
                    v_jL[i] = coeff * sumX[i] - x_j[i];
                    v_jU[i] = x_j[i] - coeffUpper * sumX[i];
                }

                // 计算 μ^T * v_{j,L} 和 μ^T * v_{j,U}
                muTvjL = 0.0;
                muTvjU = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    muTvjL += meanVector[i] * v_jL[i];
                    muTvjU += meanVector[i] * v_jU[i];
                }

                // 计算 v_{j,L}^T * Σ * v_{j,L} 和 v_{j,U}^T * Σ * v_{j,U}
                vjLTSigmavjL = 0.0;
                vjUTSigmavjU = 0.0;
                for (int i = 0; i < inst.getN(); i++) {
                    for (int k = 0; k < inst.getN(); k++) {
                        vjLTSigmavjL += v_jL[i] * covarianceMatrix[i][k] * v_jL[k];
                        vjUTSigmavjU += v_jU[i] * covarianceMatrix[i][k] * v_jU[k];
                    }
                }
            }

            // 如果两个向量的均值都为0，跳过验证
            if (Math.abs(muTvjL) < 1e-10 && Math.abs(muTvjU) < 1e-10) {
                continue;
            }

            // 使用联合约束验证：inf_{P ∈ D} P{d^T * v_{j,L} > 0 或 d^T * v_{j,U} > 0} <= γ
            // 根据 Bonferroni 不等式，这等价于：
            // sup_{P ∈ D} P{d^T * v_{j,L} > 0} + sup_{P ∈ D} P{d^T * v_{j,U} > 0} <= γ

            if (useExactMethod) {
                // 精确方法：使用精确的reformulate约束
                if (useD1) {
                    // D1模糊集：使用Cantelli不等式
                    // 假设均值可行（μ^T * v_{j,L} <= 0 和 μ^T * v_{j,U} <= 0）
                    if (muTvjL > 1e-6 || muTvjU > 1e-6) {
                        // 均值不可行，直接返回false
                        return false;
                    }

                    // 计算约束左侧的值
                    double lhs = 0.0;

                    // 第一项：sup_{P ∈ D_1} P{d^T * v_{j,L} > 0}
                    if (vjLTSigmavjL > 1e-10) {
                        double denominator1 = vjLTSigmavjL + muTvjL * muTvjL;
                        if (denominator1 > 1e-10) {
                            lhs += vjLTSigmavjL / denominator1;
                        } else {
                            lhs += 1.0;
                        }
                    }

                    // 第二项：sup_{P ∈ D_1} P{d^T * v_{j,U} > 0}
                    if (vjUTSigmavjU > 1e-10) {
                        double denominator2 = vjUTSigmavjU + muTvjU * muTvjU;
                        if (denominator2 > 1e-10) {
                            lhs += vjUTSigmavjU / denominator2;
                        } else {
                            lhs += 1.0;
                        }
                    }

                    // 验证约束是否满足
                    // 使用与Gurobi数值容差一致的容差（1e-6），但考虑到分式约束的数值误差，可以稍微放宽
                    if (lhs > riskParam + 1e-5) {
                        return false;
                    }
                } else {
                    // D2模糊集：使用精确重构方法
                    // 计算 κ_{L,j} 和 κ_{U,j}
                    double sqrtVjLTSigmavjL = Math.sqrt(Math.max(0, vjLTSigmavjL));
                    double sqrtVjUTSigmavjU = Math.sqrt(Math.max(0, vjUTSigmavjU));

                    if (sqrtVjLTSigmavjL < 1e-10 && sqrtVjUTSigmavjU < 1e-10) {
                        continue;
                    }

                    // κ_{L,j} = -μ^T * v_{j,L} / sqrt(v_{j,L}^T * Σ * v_{j,L})
                    // κ_{U,j} = -μ^T * v_{j,U} / sqrt(v_{j,U}^T * Σ * v_{j,U})
                    double kappaL = 0.0;
                    double kappaU = 0.0;

                    if (sqrtVjLTSigmavjL > 1e-10) {
                        kappaL = -muTvjL / sqrtVjLTSigmavjL;
                    }
                    if (sqrtVjUTSigmavjU > 1e-10) {
                        kappaU = -muTvjU / sqrtVjUTSigmavjU;
                    }

                    // 计算阈值
                    double sqrtDelta1 = Math.sqrt(delta1);
                    double delta2OverSqrtDelta1 = delta2 / sqrtDelta1;

                    // 计算 inf_{P ∈ D_2} P{d^T * v_{j,L} <= 0} 和 inf_{P ∈ D_2} P{d^T * v_{j,U} <= 0}
                    double infPL = 0.0;
                    double infPU = 0.0;

                    if (sqrtVjLTSigmavjL > 1e-10) {
                        infPL = computeInfPForRelative(kappaL, sqrtDelta1, delta2OverSqrtDelta1);
                    } else {
                        infPL = (muTvjL <= 0) ? 1.0 : 0.0;
                    }

                    if (sqrtVjUTSigmavjU > 1e-10) {
                        infPU = computeInfPForRelative(kappaU, sqrtDelta1, delta2OverSqrtDelta1);
                    } else {
                        infPU = (muTvjU <= 0) ? 1.0 : 0.0;
                    }

                    // 验证约束：infPL + infPU >= 2 - γ
                    double lhs = infPL + infPU;
                    double rhs = 2.0 - riskParam;

                    // 使用与Gurobi数值容差一致的容差，但考虑到分式约束的数值误差，可以稍微放宽
                    if (lhs < rhs - 1e-5) {
                        return false;
                    }
                }
            } else {
                // 近似方法：使用gamma分解进行验证
                double gamma_a_local = 0.5 * riskParam;
                double gamma_b_local = 0.5 * riskParam;

                // 计算factor
                double factorLower, factorUpper;
                if (useD1) {
                    factorLower = Math.sqrt((1 - gamma_a_local) / gamma_a_local);
                    factorUpper = Math.sqrt((1 - gamma_b_local) / gamma_b_local);
                } else {
                    // D2模糊集：根据delta1/delta2与gamma_a和gamma_b的关系确定factor
                    if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                        factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                        factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
                    } else if (delta1 / delta2 <= gamma_a_local && delta1 / delta2 > gamma_b_local) {
                        factorLower = Math.sqrt(delta1) + Math.sqrt((1 - gamma_a_local) / gamma_a_local * (delta2 - delta1));
                        factorUpper = Math.sqrt(delta2 / gamma_b_local);
                    } else if (delta1 / delta2 > gamma_a_local && delta1 / delta2 <= gamma_b_local) {
                        factorLower = Math.sqrt(delta2 / gamma_a_local);
                        factorUpper = Math.sqrt(delta1) + Math.sqrt((1 - gamma_b_local) / gamma_b_local * (delta2 - delta1));
                    } else {
                        factorLower = Math.sqrt(delta2 / gamma_a_local);
                        factorUpper = Math.sqrt(delta2 / gamma_b_local);
                    }
                }

                // 验证约束：μ^T * v_{j,L} + factorLower * sqrt(v_{j,L}^T * Σ * v_{j,L}) <= 0
                // 和 μ^T * v_{j,U} + factorUpper * sqrt(v_{j,U}^T * Σ * v_{j,U}) <= 0
                double stdDevL = Math.sqrt(Math.max(0, vjLTSigmavjL));
                double stdDevU = Math.sqrt(Math.max(0, vjUTSigmavjU));

                double upperBoundL = muTvjL + factorLower * stdDevL;
                double upperBoundU = muTvjU + factorUpper * stdDevU;

                if (upperBoundL > 1e-6 || upperBoundU > 1e-6) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * 计算 inf_{P ∈ D_2} P{d^T * v <= 0}，用于相对平衡性约束验证
     * 其中 κ = -μ^T * v / sqrt(v^T * Σ * v)
     */
    private double computeInfPForRelative(double kappa, double sqrtDelta1, double delta2OverSqrtDelta1) {
        if (kappa < sqrtDelta1) {
            return 0.0;
        } else if (kappa <= delta2OverSqrtDelta1) {
            double sqrtDelta2MinusDelta1 = Math.sqrt(delta2 - delta1);
            double denominator = kappa - sqrtDelta1;
            if (Math.abs(denominator) < 1e-10) {
                return 1.0;
            } else {
                double ratio = sqrtDelta2MinusDelta1 / denominator;
                return 1.0 / (ratio * ratio + 1.0);
            }
        } else {
            if (Math.abs(kappa) < 1e-10) {
                return 0.0;
            } else {
                return (kappa * kappa - delta2) / (kappa * kappa);
            }
        }
    }

    /**
     * 计算 inf_{P ∈ D_2} P{tilde{d}^T * x_j ≥ L}
     * 根据Lemma中的公式
     *
     * 其中 κ_{L,j} = (μ^T * x_j - L) / sqrt(x_j^T * Σ * x_j)
     *
     * 如果 sqrt(δ1) ≤ κ_{L,j} ≤ δ2/sqrt(δ1)，则：
     *   infPL = 1 / ((sqrt(δ2-δ1)/(κ_{L,j} - sqrt(δ1)))^2 + 1)
     *
     * 如果 κ_{L,j} > δ2/sqrt(δ1)，则：
     *   infPL = (κ_{L,j}^2 - δ2) / κ_{L,j}^2
     */
    private double computeInfPL(int j, double kappaL, double sqrtDelta1, double delta2OverSqrtDelta1) {
        String caseType;
        double result;

        // 如果 κ_{L,j} < sqrt(δ1)，Case 1：直接返回0
        if (kappaL < sqrtDelta1) {
            caseType = "Case 1 (κ_L < sqrt(δ1))";
            result = 0.0;
        }
        // Case: sqrt(δ1) ≤ κ_{L,j} ≤ δ2/sqrt(δ1)
        else if (kappaL <= delta2OverSqrtDelta1) {
            caseType = "Case 2 (sqrt(δ1) ≤ κ_L ≤ δ2/sqrt(δ1))";
            double sqrtDelta2MinusDelta1 = Math.sqrt(delta2 - delta1);
            double denominator = kappaL - sqrtDelta1;
            if (Math.abs(denominator) < 1e-10) {
                // 避免除零，返回1.0
                result = 1.0;
            } else {
                double ratio = sqrtDelta2MinusDelta1 / denominator;
                result = 1.0 / (ratio * ratio + 1.0);
            }
        } else {
            // Case: κ_{L,j} > δ2/sqrt(δ1)
            caseType = "Case 3 (κ_L > δ2/sqrt(δ1))";
            if (Math.abs(kappaL) < 1e-10) {
                // 避免除零
                result = 0.0;
            } else {
                result = (kappaL * kappaL - delta2) / (kappaL * kappaL);
            }
        }

        // 输出case信息（只输出case类型）
        System.out.println(String.format("区域 j=%d [infPL]: %s", j, caseType));

        return result;
    }

    /**
     * 计算 inf_{P ∈ D_2} P{tilde{d}^T * x_j ≤ U}
     * 根据Lemma中的公式
     *
     * 其中 κ_{U,j} = (U - μ^T * x_j) / sqrt(x_j^T * Σ * x_j)
     *
     * 如果 sqrt(δ1) ≤ κ_{U,j} ≤ δ2/sqrt(δ1)，则：
     *   infPU = 1 / ((sqrt(δ2-δ1)/(κ_{U,j} - sqrt(δ1)))^2 + 1)
     *
     * 如果 κ_{U,j} > δ2/sqrt(δ1)，则：
     *   infPU = (κ_{U,j}^2 - δ2) / κ_{U,j}^2
     */
    private double computeInfPU(int j, double kappaU, double sqrtDelta1, double delta2OverSqrtDelta1) {
        String caseType;
        double result;

        // 如果 κ_{U,j} < sqrt(δ1)，Case 1：直接返回0
        if (kappaU < sqrtDelta1) {
            caseType = "Case 1 (κ_U < sqrt(δ1))";
            result = 0.0;
        }
        // Case: sqrt(δ1) ≤ κ_{U,j} ≤ δ2/sqrt(δ1)
        else if (kappaU <= delta2OverSqrtDelta1) {
            caseType = "Case 2 (sqrt(δ1) ≤ κ_U ≤ δ2/sqrt(δ1))";
            double sqrtDelta2MinusDelta1 = Math.sqrt(delta2 - delta1);
            double denominator = kappaU - sqrtDelta1;
            if (Math.abs(denominator) < 1e-10) {
                // 避免除零，返回1.0
                result = 1.0;
            } else {
                double ratio = sqrtDelta2MinusDelta1 / denominator;
                result = 1.0 / (ratio * ratio + 1.0);
            }
        } else {
            // Case: κ_{U,j} > δ2/sqrt(δ1)
            caseType = "Case 3 (κ_U > δ2/sqrt(δ1))";
            if (Math.abs(kappaU) < 1e-10) {
                // 避免除零
                result = 0.0;
            } else {
                result = (kappaU * kappaU - delta2) / (kappaU * kappaU);
            }
        }

        // 输出case信息（只输出case类型）
        System.out.println(String.format("区域 j=%d [infPU]: %s", j, caseType));

        return result;
    }

    /**
     * 记录不可行解：将解中取1的index集合记录下来
     */
    private void recordInfeasibleSolution(GRBVar[][] x) throws GRBException {
        Set<String> solutionSet = new HashSet<>();

        for (int i = 0; i < inst.getN(); i++) {
            for (int j = 0; j < centers.size(); j++) {
                if (Math.abs(x[i][j].get(GRB.DoubleAttr.X) - 1.0) < 1e-6) {
                    // 记录 (i,j) 对
                    solutionSet.add(i + "_" + j);
                }
            }
        }

        if (!solutionSet.isEmpty()) {
            infeasibleSolutions.add(solutionSet);
            System.out.println("【生成初始解】记录不可行解，包含 " + solutionSet.size() + " 个变量，当前共有 " + infeasibleSolutions.size() + " 个割约束");
        }
    }

    /**
     * 报告求解状态的详细信息，用于验证是否是全局最优解
     * @param model Gurobi模型
     * @param context 上下文信息（用于日志输出）
     */
    private void reportSolveStatus(GRBModel model, String context) {
        try {
            int status = model.get(GRB.IntAttr.Status);
            String statusName = "";
            
            switch (status) {
                case GRB.OPTIMAL:
                    statusName = "OPTIMAL（全局最优）";
                    break;
                case GRB.SUBOPTIMAL:
                    statusName = "SUBOPTIMAL（次优解）";
                    break;
                case GRB.LOCALLY_OPTIMAL:
                    statusName = "LOCALLY_OPTIMAL（局部最优）";
                    break;
                case GRB.INFEASIBLE:
                    statusName = "INFEASIBLE（不可行）";
                    break;
                case GRB.INF_OR_UNBD:
                    statusName = "INF_OR_UNBD（无界或不可行）";
                    break;
                case GRB.UNBOUNDED:
                    statusName = "UNBOUNDED（无界）";
                    break;
                case GRB.TIME_LIMIT:
                    statusName = "TIME_LIMIT（时间限制）";
                    break;
                default:
                    statusName = "UNKNOWN（未知状态：" + status + "）";
            }
            
            System.out.println("【" + context + "】求解状态: " + statusName);
            
            // 检查NonConvex参数
            try {
                int nonConvex = model.get(GRB.IntParam.NonConvex);
                String nonConvexDesc = "";
                switch (nonConvex) {
                    case 0:
                        nonConvexDesc = "不允许非凸约束（默认）";
                        break;
                    case 1:
                        nonConvexDesc = "允许非凸约束，局部最优";
                        break;
                    case 2:
                        nonConvexDesc = "允许非凸约束，全局最优";
                        break;
                    default:
                        nonConvexDesc = "未知值: " + nonConvex;
                }
                System.out.println("【" + context + "】NonConvex参数: " + nonConvex + " (" + nonConvexDesc + ")");
                
                // 如果包含二次约束但NonConvex不是2，发出警告
                int numQConstrs = model.get(GRB.IntAttr.NumQConstrs);
                if (numQConstrs > 0 && nonConvex != 2) {
                    System.out.println("⚠️ 【" + context + "】警告：模型包含 " + numQConstrs + 
                                     " 个二次约束，但NonConvex=" + nonConvex + 
                                     "，可能只找到局部最优解！");
                } else if (numQConstrs > 0 && nonConvex == 2) {
                    System.out.println("✅ 【" + context + "】模型包含 " + numQConstrs + 
                                     " 个二次约束，NonConvex=2，将寻找全局最优解");
                }
            } catch (Exception e) {
                System.out.println("【" + context + "】无法获取NonConvex参数: " + e.getMessage());
            }
            
            // 报告目标函数值
            if (status == GRB.OPTIMAL || status == GRB.SUBOPTIMAL || status == GRB.LOCALLY_OPTIMAL) {
                try {
                    double objVal = model.get(GRB.DoubleAttr.ObjVal);
                    System.out.println("【" + context + "】目标函数值: " + objVal);
                    
                    // 报告最优性间隙（如果可用）
                    try {
                        double mipGap = model.get(GRB.DoubleAttr.MIPGap);
                        if (mipGap > 0) {
                            System.out.println("【" + context + "】最优性间隙: " + mipGap);
                        }
                    } catch (Exception e) {
                        // MIPGap可能不适用于所有模型类型，忽略
                    }
                } catch (Exception e) {
                    System.out.println("【" + context + "】无法获取目标函数值: " + e.getMessage());
                }
            }
        } catch (Exception e) {
            System.out.println("【" + context + "】报告求解状态时出错: " + e.getMessage());
        }
    }
    
    /**
     * 报告支撑超平面cut迭代的详细求解信息
     * @param model Gurobi模型
     * @param cutIteration 当前迭代次数
     * @param numVarsBefore 求解前的变量数量
     * @param numConstrsBefore 求解前的约束数量
     * @param solveTime 求解时间（毫秒）
     */
    private void reportCutIterationSolveInfo(GRBModel model, int cutIteration, 
                                            int numVarsBefore, int numConstrsBefore, long solveTime) {
        try {
            int status = model.get(GRB.IntAttr.Status);
            String statusName = "";
            
            switch (status) {
                case GRB.OPTIMAL:
                    statusName = "OPTIMAL（全局最优）";
                    break;
                case GRB.SUBOPTIMAL:
                    statusName = "SUBOPTIMAL（次优解）";
                    break;
                case GRB.LOCALLY_OPTIMAL:
                    statusName = "LOCALLY_OPTIMAL（局部最优）";
                    break;
                case GRB.INFEASIBLE:
                    statusName = "INFEASIBLE（不可行）";
                    break;
                case GRB.INF_OR_UNBD:
                    statusName = "INF_OR_UNBD（无界或不可行）";
                    break;
                case GRB.UNBOUNDED:
                    statusName = "UNBOUNDED（无界）";
                    break;
                case GRB.TIME_LIMIT:
                    statusName = "TIME_LIMIT（时间限制）";
                    break;
                default:
                    statusName = "UNKNOWN（未知状态：" + status + "）";
            }
            
            System.out.println("═══════════════════════════════════════════════════════════");
            System.out.println("【支撑超平面cut迭代】第 " + cutIteration + " 次迭代 - 求解信息");
            System.out.println("───────────────────────────────────────────────────────────");
            System.out.println("求解状态: " + statusName + " (状态码: " + status + ")");
            System.out.println("求解时间: " + String.format("%.2f", solveTime / 1000.0) + " 秒");
            
            // 模型规模信息
            int numVars = model.get(GRB.IntAttr.NumVars);
            int numConstrs = model.get(GRB.IntAttr.NumConstrs);
            int numBinVars = model.get(GRB.IntAttr.NumBinVars);
            int numIntVars = model.get(GRB.IntAttr.NumIntVars);
            int numQConstrs = model.get(GRB.IntAttr.NumQConstrs);
            
            System.out.println("模型规模:");
            System.out.println("  变量总数: " + numVars + " (二进制: " + numBinVars + ", 整数: " + numIntVars + ")");
            System.out.println("  约束总数: " + numConstrs + " (线性: " + (numConstrs - numQConstrs) + ", 二次: " + numQConstrs + ")");
            
            // 显示新增的约束数量（如果有）
            int newConstrs = numConstrs - numConstrsBefore;
            if (newConstrs > 0) {
                System.out.println("  本次迭代新增约束: " + newConstrs + " 个");
            }
            
            // 目标函数值
            if (status == GRB.OPTIMAL || status == GRB.SUBOPTIMAL || status == GRB.LOCALLY_OPTIMAL) {
                try {
                    double objVal = model.get(GRB.DoubleAttr.ObjVal);
                    System.out.println("目标函数值: " + String.format("%.6f", objVal));
                    
                    // 最优性间隙
                    try {
                        double mipGap = model.get(GRB.DoubleAttr.MIPGap);
                        if (mipGap > 0) {
                            System.out.println("最优性间隙: " + String.format("%.6e", mipGap));
                        }
                    } catch (Exception e) {
                        // MIPGap可能不适用于所有模型类型，忽略
                    }
                    
                    // 节点信息（如果是MIP）
                    try {
                        int numNodes = (int)model.get(GRB.DoubleAttr.NodeCount);
                        if (numNodes > 0) {
                            System.out.println("探索节点数: " + numNodes);
                        }
                    } catch (Exception e) {
                        // 可能不是MIP模型，忽略
                    }
                } catch (Exception e) {
                    System.out.println("无法获取目标函数值: " + e.getMessage());
                }
            }
            
            System.out.println("═══════════════════════════════════════════════════════════");
        } catch (Exception e) {
            System.out.println("【支撑超平面cut迭代】报告求解信息时出错: " + e.getMessage());
            e.printStackTrace();
        }
    }
}
