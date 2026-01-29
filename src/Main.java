import com.gurobi.gurobi.GRBException;

import java.io.FileNotFoundException;
import java.io.IOException;

import java.io.File;
import java.io.IOException;

public class Main {


    public static void main(String[] args) throws IOException, GRBException {



//        Instance ins = new Instance("./instances/gen_500.dat");
//        Algo algo = new Algo(ins);
//        algo.run("DU280-05.dat");
        File dir = new File("./instances"); // 获取文件夹路径
        File[] files = dir.listFiles(); // 获取文件夹中的所有文件

        for (File file : files) { // 遍历每个文件
            if (file.isFile()) { // 如果是文件
                String filename = file.getName();
                if (filename.endsWith(".dat")) { // 如果是数据文件
                    System.out.println("------" + filename + "-------");
                    Instance ins = new Instance(file.getPath());
                    int n = ins.getN(); // 获取区域个数
                    Area[] areas = ins.getAreas(); // 获取区域信息
                    int[][] edges = ins.getEdges(); // 获取边集合

                    System.out.println(filename);
                    Algo algo = new Algo(ins);
                    long startTime = System.currentTimeMillis(); // 获取开始时间
                    algo.run(filename);

                    System.out.println();
//                    // 将结果写入输出文件
//                    String outputFilePath = "./output/" + filename.replace(".dat", ".sol");
//                    result.writeToFile(outputFilePath);
                }
            }
        }
    }
}
