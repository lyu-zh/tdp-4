import java.lang.reflect.Constructor;
import java.lang.reflect.Method;

/**
 * COPT Java API 自检工具（反射版）。
 *
 * 目的：
 * 1) 检查是否能加载 COPT Java 类；
 * 2) 检查是否能创建环境与模型（间接验证 license 与 native 库）；
 * 3) 尝试求解一个最小模型。
 *
 * 说明：
 * - 使用反射，避免在未配置 COPT jar 时编译失败。
 * - 若你已经把 COPT jar 配好，这个类可直接运行给出诊断信息。
 */
public class CoptApiSelfCheck {

    public static void main(String[] args) {
        System.out.println("=== COPT Java API Self Check ===");
        System.out.println("java.version = " + System.getProperty("java.version"));
        System.out.println("java.library.path = " + System.getProperty("java.library.path"));
        System.out.println("COPT_HOME = " + System.getenv("COPT_HOME"));
        System.out.println("LD_LIBRARY_PATH = " + System.getenv("LD_LIBRARY_PATH"));
        System.out.println("PATH contains copt? = " + containsIgnoreCase(System.getenv("PATH"), "copt"));

        try {
            // COPT 8 的 Java 包名通常为 copt.*
            Class<?> constsClass = tryLoadClass("copt.Consts");
            if (constsClass == null) {
                // 兼容某些发行版使用 com.copt.*
                constsClass = tryLoadClass("com.copt.COPT");
            }
            if (constsClass == null) {
                throw new RuntimeException("未找到 COPT 常量类（copt.Consts 或 com.copt.COPT），请先把 copt.jar 加入 classpath。");
            }
            System.out.println("[OK] Loaded constants class: " + constsClass.getName());

            Class<?> envClass = firstNonNull(tryLoadClass("copt.Envr"), tryLoadClass("com.copt.Envr"));
            Class<?> modelClass = firstNonNull(tryLoadClass("copt.Model"), tryLoadClass("com.copt.Model"));
            Class<?> exprClass = firstNonNull(tryLoadClass("copt.Expr"), tryLoadClass("com.copt.Expr"));
            Class<?> intAttrClass = firstNonNull(tryLoadClass("copt.IntAttr"), tryLoadClass("com.copt.IntAttr"));
            Class<?> dblAttrClass = firstNonNull(tryLoadClass("copt.DblAttr"), tryLoadClass("com.copt.DblAttr"));
            Class<?> dblInfoClass = firstNonNull(tryLoadClass("copt.DblInfo"), tryLoadClass("com.copt.DblInfo"));
            if (envClass == null || modelClass == null || exprClass == null) {
                throw new RuntimeException("COPT 关键类缺失（Envr/Model/Expr），请检查 jar 版本。");
            }
            if (intAttrClass == null || dblAttrClass == null || dblInfoClass == null) {
                throw new RuntimeException("COPT 属性类缺失（IntAttr/DblAttr/DblInfo），请检查 jar 版本。");
            }
            System.out.println("[OK] Loaded classes: Envr/Model/Expr");

            // new Envr()
            Object env = newInstance(envClass);
            System.out.println("[OK] Created Envr");

            // env.createModel("smoke")
            Method createModel = envClass.getMethod("createModel", String.class);
            Object model = createModel.invoke(env, "smoke");
            System.out.println("[OK] Created Model");

            // 常量（反射获取）
            char CONTINUOUS = getStaticChar(constsClass, "CONTINUOUS");
            double INFINITY = getStaticDouble(constsClass, "INFINITY");
            char LESS_EQUAL = getStaticChar(constsClass, "LESS_EQUAL");
            int MAXIMIZE = getStaticInt(constsClass, "MAXIMIZE");

            Class<?> varClass = firstNonNull(tryLoadClass("copt.Var"), tryLoadClass("com.copt.Var"));
            if (varClass == null) {
                throw new RuntimeException("未找到 Var 类（copt.Var/com.copt.Var），请检查 COPT 版本。");
            }

            // model.addVar(lb, ub, obj, vtype, name)
            Method addVar = modelClass.getMethod("addVar", double.class, double.class, double.class, char.class, String.class);
            Object x = addVar.invoke(model, 0.0, INFINITY, 0.0, CONTINUOUS, "x");
            Object y = addVar.invoke(model, 0.0, INFINITY, 0.0, CONTINUOUS, "y");

            // lhs = x + 2y
            Object lhs = newInstance(exprClass);
            // COPT 8: Expr.addTerm(Var, double)
            Method addTerm = exprClass.getMethod("addTerm", varClass, double.class);
            addTerm.invoke(lhs, x, 1.0);
            addTerm.invoke(lhs, y, 2.0);

            // model.addConstr(lhs, <=, 4.0, "c1")
            Method addConstr = modelClass.getMethod("addConstr", exprClass, char.class, double.class, String.class);
            addConstr.invoke(model, lhs, LESS_EQUAL, 4.0, "c1");

            // obj = 3x + y
            Object obj = newInstance(exprClass);
            addTerm.invoke(obj, x, 3.0);
            addTerm.invoke(obj, y, 1.0);

            Method setObjective = modelClass.getMethod("setObjective", exprClass, int.class);
            setObjective.invoke(model, obj, MAXIMIZE);

            // solve
            Method solve = modelClass.getMethod("solve");
            solve.invoke(model);
            System.out.println("[OK] Model solved");

            // 输出状态和结果
            Method getIntAttr = modelClass.getMethod("getIntAttr", String.class);
            Method getDblAttr = modelClass.getMethod("getDblAttr", String.class);
            Method getVarInfo = varClass.getMethod("get", String.class);
            String lpStatusName = (String) intAttrClass.getField("LpStatus").get(null);
            String lpObjName = (String) dblAttrClass.getField("LpObjVal").get(null);
            String valName = (String) dblInfoClass.getField("Value").get(null);
            int status = (Integer) getIntAttr.invoke(model, lpStatusName);
            double objVal = (Double) getDblAttr.invoke(model, lpObjName);
            double xVal = (Double) getVarInfo.invoke(x, valName);
            double yVal = (Double) getVarInfo.invoke(y, valName);

            System.out.println("Status = " + status);
            System.out.println("ObjVal = " + objVal);
            System.out.println("x = " + xVal);
            System.out.println("y = " + yVal);

            // dispose
            disposeIfExists(model, modelClass);
            disposeIfExists(env, envClass);

            System.out.println("[PASS] COPT Java API 基本可用。");
        } catch (Throwable t) {
            System.out.println("[FAIL] COPT Java API 自检失败：");
            t.printStackTrace(System.out);
            System.exit(1);
        }
    }

    private static Class<?> tryLoadClass(String name) {
        try {
            return Class.forName(name);
        } catch (Throwable ignored) {
            return null;
        }
    }

    private static Object newInstance(Class<?> cls) throws Exception {
        Constructor<?> ctor = cls.getDeclaredConstructor();
        ctor.setAccessible(true);
        return ctor.newInstance();
    }

    private static int getStaticInt(Class<?> cls, String field) throws Exception {
        return cls.getField(field).getInt(null);
    }

    private static char getStaticChar(Class<?> cls, String field) throws Exception {
        return cls.getField(field).getChar(null);
    }

    private static double getStaticDouble(Class<?> cls, String field) throws Exception {
        return cls.getField(field).getDouble(null);
    }

    private static <T> T firstNonNull(T a, T b) {
        return a != null ? a : b;
    }

    private static void disposeIfExists(Object obj, Class<?> cls) {
        try {
            Method dispose = cls.getMethod("dispose");
            dispose.invoke(obj);
        } catch (Throwable ignored) {
        }
    }

    private static boolean containsIgnoreCase(String s, String key) {
        if (s == null || key == null) return false;
        return s.toLowerCase().contains(key.toLowerCase());
    }
}

