public class Lab2Variant5 {

    // Аналитическое решение φ(x) = 2 / (1 + x^2)
    private static double exactSolution(double x) {
        return 2.0 / (1.0 + x * x);
    }

    // Правая часть для y'' = f(x, y, y')
    // Выведено из исходного уравнения yy' + 2x(y')^2 = xyy''
    private static double secondDerivative(double x, double y, double yp) {
        return (y * yp + 2.0 * x * yp * yp) / (x * y);
    }

    // Один шаг метода Рунге–Кутты 4-го порядка для системы:
    // y1' = y2
    // y2' = f(x, y1, y2)
    private static double[] rk4Step(double x, double y, double yp, double h) {
        // k1
        double k1y = h * yp;
        double k1yp = h * secondDerivative(x, y, yp);

        // k2
        double k2y = h * (yp + 0.5 * k1yp);
        double k2yp = h * secondDerivative(
                x + 0.5 * h,
                y + 0.5 * k1y,
                yp + 0.5 * k1yp
        );

        // k3
        double k3y = h * (yp + 0.5 * k2yp);
        double k3yp = h * secondDerivative(
                x + 0.5 * h,
                y + 0.5 * k2y,
                yp + 0.5 * k2yp
        );

        // k4
        double k4y = h * (yp + k3yp);
        double k4yp = h * secondDerivative(
                x + h,
                y + k3y,
                yp + k3yp
        );

        double yNext = y + (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0;
        double ypNext = yp + (k1yp + 2.0 * k2yp + 2.0 * k3yp + k4yp) / 6.0;

        return new double[] { yNext, ypNext };
    }

    public static void main(String[] args) {
        // yy' + 2x (y')^2 = xyy'', y(1)=1, y'(1)=-1

        double x0 = 1.0;
        double y0 = 1.0;
        double yp0 = -1.0;

        double b = 2.0;   // правый конец [x0, b]
        int n = 10;       // число шагов

        double h = (b - x0) / n;

        double x = x0;
        double y = y0;
        double yp = yp0;

        double maxAbsError = 0.0;
        double maxRelError = 0.0;

        System.out.printf("%-4s %-10s %-15s %-15s %-15s %-15s%n",
                "k", "x_k", "y_num", "y_exact", "absError", "relError");

        for (int k = 0; k <= n; k++) {
            double yExact = exactSolution(x);
            double absError = Math.abs(y - yExact);
            double relError = Math.abs(absError / yExact);

            if (absError > maxAbsError) {
                maxAbsError = absError;
            }
            if (relError > maxRelError) {
                maxRelError = relError;
            }

            System.out.printf("%-4d %-10.5f %-15.8f %-15.8f %-15.8e %-15.8e%n",
                    k, x, y, yExact, absError, relError);

            if (k < n) {
                double[] next = rk4Step(x, y, yp, h);
                y = next[0];
                yp = next[1];
                x = x0 + (k + 1) * h;
            }
        }

        System.out.println();
        System.out.printf("Max absolute error   Δ = %.8e%n", maxAbsError);
        System.out.printf("Max relative error   δ = %.8e%n", maxRelError);
    }
}
