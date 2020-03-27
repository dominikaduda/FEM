import static java.lang.Math.abs;

public class JacobiMethod {

    private static Double A[][];
    private static double M[][];
    private static double N[];
    private static double b[];
    private static double x1[];
    private static double x2[];

    public double[] getX1() {
        return x1;
    }

    public JacobiMethod(Global_data global_data, double[][] H, double[] P) {
        jacobiMethod(global_data, H, P);
    }

    public void jacobiMethod(Global_data global_data, double[][] H, double[] P) {

        int num = (int) ((global_data.nH) * (global_data.nH));
        int iter;
        int i, j, k;

        A = new Double[num][num];
        M = new double[num][num];
        N = new double[num];
        b = new double[num];
        x1 = new double[num];
        x2 = new double[num];

        // Get values of A
        for (i = 0; i < num; i++)
            for (j = 0; j < num; j++) {
                A[i][j] = H[i][j];
                if ((i == j) && (A[i][j] == 0)) {
                    System.out.println("Wartosci na przekatnej musza byc rozne od 0");
                    return;
                }
            }

// Get values of b
        for (i = 0; i < num; i++) {
            b[i] = P[i];
        }

// Calculate N = D^-1
        for (i = 0; i < num; i++)
            N[i] = 1 / A[i][i];

// Calculate M = -D^-1 (L + U)
        for (i = 0; i < num; i++)
            for (j = 0; j < num; j++)
                if (i == j)
                    M[i][j] = 0;
                else
                    M[i][j] = -(A[i][j] * N[i]);

//Initialize x
        for (i = 0; i < num; i++)
            x1[i] = 0;
        iter = 100;

        for (k = 0; k < iter; k++) {
            for (i = 0; i < num; i++) {
                x2[i] = N[i] * b[i];
                for (j = 0; j < num; j++)
                    x2[i] += M[i][j] * x1[j];
            }
            for (i = 0; i < num; i++) {
                x1[i] = x2[i];
                if (abs(x1[i] = x2[i]) < 0.001) iter = k - 1;
            }
        }
    }
}


