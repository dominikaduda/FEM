import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.sqrt;

public class Element {

    boolean area1 = false;
    boolean area2 = false;
    boolean area3 = false;
    boolean area4 = false;

    public int[] nodeID = new int[4];
    public static int ID = 1;
    double[][] bigJakobiesMatrix;
    double[] bigDetJ;
    double[][] bigReverseJakobiesMatrix;
    double[][] bigdNdX;
    double[][] bigdNdY;
    double[][] H;
    double[][] bigC;
    double[] P;
    double[][] H_BC;

    Siatka_Grid grid;
    public Siatka_Grid getGrid() {
        return grid;
    }
    public void setGrid(Siatka_Grid grid) {
        this.grid = grid;
    }

    public Element(double nL, double nH, Global_data globaldata) {
        if (ID % nH == 0) ID++;
        nodeID[0] = ID;
        nodeID[3] = ID + 1;
        nodeID[1] = ID + (int) (nH);
        nodeID[2] = ID + (int) (nH) + 1;

        ID++;
        prepareAreas(nL, nH);
        H_BC = new double[4][4];
        H_BC = prepareHBC( globaldata);
        P = new double[4];
        P=prepareP(globaldata);

    }

    public void prepareAreas(double nL, double nH) {
        for (int i = 0; i < nodeID.length; i++) {
            for (int j = 0; j < nodeID.length; j++) {
                if (nodeID[i] == (1 + nH * j)) area1 = true;
                if (nodeID[i] >= (nH * nL - nH + 1)) area2 = true;
                if (nodeID[i] == (nH * j)) area3 = true;
            }
            if (nodeID[i] <= nH) area4 = true;
        }

    }


    public void setBigJakobiesMatrix(Jakobian[] jakobiesMatrix) {
        List<Double> listOfValuesJakobians = new ArrayList<>();
        List<Double> listOfValuesReverseJakobians = new ArrayList<>();
        bigJakobiesMatrix = new double[jakobiesMatrix.length][4];
        bigReverseJakobiesMatrix = new double[jakobiesMatrix.length][4];
        bigDetJ = new double[jakobiesMatrix.length];

        for (int i = 0; i < jakobiesMatrix.length; i++) {
            for (int j = 0; j < jakobiesMatrix[0].macierzJakobiego.length; j++) {
                for (int g = 0; g < jakobiesMatrix[0].macierzJakobiego[0].length; g++) {
                    listOfValuesJakobians.add(jakobiesMatrix[i].macierzJakobiego[j][g]);
                    listOfValuesReverseJakobians.add(jakobiesMatrix[i].macierzJakobiegoOdwr[j][g]);

                }
            }
            bigDetJ[i] = jakobiesMatrix[0].detJ;
        }

        int index = 0;
        for (int i = 0; i < bigJakobiesMatrix.length; i++) {
            for (int j = 0; j < bigJakobiesMatrix[0].length; j++) {
                bigJakobiesMatrix[i][j] = listOfValuesJakobians.get(index);
                bigReverseJakobiesMatrix[i][j] = listOfValuesReverseJakobians.get(index);
                index++;
            }
        }

        bigdNdX = new double[jakobiesMatrix.length][4];
        ;
        for (int i = 0; i < bigdNdX.length; i++) {
            for (int j = 0; j < bigdNdX[0].length; j++) {
                bigdNdX[i][j] = jakobiesMatrix[i].dNdX[j];
            }
        }

        bigdNdY = new double[jakobiesMatrix.length][4];
        ;
        for (int i = 0; i < bigdNdY.length; i++) {
            for (int j = 0; j < bigdNdY[0].length; j++) {
                bigdNdY[i][j] = jakobiesMatrix[i].dNdY[j];
            }
        }


        H = new double[4][4];
        double suma = 0;
        for (int i = 0; i < H.length; i++)
            for (int j = 0; j < H[0].length; j++) {
                for (int g = 0; g < jakobiesMatrix.length; g++) {
                    suma = suma + jakobiesMatrix[g].KsumadetJ[i][j];
                }
                H[i][j] = suma;
                suma = 0;
            }


    }


    public void setBigC(Jakobian[] jakobiesMatrix) {
        bigC = new double[4][4];
        for (int i = 0; i < bigC.length; i++) {
            for (int j = 0; j < bigC[0].length; j++) {
                for (int g = 0; g < jakobiesMatrix.length; g++) {
                    bigC[i][j] += jakobiesMatrix[g].C[i][j];
                }
            }
        }
    }
    public double[][] prepareHBC(Global_data global_data) {

        double ksi1;
        double ksi2;
        double eta1;
        double eta2;
        double[][] for1pc1 = new double[4][4];
        double[][] for2pc1 = new double[4][4];
        double[][] for1pc2 = new double[4][4];
        double[][] for2pc2 = new double[4][4];
        double[][] for1pc3 = new double[4][4];
        double[][] for2pc3 = new double[4][4];
        double[][] for1pc4 = new double[4][4];
        double[][] for2pc4 = new double[4][4];
        double[][] sum1 = new double[4][4];
        double[][] sum2 = new double[4][4];
        double[][] sum3 = new double[4][4];
        double[][] sum4 = new double[4][4];

        double[][] sum = new double[4][4];

        if (area1) {
            ksi1 = ((-1) / sqrt(3));
            eta1 = -1;
            for1pc1 = fillTable(ksi1, eta1, global_data.alfa);
            ksi2 = (1 / sqrt(3));
            eta2 = -1;
            for2pc1 = fillTable(ksi2, eta2, global_data.alfa);

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {   // detJ = dlugosc/2
                    sum1[i][j] = global_data.L / (global_data.nL - 1) / 2 * (for1pc1[i][j] + for2pc1[i][j]);
                }
            }
        }


        if (area2) {
            ksi1 = 1;
            eta1 = (-1) / sqrt(3);
            for1pc2 = fillTable(ksi1, eta1, global_data.alfa);
            ksi2 = 1;
            eta2 = 1 / sqrt(3);
            for2pc2 = fillTable(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    sum2[i][j] = global_data.H / (global_data.nH - 1) / 2 * (for1pc2[i][j] + for2pc2[i][j]);
                }
            }
        }


        if (area3) {
            ksi1 = 1 / sqrt(3);
            eta1 = 1;
            for1pc3 = fillTable(ksi1, eta1, global_data.alfa);
            ksi2 = (-1) / sqrt(3);
            eta2 = 1;
            for2pc3 = fillTable(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    sum3[i][j] = global_data.L / (global_data.nL - 1) / 2 * (for1pc3[i][j] + for2pc3[i][j]);
                }
            }
        }
        if (area4) {
            ksi1 = -1;
            eta1 = 1 / sqrt(3);
            for1pc4 = fillTable(ksi1, eta1, global_data.alfa);
            ksi2 = -1;
            eta2 = -1 / sqrt(3);
            for2pc4 = fillTable(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    sum4[i][j] = global_data.H / (global_data.nH - 1) / 2 * (for1pc4[i][j] + for2pc4[i][j]);
                }
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                sum[i][j] = sum1[i][j] + sum2[i][j] + sum3[i][j] + sum4[i][j];
            }
        }
        return sum;
    }


    public double[][] fillTable(double ksi, double eta, double alfa) {
        double[] helpShape = new double[4];
        double[][] filled = new double[4][4];
        helpShape[0] = 0.25 * (1 - ksi) * (1 - eta);
        helpShape[1] = 0.25 * (1 + ksi) * (1 - eta);
        helpShape[2] = 0.25 * (1 + ksi) * (1 + eta);
        helpShape[3] = 0.25 * (1 - ksi) * (1 + eta);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                filled[i][j] = helpShape[i] * helpShape[j] * alfa;
            }
        }
        return filled;
    }

    public double[] prepareP(Global_data global_data) {

        double ksi1;
        double ksi2;
        double eta1;
        double eta2;

        double[] P1_1 = new double[4];
        double[] P2_1 = new double[4];
        double[] P3_1 = new double[4];
        double[] P4_1 = new double[4];
        double[] P1_2 = new double[4];
        double[] P2_2 = new double[4];
        double[] P3_2 = new double[4];
        double[] P4_2 = new double[4];
        double[] sum1 = new double[4];
        double[] sum2 = new double[4];
        double[] sum3 = new double[4];
        double[] sum4 = new double[4];

        double[] P = new double[4];

        if (area1) {
            ksi1 = ((-1) / sqrt(3));
            eta1 = -1;
            P1_1 = fillP(ksi1, eta1, global_data.alfa);
            ksi2 = (1 / sqrt(3));
            eta2 = -1;
            P1_2 = fillP(ksi2, eta2, global_data.alfa);

            for (int i = 0; i < 4; i++) {
                sum1[i] = global_data.L / (global_data.nL - 1) / 2 * (P1_1[i] + P1_2[i]);
            }
        }


        if (area2) {
            ksi1 = 1;
            eta1 = (-1) / sqrt(3);
            P2_1 = fillP(ksi1, eta1, global_data.alfa);
            ksi2 = 1;
            eta2 = 1 / sqrt(3);
            P2_2= fillP(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                sum2[i] = global_data.L / (global_data.nL - 1) / 2 * (P2_1[i] + P2_2[i]);

            }
        }


        if (area3) {
            ksi1 = 1 / sqrt(3);
            eta1 = 1;
            P3_1 = fillP(ksi1, eta1, global_data.alfa);
            ksi2 = (-1) / sqrt(3);
            eta2 = 1;
            P3_2 = fillP(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                sum3[i] = global_data.L / (global_data.nL - 1) / 2 * (P3_1[i] + P3_2[i]);

            }
        }

        if (area4) {
            ksi1 = -1;
            eta1 = 1 / sqrt(3);
            P4_1 = fillP(ksi1, eta1, global_data.alfa);
            ksi2 = -1;
            eta2 = -1 / sqrt(3);
            P4_2= fillP(ksi2, eta2, global_data.alfa);
            for (int i = 0; i < 4; i++) {
                sum4[i] = global_data.L / (global_data.nL - 1) / 2 * (P4_1[i] + P4_2[i]);

            }
        }

        for (int i = 0; i < 4; i++) P[i] = (sum1[i] + sum2[i] + sum3[i] + sum4[i])*global_data.ambTemp;
        return P;
    }


    public double[] fillP(double ksi, double eta, double alfa) {
        double[] helpShape = new double[4];
        double[] filled = new double[4];
        helpShape[0] = 0.25 * (1 - ksi) * (1 - eta);
        helpShape[1] = 0.25 * (1 + ksi) * (1 - eta);
        helpShape[2] = 0.25 * (1 + ksi) * (1 + eta);
        helpShape[3] = 0.25 * (1 - ksi) * (1 + eta);
        for (int i = 0; i < 4; i++) {
            filled[i] = helpShape[i] * alfa;
        }
        return filled;
    }

    /*********************** printowanie *****************************/


    public void printBigJakobiesMatrix() {
        System.out.println();
        System.out.println("Duza macierz Jakobiego");
        for (int i = 0; i < bigJakobiesMatrix.length; i++) {
            for (int j = 0; j < bigJakobiesMatrix[0].length; j++) {
                System.out.print(bigJakobiesMatrix[i][j] + " ");
                if (j == bigJakobiesMatrix[0].length - 1) System.out.println();
            }
        }
    }


    public void printBigReverseJakobiesMatrix() {
        System.out.println();
        System.out.println("Duza odwrotna macierz Jakobiego");
        for (int i = 0; i < bigReverseJakobiesMatrix.length; i++) {
            for (int j = 0; j < bigReverseJakobiesMatrix[0].length; j++) {
                System.out.print(bigReverseJakobiesMatrix[i][j] + " ");
                if (j == bigReverseJakobiesMatrix[0].length - 1) System.out.println();
            }
        }
    }

    public void printBigDetJ() {
        System.out.println();
        System.out.println("Macierz wyznaczników");

        for (int i = 0; i < bigDetJ.length; i++)
            System.out.print(bigDetJ[i] + " ");

    }

    public void printBigdNdX() {
        System.out.println();
        System.out.println("Duże dNdX");
        for (int i = 0; i < bigdNdX.length; i++) {
            for (int j = 0; j < bigdNdX[0].length; j++) {
                System.out.print(bigdNdX[i][j] + "   ");
                if (j == (bigdNdX[0].length - 1))
                    System.out.println();
            }
        }
    }

    public void printBigdNdY() {
        System.out.println();
        System.out.println("Duże dNdY");
        for (int i = 0; i < bigdNdY.length; i++) {
            for (int j = 0; j < bigdNdY[0].length; j++) {
                System.out.print(bigdNdY[i][j] + "   ");
                if (j == (bigdNdY[0].length - 1))
                    System.out.println();
            }
        }
    }

    public void printH() {
        System.out.println();
        System.out.println("Macierz H ");
        for (int i = 0; i < H.length; i++) {
            for (int j = 0; j < H[0].length; j++) {
                System.out.print(H[i][j] + "   ");
                if (j == (H[0].length - 1))
                    System.out.println();
            }
        }
    }

    public void printBigC() {
        System.out.println();
        System.out.println("Macierz big C ");
        for (int i = 0; i < bigC.length; i++) {
            for (int j = 0; j < bigC[0].length; j++) {
                System.out.print(bigC[i][j] + "   ");
                if (j == (bigC[0].length - 1))
                    System.out.println();
            }
        }
    }

    public void printHBC() {
        System.out.println();
        System.out.println("Macierz H_BC");
        for (int i = 0; i < H_BC[0].length; i++) {
            for (int j = 0; j < H_BC.length; j++) {
                System.out.print(H_BC[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }

    public void printP(){
        System.out.println();
        System.out.println("Wektor P: ");

        for (int i =0 ; i< P.length; i++)
            System.out.print(P[i] + " ");
        System.out.println();

    }

}
