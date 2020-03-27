import static java.lang.Math.sqrt;

public class Jakobian {

    double dYdEta;
    double dYdKsi;
    double dXdKsi;
    double dXdEta;
    double[][] macierzJakobiego;
    double detJ;
    int numberElement;
    int punktCalkowania;
    double[][] macierzJakobiegoOdwr;
    double dNdX[];
    double dNdY[];
    double dNdXdNdXT[][];
    double dNdYdNdYT[][];
    double KsumadetJ[][];
    double C[][];
    //double H_BC[][];
    //double[] P;

    public Jakobian(Element[] elements, Node[] nodes, ElementUni elementuni,
                    Global_data globaldata, int numberElement, int punktCalkowania) {

        this.punktCalkowania = punktCalkowania;
        this.numberElement = numberElement;
        dYdEta = 0;
        dXdKsi = 0;
        dXdEta = 0;
        dYdKsi = 0;


        for (int i = 0; i < 4 /*(globaldata.PointQuantity * globaldata.PointQuantity)*/; i++) {
            dYdEta = dYdEta + nodes[elements[numberElement].nodeID[i]].getY() * elementuni.dNdEta[punktCalkowania][i];
            dYdKsi = dYdKsi + nodes[elements[numberElement].nodeID[i]].getY() * elementuni.dNdKsi[punktCalkowania][i];
            dXdKsi = dXdKsi + nodes[elements[numberElement].nodeID[i]].getX() * elementuni.dNdKsi[punktCalkowania][i];
            dXdEta = dXdEta + nodes[elements[numberElement].nodeID[i]].getX() * elementuni.dNdEta[punktCalkowania][i];
        }

        macierzJakobiego = new double[2][2];
        macierzJakobiego[0][0] = dXdKsi;
        macierzJakobiego[0][1] = dYdKsi;
        macierzJakobiego[1][0] = dXdEta;
        macierzJakobiego[1][1] = dYdEta;

        detJ = macierzJakobiego[0][0] * macierzJakobiego[1][1] - macierzJakobiego[0][1] * macierzJakobiego[1][0];

        //MACIERZ DOPELNIEN
        macierzJakobiegoOdwr = new double[2][2];
        macierzJakobiegoOdwr[0][0] = (1 / detJ) * dYdEta;
        macierzJakobiegoOdwr[0][1] = (1 / detJ) * (-1) * dYdKsi;
        macierzJakobiegoOdwr[1][0] = (1 / detJ) * (-1) * dXdEta;
        macierzJakobiegoOdwr[1][1] = (1 / detJ) * dXdKsi;

        //DNdX dla punktu calkowania
        dNdX = new double[4];
        for (int i = 0; i < dNdX.length; i++) {
            dNdX[i] = (1 / detJ) * (dYdEta * elementuni.dNdKsi[punktCalkowania][i] + (-1) * dXdEta
                    * elementuni.dNdEta[punktCalkowania][i]);
        }
        //DNdY dla punktu calkowania
        dNdY = new double[4];
        for (int i = 0; i < dNdY.length; i++) {
            dNdY[i] = (1 / detJ) * ((-1) * dXdEta * elementuni.dNdKsi[punktCalkowania][i] + dXdKsi
                    * elementuni.dNdEta[punktCalkowania][i]);
        }
        //dNdX*transponowane dNdX dla punktu calkowania
        dNdXdNdXT = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                dNdXdNdXT[i][j] = dNdX[i] * dNdX[j];
            }
        }
        //dNdY*transponowane dNdY dla punktu calkowania
        dNdYdNdYT = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                dNdYdNdYT[i][j] = dNdY[i] * dNdY[j];
            }
        }

        //MACIERZ H
        KsumadetJ = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                KsumadetJ[i][j] = globaldata.K * (dNdXdNdXT[i][j] * globaldata.gauss[1][0] + dNdYdNdYT[i][j] * globaldata.gauss[1][1]) * detJ;
            }
        }


        //macierz C
        C = new double[4][4];
        for (int i = 0; i < elementuni.shapeFunction[0].length; i++) {
            for (int j = 0; j < elementuni.shapeFunction.length; j++) {

                C[i][j] = detJ * globaldata.c * globaldata.ro
                        * elementuni.shapeFunction[punktCalkowania][i] * elementuni.shapeFunction[punktCalkowania][j];
            }
        }

    }

 //************************************* tu juz samo wypisywanie ************************************//
    public void printMacierzJakobiego() {

        System.out.println();
        System.out.println("Element nr: " + numberElement);
        System.out.println("W punkcie calkowania nr: " + punktCalkowania);

        System.out.println("Macierz Jakobiego");

        for (int i = 0; i < macierzJakobiego.length; i++) {
            for (int j = 0; j < macierzJakobiego[0].length; j++) {
                System.out.print(macierzJakobiego[i][j] + " ");
                if (j == 1) System.out.println();
            }
        }
        System.out.println();
        System.out.println("DET J");
        System.out.println(detJ);
    }

    public void printMacierzJakobiegoOdwr() {

        System.out.println();
        System.out.println("Macierz Jakobiego Odwrotna");

        for (int i = 0; i < macierzJakobiegoOdwr.length; i++) {
            for (int j = 0; j < macierzJakobiegoOdwr[0].length; j++) {
                System.out.print(macierzJakobiegoOdwr[i][j] + " ");
                if (j == 1) System.out.println();
            }
        }
    }

    public void printdNdX() {
        System.out.println("dNdX dla " + punktCalkowania + " punktu calkowania");
        for (int i = 0; i < dNdX.length; i++)
            System.out.print(dNdX[i]);
    }

    public void printdNdY() {
        System.out.println("dNdY dla " + punktCalkowania + " punktu calkowania");
        for (int i = 0; i < dNdY.length; i++)
            System.out.print(dNdY[i]);

    }

    public void printdNdXdNdXT() {
        System.out.println();
        System.out.println("dNdX * transponowana dNdX dla punktu calkowania nr:  " + punktCalkowania);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.print(dNdXdNdXT[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }

    public void printdNdYdNdYT() {
        System.out.println();
        System.out.println("dNdY * transponowana dNdY dla punktu calkowania nr:  " + punktCalkowania);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.print(dNdYdNdYT[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }

    public void printKsumaDetJ() {
        System.out.println();
        System.out.println("K*( dNdX* transponowane dNdX + dNdY* transponowane dNdY ) * det J" +
                " dla punktu calkowania nr:  " + punktCalkowania);
        for (int i = 0; i < KsumadetJ[0].length; i++) {
            for (int j = 0; j < KsumadetJ.length; j++) {
                System.out.print(KsumadetJ[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }

    public void printC() {
        System.out.println();
        System.out.println("Macierz C dla punktu calkowania nr " + punktCalkowania);
        for (int i = 0; i < C[0].length; i++) {
            for (int j = 0; j < C.length; j++) {
                System.out.print(C[i][j] + " ");
                if (j == 3) System.out.println();
            }
        }
    }



    public double[][] getMacierzJakobiego() {
        return macierzJakobiego;
    }

    public double[][] getMacierzJakobiegoOdwr() {
        return macierzJakobiegoOdwr;
    }

    public double getDetJ() {
        return detJ;
    }

    public double[][] getC() {
        return C;
    }



}
