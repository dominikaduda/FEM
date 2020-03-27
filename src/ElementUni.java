import java.awt.*;

public class ElementUni {

    double ksi;
    double eta;
    public double[][] shapeFunction;
    public double[][] dNdKsi;
    public double[][] dNdEta;
    int PointQuantity;

    public ElementUni(Global_data global_data) {
        PointQuantity = global_data.PointQuantity;
        shapeFunction = new double[PointQuantity * PointQuantity][4];
        dNdKsi = new double[PointQuantity * PointQuantity][4];
        dNdEta = new double[PointQuantity * PointQuantity][4];
        /*********************WYPELNIANIE FUNKCJI KSZTALTU*******************/
        for (int i = 0; i < shapeFunction.length; i++) {
            ksi = global_data.points[i][0];
            eta = global_data.points[i][1];
            for (int j = 0; j < shapeFunction[0].length /*4*/ ; j++) {

                if (j == 0) shapeFunction[i][j] = 0.25 * (1 - ksi) * (1 - eta);
                if (j == 1) shapeFunction[i][j] = 0.25* (1 + ksi) * (1 - eta);
                if (j == 2) shapeFunction[i][j] = 0.25 * (1 + ksi) * (1 + eta);
                if (j == 3) shapeFunction[i][j] = 0.25  * (1 - ksi) * (1 + eta);

            }
        }

        /***************WYPELNIANIE TABLICY POCHODNYCH PO KSI********************/
        for ( int i = 0; i < dNdKsi.length; i++) {
            for (int j = 0; j < dNdKsi[0].length /*4*/; j++) {
                if (j == 0) dNdKsi[i][j] = -0.25 * (1 - global_data.points[i][1]);
                if (j == 1) dNdKsi[i][j] = 0.25 * (1 - global_data.points[i][1]);
                if (j == 2) dNdKsi[i][j] = 0.25 * (1 + global_data.points[i][1]);
                if (j == 3) dNdKsi[i][j] = -0.25 * (1 + global_data.points[i][1]);
            }
        }

        /****************WYPELNIANIE TABLICY POCHODNYCH PO ETA********************/
        for ( int i = 0; i < dNdEta.length; i++) {
            for (int j = 0; j < dNdEta[0].length /*4*/; j++) {
                if (j == 0) dNdEta[i][j] = -0.25 * (1 - global_data.points[i][0]);
                if (j == 1) dNdEta[i][j] = -0.25 * (1 + global_data.points[i][0]);
                if (j == 2) dNdEta[i][j] = 0.25 * (1+ global_data.points[i][0]);
                if (j == 3) dNdEta[i][j] = 0.25 * (1- global_data.points[i][0]);
            }
        }



    }
    public void printShapeFunction(){
        System.out.println("");
        System.out.println("SHAPE FUNCTIONS:");
        for (int i = 0; i < shapeFunction.length; i++) {
            for (int j = 0; j < shapeFunction[0].length /*4*/ ; j++) {
                System.out.print(shapeFunction[i][j] + "   ");
                if (j == shapeFunction[0].length -1 ) System.out.println("");
            }

        }
    }

    public void printdNdKSI(){
        System.out.println("");
        System.out.println("dNdKSI: ");
        for (int i = 0; i < dNdKsi.length; i++) {
            for (int j = 0; j < dNdKsi[0].length /*4*/ ; j++) {
                System.out.print(dNdKsi[i][j] + "   ");
                if (j == dNdKsi[0].length -1 ) System.out.println("");
            }

        }
    }

    public void printdNdEta(){
        System.out.println("");
        System.out.println("dNdEta: ");
        for (int i = 0; i < dNdEta.length; i++) {
            for (int j = 0; j < dNdEta[0].length /*4*/ ; j++) {
                System.out.print(dNdEta[i][j] + "   ");
                if (j == dNdEta[0].length -1 ) System.out.println("");
            }

        }


    }
}
