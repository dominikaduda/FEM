import java.util.Scanner;
import java.nio.file.Paths;
import java.nio.file.Path;
import java.io.IOException;

public class Main {

    public static void main(String[] args) {

        int i=0;
        Global_data dataByUser = new Global_data();

        /******************************** TABLICA ELEMENTOW ************************************/
        double iloczyn_double = (dataByUser.nL - 1) * (dataByUser.nH - 1);
        int iloczyn = (int) (iloczyn_double);
        Element[] elements = new Element[iloczyn];
        //Element[] elements = new Element[(dataByUser.nL-1)*(dataByUser.nH-1)];
        for (i = 0; i < ((dataByUser.nL - 1) * (dataByUser.nH - 1)); i++) {
            elements[i] = new Element(dataByUser.nL, dataByUser.nH, dataByUser);

            System.out.println("WEZLY ELEMENTU NR " + i);
            System.out.println(elements[i].nodeID[0]);
            System.out.println(elements[i].nodeID[1]);
            System.out.println(elements[i].nodeID[2]);
            System.out.println(elements[i].nodeID[3]);

            System.out.println("**************************************");
        }

        /****************************TABLICA WEZLOW****************************************/

        double x = 0;
        double y = 0;
        int id = 1;

        Node[] nodes = new Node[(int) (dataByUser.nL * dataByUser.nH)];
        for (i = 0; i < dataByUser.nL * dataByUser.nH; i++) {
            nodes[i] = new Node(x, y, id, dataByUser.L, dataByUser.H);
            if (y < dataByUser.H) {
                y = y + (dataByUser.H / (dataByUser.nH - 1));
            } else {
                y = 0;
                x = x + (dataByUser.L / (dataByUser.nL - 1));
            }

            id++;
        }

        /**************************GRID-pomocniczo*************************/
        Siatka_Grid myGrid = new Siatka_Grid(dataByUser.nL,dataByUser.nH,elements,nodes);


        /************************************* SPRAWDZANIE ID NODE 0 ******************************/

        System.out.println("");
        System.out.println("WSPOLRZEDNE NODE 0");
        System.out.println(nodes[0].id);
        System.out.println("X: " + nodes[0].x);
        System.out.println("Y: " + nodes[0].y);


        /**********************TABELKA GAUSSA Z PUNKTAMI I WAGAMI******************************/

        System.out.println("");
        System.out.println("TABELKA GAUSSA Z PUNKTAMI I WAGAMI");
        for (i = 0; i < dataByUser.gauss.length; i++) {
            for (int j = 0; j < dataByUser.gauss[0].length; j++) {
                System.out.print(dataByUser.gauss[i][j] + " ");
                if (j == dataByUser.gauss[0].length - 1) System.out.println("");
            }
        }

        /*************************TABELA ZE WSPOLRZEDNYMI PUNKTOW ***************************/

        System.out.println("");
        System.out.println("TABELA ZE WSPOLRZEDNYMI PUNKTOW:");

        for (i = 0; i < dataByUser.points.length; i++) {
            for (int j = 0; j < dataByUser.points[0].length; j++) {
                System.out.print(dataByUser.points[i][j] + " ");
                if (j == dataByUser.points[0].length - 1) System.out.println("");
            }

        }

        /*************************ELEMENT UNIWERSALNY **********************************/
        ElementUni elementuni = new ElementUni(dataByUser);
        elementuni.printShapeFunction();
        elementuni.printdNdKSI();
        elementuni.printdNdEta();



        /*****************************JAKOBIAN*******************************************/

        Jakobian[] jakobianElementZero= new Jakobian[dataByUser.PointQuantity*dataByUser.PointQuantity]; //4 bo taka ilosc punktow calkowania
        for(i=0; i<jakobianElementZero.length; i++) {
            jakobianElementZero[i] = new Jakobian(elements, nodes,
                    elementuni, dataByUser, 0, i);
            jakobianElementZero[i].printdNdXdNdXT();
            jakobianElementZero[i].printKsumaDetJ();
            /*jakobianElementZero[i].printMacierzJakobiego();
            jakobianElementZero[i].printdNdY();*/
        }

        for (i =0; i < elements.length; i++ ){
        elements[i].setBigJakobiesMatrix(jakobianElementZero);
        elements[0].printBigJakobiesMatrix();
        elements[0].printBigReverseJakobiesMatrix();
        elements[0].printBigDetJ();
        elements[0].printBigdNdX();
        elements[0].printBigdNdY();
        elements[0].printH();

        jakobianElementZero[0].printC();

        elements[i].setBigC(jakobianElementZero);
        elements[0].printBigC();
        elements[0].printHBC();
        elements[0].printP();}

        Aggregation agr = new Aggregation(dataByUser, elements);
        //agr.printGlobalC();
        //agr.printGlobalH();
        //elements[8].printP();

       //agr.printGlobalP();

        //elements[0].printHBC();

    }

}

