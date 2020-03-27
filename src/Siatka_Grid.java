public class Siatka_Grid {
    double tmp_nH;
    double tmp_nL;
    double ne=(tmp_nL-1)*(tmp_nH-1);
    double nh=tmp_nL*tmp_nH;
    Element[] elementGrid;
    Node[] nodeGrid;

    public Siatka_Grid(double nL, double nH, Element[] element, Node[] node)
    {
        tmp_nH=nH;
        tmp_nL=nL;
        elementGrid=element;
        nodeGrid=node;
    }
}
