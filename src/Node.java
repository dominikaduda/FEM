public class Node {
    double x;
    double y;

    //double t;
    int id;
    int bc;


    public Node (double x, double y, int id, double width, double height)
    {   this.id=id;
        this.x=x;
        this.y=y;
        setBc(width, height);
    }

    private void setBc(double width, double height){

        //most lower / most upper / most left / most right edge
        if(this.x==0 || this.x==width || this.y==0 || this.y==height) this.bc=1;
        else this.bc=0;
    }

    public double getX(){
        return x;
    }

    public double getY(){
        return y;
    }

    public int getBc() {
        return bc;
    }

}
