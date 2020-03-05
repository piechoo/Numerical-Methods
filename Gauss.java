//Autor: Piotr Piechowicz
public class Gauss {
    public double tab[][] = new double[128][128];
    public double e[]=new double[128];
    public double x1[]=new double[128];
    public double x2[]=new double[128];
    public int cont=0;
    public Gauss(){
    for(int i=0; i<128;i++){//wypełnienie macierzy
        e[i]=1.0;
        x1[i]=0.0;
        tab[i][i]=4.0;
        if(i>3)
            tab[i][i-4]=1.0;
        if(i>0)
            tab[i][i-1]=1.0;
        if(i<124)
            tab[i][i+4]=1.0;
        if(i<127)
            tab[i][i+1]=1.0;
    }}
    double[] gaus(double x1[],double x2[],double e[],double a[][])
    {
        while(true)                     //Uwzgledniam strukturę macierzy przez unikniecie zbednych mnozen przez 0
        {
            x2[0]=(e[0]-a[0][1]*x1[1]-a[0][4]*x1[4])/a[0][0];
            x2[1]=(e[1]-a[1][0]*x2[0]-a[1][2]*x1[2]-a[1][5]*x1[5])/a[1][1];
            x2[2]=(e[2]-a[2][1]*x2[1]-a[2][3]*x1[3]-a[2][6]*x1[6])/a[2][2];
            x2[3]=(e[3]-a[3][2]*x2[2]-a[3][4]*x1[4]-a[3][7]*x1[7])/a[3][3];
            for(int j=4;j<124;j++){
                x2[j]=(e[j]-a[j][j-4]*x2[j-4]-a[j][j-1]*x2[j-1]-a[j][j+1]*x1[j+1]-a[j][j+4]*x1[j+4])/a[j][j];
            }
            x2[124]=(e[124]-a[124][120]*x2[120]-a[124][123]*x2[123]-a[124][125]*x1[125])/a[124][124];
            x2[125]=(e[125]-a[125][121]*x2[121]-a[125][124]*x2[124]-a[125][126]*x1[126])/a[125][125];
            x2[126]=(e[126]-a[126][122]*x2[122]-a[126][125]*x2[125]-a[126][127]*x1[127])/a[126][126];
            x2[127]=(e[127]-a[127][123]*x2[123]-a[127][126]*x2[126])/a[127][127];
            cont++;
            System.out.println(skal(dod(x2,x1,0),(dod(x2,x1,0))));
            if(end(x1,x2)==true)
                break;
            System.arraycopy(x2, 0, x1, 0, x1.length);
        }
        return x2;

    }
    boolean end(double x1[],double x2[]) // sprawdzenie czy wszystkie rozwiązania x(k+1)= x(k)
    {
        int counter=0;
        for (int i=0; i<128;i++)
        {
            if(x1[i]==x2[i])
                counter++;
        }
        if(counter==128)
            return true;
        else
            return false;
    }
    public static void main(String[] args)
    {
        double roz[]=new double [128];
        Gauss g1= new Gauss();
        roz=g1.gaus(g1.x1,g1.x2,g1.e,g1.tab);
        //wypisanie rozwiazan
        for(int i=0;i<128;i++){
            System.out.println("x nr:" + i +" "+ roz[i]);
        }

        //sprawdzenie
        for(int j=4;j<124;j++){//sprawdzenie
            System.out.println((4*roz[j]+roz[j-4]+roz[j-1]+roz[j+1]+roz[j+4]));
        }
        System.out.println("Wykonano "+g1.cont+" krokow");

    }
    double []dod(double a[],double b[],int w)//dodawanie i odejmowanie wektora
    {
        double r[]= new double[a.length];
        for(int i=0;i<a.length;i++)
        {
            if(w>1)
                r[i]=b[i]+a[i];
            else
                r[i]=a[i]-b[i];
        }
        return r;
    }
    double skal(double a[],double b[])//iloczyn skalarny wektorów
    {
        double sum=0.0;
        for(int i=0;i<a.length;i++)
        {
            sum+=a[i]*b[i];
        }
        return sum;
    }
}
