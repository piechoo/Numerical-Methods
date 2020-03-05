//Autor:Piotr Piechowicz
import java.lang.Math;
public class Eigenvalues {


    public static void main(String[]args)
    {

        double[][] A1 = { { 19, 13, 10, 10, 13, -17 },
                { 13, 13, 10, 10, -11, 13 },
                { 10, 10, 10, -2, 10, 10 },
                { 10, 10, -2, 10, 10, 10 },
                { 13, -11, 10, 10, 13, 13 },
                { -17, 13, 10, 10, 13, 19 },
        };
        double []y1=new double[A1.length];
        y1[0]=1;
       
        double []y2 = { 1, 1, 1, 1, 1, 1 };
        double []wyn=wart_wlasna(A1,y1,y2);
        System.out.println("\nWartosc wlasna 1= "+(wyn[0]/12));
        System.out.println("Wartosc wlasna 2= "+(wyn[1]/12));

    }


    static double skal(double a[],double b[])//obliczanie iloczynu skalarnego wektorów
    {
        double sum=0.0;
        for(int i=0;i<a.length;i++)
        {
            sum+=a[i]*b[i];
        }
        return sum;
    }


    static double []mnoz(double a[],double x)//mnozenie wektora przez skalar
    {
        double r[]= new double[a.length];
        for(int i=0;i<a.length;i++)
        {
            r[i]=x*a[i];
        }
        return r;
    }

    static double []dod(double a[],double b[],int w)//dodawanie i odejmowanie wektora
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


    static double[] AX(double[][]a,double[]x)//mnozenie macierzy przez wektor
    {
        int rows =a.length;
        int columns = a[0].length;

        double[] result = new double[rows];

        for (int row = 0; row < rows; row++) {
            double sum = 0;
            for (int column = 0; column < columns; column++) {
                sum += a[row][column]
                        * x[column];
            }
            result[row] = sum;
        }
        return result;
    }


    static double[] wart_wlasna(double[][]A,double[]K1,double[]K2)
    {
        double E=1;
        double[]wyn=new double[2];

        double []y1=new double[A.length];
        double []z1=new double[A.length];
        double []z2=new double[A.length];
        double []y2=new double[A.length];
        double []k1=K1;
        double []k2=K2;

        for(int i=0;i<A.length;i++) {
            z2[i]=0;
            z1[i] = 0;
        }
        while(E>1E-15)//maksymalna dokladnosc dla ktorej program znajduje wartosc wlasna ( dla większej dokładnosci program dziala w nieskonczonosc)
        {
            //obliczanie wartosci 1
            z1 = AX(A, k1);//mnozenie macierzy A i K
            y1=mnoz(z1,(1/Math.sqrt(skal(z1,z1))));//obliczanie Y[k+1]
            E=Math.sqrt(skal(dod(y1,k1,0),dod(y1,k1,0)));//norma wektora(y[k+1]-y[k])
            System.arraycopy(y1,0,k1,0,y1.length);

        }
        E=1.0;
        //obliczanie drugiego wektora orgtogonalnego do pierwszego
        double sum=0;
        for(int g=0;g<K2.length-1;g++)
            sum+=z1[g];
        sum=-sum;
        k2[K2.length-1]=sum/z1[K2.length-1];

        //normalizacja wektora ortogonalnego
        double norma=Math.sqrt(skal(k2,k2));
        k2=mnoz(k2,1/norma);


        while(E>1E-15)//maksymalna dokladnosc dla ktorej program znajduje wartosc wlasna ( dla większej dokładnosci program dziala w nieskonczonosc) - zwiazane z dokladnoscia double w java
        {
            //obliczanie wartosci 2

            z2 = AX(A, k2);//mnozenie macierzy A i K

            z2=dod(z2,(mnoz(y1,skal(y1,z2))),0);

            y2=mnoz(z2,(1/Math.sqrt(skal(z2,z2))));//obliczanie Y[k+1]

            E=Math.sqrt(skal(dod(y2,k2,0),dod(y2,k2,0)));//norma wektora(y[k+1]-y[k])
            System.arraycopy(y2,0,k2,0,y2.length);
        }

        wyn[0]=Math.sqrt(skal(z1,z1));
        wyn[1]=Math.sqrt(skal(z2,z2));
        System.out.println("Wektor wlasny 1:");
        for(int j=0;j<z1.length;j++) {
            y1[j] = y1[j] / y1[z1.length - 1];
            System.out.print(y1[j]+" ");
        }
        System.out.println("\nWektor wlasny 2:");
        for(int j=0;j<z1.length;j++) {
            y2[j]=y2[j]/y2[z1.length-1];
            System.out.print(y2[j] + " ");
        }

        //roznica normy A*e[1]-lambda1*e[1]
        double[]g=dod(AX(A,y1),mnoz(y1,wyn[0]),0);
       // System.out.println("\nroznica normy A*e[1]-lambda1*e[1] "+Math.sqrt(skal(g,g)));

        //roznica normy A*e[2]-lambda2*e[2]
        g=dod(AX(A,y2),mnoz(y2,wyn[1]),0);
        //System.out.println("roznica normy A*e[2]-lambda2*e[2] "+Math.sqrt(skal(g,g)));


        return wyn;
    }
}
