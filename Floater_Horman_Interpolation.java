//Autor : Piotr Piechowicz

public class Horman {
    public static void main(String arr[])
    {
        int siz=65;
        double[] x=new double[siz];
        double[] fx = new double[siz];
        double[] roz;
        double st=-1.0;
        for(int i=0;i<siz;i++)
        {
            x[i]=st+i*1.0/32.0;
            fx[i]=1/(1+5*Math.pow(x[i],2));
        }

        floater(x,fx,3);
        
    }
    public static void floater(double x[],double fx[],int d)
    {
        int k=x.length;
        double [] wagi=new double[k];
        double suma;
        for(int i=0;i<k;i++)//poniewaz wezly roownoodlegle mozna skorzystac z uproszczonego wzoru na wagi
        {
            suma=0;
            for(int j=i-d;j<i;j++)
            {
                suma+=Newton(d,i-j);
            }
            int a=i-d;
            if(a<0)
                a=(-1*a);

            if((a)%2==1)
                suma=suma*(-1);
			
            wagi[i]=suma;
            System.out.println(i+" "+wagi[i]);
        }
        System.out.print("(");
        for(int z=0;z<k;z++)
        {
            if((z==k-1)) {

             if(x[z]>=0)
                System.out.print("(" + (wagi[z] * fx[z]) + "/(x-" + x[z] + ")) ");
             else
                 System.out.print("("+(wagi[z]*fx[z])+"/(x+"+(-x[z])+")) ");
            }
            else
            {
                if(x[z]>=0)
                    System.out.print("(" + (wagi[z] * fx[z]) + "/(x-" + x[z] + ")) +");
                else
                    System.out.print("("+(wagi[z]*fx[z])+"/(x+"+(-x[z])+")) +");
            }

        }
        System.out.print(")/\n(");
        for(int w=0;w<k;w++)
        {
            if(w==k-1) {
                if (x[w] >= 0)
                    System.out.print("(" + (wagi[w]) + "/(x-" + x[w] + ")) ");

                else
                    System.out.print("(" + (wagi[w]) + "/(x+"+(-x[w]) + ")) ");
            }
            else
            {
                if (x[w] >= 0)
                    System.out.print("(" + (wagi[w]) + "/(x-" + x[w] + ")) + ");

                else
                    System.out.print("(" + (wagi[w]) + "/(x+"+(-x[w]) + ")) + ");
            }
        }
        System.out.print(")");
    }

    public static long Newton( int n, int k )
    {
        long  Wynik = 1;       // Deklaracja zmiennych
        int i;

        for(i = 1; i <= k; i++) // Od 1 do k wykonujemy :
        {
            Wynik = Wynik * ( n - i + 1 ) / i;      // Obliczanie ze wzoru iteracyjnego
        }

        return Wynik;   // Zwróć Wynik
    }

}
