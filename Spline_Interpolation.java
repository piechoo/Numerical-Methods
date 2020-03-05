//Autor: Piotr Piechowicz
import java.lang.Math;

public class Splajn {

    public static void main(String arr[])
    {
        int siz=65;
        double[] x=new double[siz];
        double[] fx = new double[siz];
        double[] dzet=new double[siz-2];
        double[][] troj_dzet=new double[siz-2][siz-2];
        double[] roz_dzet=new double[siz-2];
        double[] dzety=new double[siz];
        double st=-1.0;
        double[][] funkcje=new double[siz-1][4];

        for(int i=0;i<siz;i++)
        {
            x[i]=st+(i*1.0/32.0);//wypełnienie wartosci x
            fx[i]=1/(1+5*Math.pow(x[i],2));//wypełnienie wartosci fx

            i=i-2;//przesunięcie żeby dane potrzebne do obliczenia ksi byly wypełnione
            if((i<siz-3)&&(i>0))//wypelnianie wartosci macierzy trojdiagonalnej do znalezienia ksi( x równoodległe wiec wystarczy diagonala z 4 i dwa pasma z 1)
            {
                troj_dzet[i][i]=4.0;
                troj_dzet[i][i-1]=1.0;
                troj_dzet[i][i+1]=1.0;
                roz_dzet[i]=6.0*(fx[i]-2*fx[i+1]+fx[i+2])/Math.pow(0.03125,2);
            }
            i++;
            i++;
        }
        troj_dzet[0][0]=4.0;
        roz_dzet[0]=6.0/Math.pow(0.03125,2)*(fx[0]-2*fx[1]+fx[2]);
        roz_dzet[siz-3]=6.0/Math.pow(0.03125,2)*(fx[siz-3]-2*fx[siz-2]+fx[siz-1]);
        troj_dzet[0][1]=1.0;
        troj_dzet[siz-3][siz-3]=4.0;
        troj_dzet[siz-3][siz-4]=1.0;

        double[][] mac=new double[siz][siz];
        for(int i=0;i<siz;i++)//wypełnianie macierzy wartosciami x i f(x)
        {
            for (int j=0;j<siz;j++)
            {
                mac[i][j]=Math.pow(x[i], j);
            }
        }
        dzet=thomas(troj_dzet,roz_dzet);
        dzety[0]=0;//ksi[0]=0;
        dzety[siz-1]=0;//ksi[n]=0;

        for(int i=0;i<siz-2;i++)
        {
            dzety[i+1]=dzet[i];
        }

        funkcje=fcje(x,dzety,fx);

        for(int i=0;i<siz-1;i++)
        {
           System.out.println(i+".    "+x[i]+","+x[i+1]);
           System.out.println(i+". "+funkcje[i][3]+"*x^3+"+funkcje[i][2]+"*x^2+"+funkcje[i][1]+"*x+"+funkcje[i][0]);
        }


    }


    public static double[] thomas(double[][]A,double[]D) {
        int N = A.length;
        double[] bet=new double[N];
        double[] gam=new double[N];
        double[] X=new double[N];

        bet[0]=-A[0][1]/A[0][0];
        gam[0]=D[0]/A[0][0];
        for(int i=1;i<N;i++){
            if(i<N-1) {
                bet[i] = -(A[i][i + 1]) / (A[i][i - 1] * bet[i - 1] + A[i][i]);
                gam[i] = (D[i] - A[i][i - 1] * gam[i - 1]) / (A[i][i - 1] * bet[i - 1] + A[i][i]);
            }
            else {//obsluga przypadku gdy pobieramy wartosc "z zza tablicy" wtedy c(n) = 0
                bet[i] = 0.0 ;
                gam[i] = (D[i] - A[i][i - 1] * gam[i - 1]) / (A[i][i - 1] * bet[i - 1] + A[i][i]);
            }

        }

        X[N-1]=gam[N-1];
        for(int j = N-2; j>-1;j--){
            X[j]=bet[j]*X[j+1]+gam[j];
        }
        return X;
    }

    public static double[][] fcje(double[] x,double[] dzet,double[] fx )
    {
        int k=fx.length;
        double [][]wyn=new double[k-1][4];
        double roznica=0.03125;
        for(int i=0;i<k-1;i++)
        {
            //C*ksi[i]
            wyn[i][3]=(-dzet[i])/(6.0*roznica);
            wyn[i][2]=(3*x[i+1]*dzet[i])/(6.0*roznica);
            wyn[i][1]=(-3*Math.pow(x[i+1],2)*dzet[i])/(6.0*roznica)+(dzet[i]*roznica)/6.0;
            wyn[i][0]=(Math.pow(x[i+1],3)*dzet[i])/(6.0*roznica)-(x[i+1]*roznica*dzet[i])/6.0;

            //D*ksi[i+1]
            wyn[i][3]+=(dzet[i+1])/(6.0*roznica);
            wyn[i][2]+=(-3*x[i]*dzet[i+1])/(6.0*roznica);
            wyn[i][1]+=(3*Math.pow(x[i],2)*dzet[i+1])/(6.0*roznica)-(dzet[i+1]*roznica)/6.0;
            wyn[i][0]+=(-Math.pow(x[i],3)*dzet[i+1])/(6.0*roznica)+(x[i]*roznica*dzet[i+1])/6.0;

            //A*f[i]
            wyn[i][1]+=(-fx[i])/roznica;
            wyn[i][0]+=(fx[i]*x[i+1])/roznica;

            //B*f[i+1]
            wyn[i][1]+=(fx[i+1])/roznica;
            wyn[i][0]+=(-fx[i+1]*x[i])/roznica;
        }


    return wyn;
    }
}
