//Autor : Piotr Piechowicz
import java.lang.Math;


public class Lagrange {
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
        double[][] mac=new double[siz][siz];
        for(int i=0;i<siz;i++)//wypełnianie macierzy
        {
            for (int j=0;j<siz;j++)
            {
                mac[i][j]=Math.pow(x[i], j);
            }
        }
        roz=luSolve(mac, fx,siz);//rozwiazanie rownania mac*x=fx faktoryzacja LU
        for(int i=0;i<siz;i++)
        {
            System.out.println(roz[i]);
        }
    }



    static double[] luSolve(double [][]mat,double []fx, int n) {
        double[][] lower = new double[n][n];
        double[][] upper = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int k = i; k < n; k++)// gorna macierz trojkatna
            {
                // sumowanie L(i, j) * U(j, k)
                double sum = 0.0;
                for (int j = 0; j < i; j++)
                    sum += (lower[i][j] * upper[j][k]);
                //obliczanie U(i, k)
                upper[i][k] = mat[i][k] - sum;
            }

            for (int k = i; k < n; k++)// dolna macierz trojkatna
            {
                if (i == k)
                    lower[i][i] = 1; // diagonala jako 1
                else {

                    // sumowanie L(k, j) * U(j, i)
                    double sum = 0.0;
                    for (int j = 0; j < i; j++)
                        sum += (lower[k][j] * upper[j][i]);

                    // obliczanie L(k, i)
                    lower[k][i] = (mat[k][i] - sum) / upper[i][i];
                }
            }
        }
        int N = mat.length;
        double[] Y = new double[n];
        double[] X = new double[n];
        int k = 0;
        double suml = 0;
        Y[0] = fx[0];
        for (int i = 1; i < n; i++) {//forwardsubstitution
            k = 0;
            suml = 0;
            while (k < i) {
                suml += lower[i][k] * Y[k];
                k++;
            }
            Y[i] = fx[i] - suml;
        }

        for (int i = n - 1; i >= 0; i--) {//backsubstitution
            suml = 0.0;
            for (int j = i + 1; j < n; j++) {
                suml += upper[i][j] * X[j];
            }
            X[i] = (Y[i] - suml) / upper[i][i];
        }
        return X;

    }

}
