
//Autor : Piotr Piechowicz
public class Eigenvectors {

    public static void main(String[] args) {
        double lamb=0.38197;
        double[][] A = { { 2-lamb, -1, 0, 0, 1 },//macierz A-Lambda*1
                { -1, 2-lamb, 1, 0, 0 },
                { 0, 1, 1-lamb, 1, 0 },
                { 0, 0, 1, 2-lamb, -1 },
                { 1, 0, 0, -1, 2-lamb },
        };
        double[] y=new double[A.length];//wektor unormowany od ktorego zaczynam przyblizanie wektora wlasnego
        y[0]=1;

        double[] e=wektor_wlasny(A,y);
        //System.out.println(Math.sqrt(skal(e,e)));//sprawdzenie czy wektor unoemowany
        System.out.println("Unormowany wektor wlasny dla lambda = "+lamb+" :");
        for (int i=0;i<y.length;i++)
            System.out.print((e[i])+" ");
        System.out.println();

    }



    static double[]wektor_wlasny(double[][]A,double[]y)
    {
        double[] k=new double[y.length];
        System.arraycopy(y,0,k,0,y.length);
        double[] Y1=new double[y.length];
        System.arraycopy(y,0,Y1,0,y.length);
        double[] e=new double[A.length];
        double E=1;
        while(E>1E-10)
        {
            //obliczanie wartosci unormowanego wektora wlasnego
            e = luSolve(A,Y1,Y1.length);//obliczanie z metoda lu

            Y1=mnoz(e,(1/Math.sqrt(skal(e,e))));//obliczanie Y[k+1]
            E=Math.sqrt(skal(dod(Y1,k,0),dod(Y1,k,0)));//norma wektora(y[k+1]-y[k])
            if(warunek(Y1,k))//jesli kolejne wektory oscyluja: (y[k+1]=-y[k]) to przerwij
                break;
            System.arraycopy(Y1,0,k,0,y.length);
        }
        return Y1;
    }

    static boolean warunek(double[]a,double[]b)//jesli kolejne wektory oscyluja: (y[k+1]=-y[k]) to przerwij
    {
        int k=0;
        for(int i=0;i<a.length;i++)
            if(a[i]==(-b[i]))
                k++;

            if(k==a.length)
                return true;
            else
                return false;
    }

    static double[] luSolve(double [][]mat,double []fx, int n) {//rozwiazywanie macierzy metoda lu
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
}
