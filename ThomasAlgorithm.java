
//Autor: Piotr Piechowicz
public class AlgThomas {//uzywam algorytmu thomasa ktory dla macierzy trojdiagonalnej ma zlozonosc o(n)

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



    public static double[] sherman(double[][] A,double[] B,double[] V) {
        int N  = A.length;
        double[] W = new double[N];
        double[] Z = thomas(A,B);
        double[] Q = thomas(A,V);
        for(int z=0;z<Z.length;z++)
        {
            System.out.println("Vector z: "+Z[z]);
        }
        for(int z=0;z<Z.length;z++)
        {
            System.out.println("Vector q: "+ Q[z]);
        }
        double gora=0;
        double dol=1;
        for(int i=0;i<N;i++){
            gora+=(V[i]*Z[i]);
            dol+=(V[i]*Q[i]);
        }
        double mnoz=gora/dol;
        for(int i=0;i<N;i++){
            W[i]=Z[i]-mnoz*Q[i];
        }
        return W;
    }

    public static void main(String[] args) {
        int N = 7;
        double[] B = { 1, 2, 3, 4, 5, 6, 7 };
        double[][] A = { { 4, 1, 0, 0, 0, 0, 0 },
                { 1, 4, 1, 0, 0, 0, 0 },
                { 0, 1, 4, 1, 0, 0, 0 },
                { 0, 0, 1, 4, 1, 0, 0 },
                { 0, 0, 0, 1, 4, 1, 0 },
                { 0, 0, 0, 0, 1, 4, 1 },
                { 0, 0, 0, 0, 0, 1, 4 },
        };
        double[][] A1 = { { 3, 1, 0, 0, 0, 0, 0 },
                { 1, 4, 1, 0, 0, 0, 0 },
                { 0, 1, 4, 1, 0, 0, 0 },
                { 0, 0, 1, 4, 1, 0, 0 },
                { 0, 0, 0, 1, 4, 1, 0 },
                { 0, 0, 0, 0, 1, 4, 1 },
                { 0, 0, 0, 0, 0, 1, 3 },
        };
        double[] V ={1,0,0,0,0,0,1};

        double[] R1 = thomas(A,B);

        double[] R2 = sherman(A1,B,V);
        for (int i = 0; i < N; i++) {
            System.out.println("X1 nr"+(i+1)+": "+ R1[i]);
        }
        System.out.println();
        for (int i = 0; i < N; i++) {
            System.out.println("X2 nr" + (i + 1) + ": " + R2[i]);
        }
    }

}