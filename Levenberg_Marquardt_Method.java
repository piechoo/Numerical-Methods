//Author : Piotr Piechowicz
import java.util.Arrays;

public class  Levenberg_Marquardt_Method {

    public static void main(String[] args) {
        double[]x={1.036,3.11351};//początkowe przybliżenie
        double lam=1.0/1024.0;//współczynnik lambda
        double []roz=leven(x,lam);
        System.out.println("x "+roz[0]);
        System.out.println("y "+roz[1]);
    }
static double []leven(double[]y,double lam) {//algorytm levenberga-Marquardta
    double[]x=new double[y.length];
    System.arraycopy(y,0,x,0,y.length);
    double lamax = 1e16;
    double lamin = 1e-16;
    double[] g = grad(x);
    int licznik=0;
    double[][] h = hes(x, lam);
    int wsk=0;
    double[] temp = dod(x, luSolve(h, g, g.length), 0);
    System.arraycopy(temp, 0, x, 0, x.length);
    double[]k=luSolve(h,g,g.length);
    while (true)
    {
        licznik++;
        System.out.println(licznik);
        if(wsk>0){
            g=grad(x);
            wsk--;
        }
        h = hes(x, lam);
        temp = dod(x, luSolve(h, g, g.length), 0);
        k=luSolve(h,g,g.length);
        if (fcja(temp) == fcja(x)){
            return temp;
        }
        if (fcja(temp) >= fcja(x)) {
            lam *= 8.0;
        }
        if (fcja(temp) < fcja(x)) {
            lam *= 1.0 / 8.0;
            System.arraycopy(temp, 0, x, 0, x.length);
            wsk++;
        }
        if (lam < lamin) {//obliczanie minimum metoda sprzeżonych gradienów
            x=sprz(x,grad(x),licznik);
            return x;


        } else if (lam > lamax) {
            System.out.println("Punkt poza basenem atrakcji minimum! Wybierz inny punkt ");
            return x;

        }
        //System.out.println("Graphics3D[{Red, PointSize[0.02], Point[{"+x[0]+", "+x[1]+", "+fcja(x)+"}]}],");
        System.out.println("x1 "+x[0]);
        System.out.println("x2 "+x[1]+"\n");
    }

}


    static double[] grad(double[]x)//wartosc gradientu funkcji dla danych zmiennych
    {
        double []g=new double[2];
        g[0]=400.0*x[0]*x[0]*x[0]-400.0*x[0]*x[1]+2.0*x[0]-2.0;
        g[1]=-200.0*x[0]*x[0]+200.0*x[1];
        return g;
    }
    static double [][]hes(double[]x,double lam)//hesjan potrzebny do metody levenberga-marquardta
    {
        double[][]h=new double[2][2];
        h[0][0]=(1+lam)*(1200*x[0]*x[0]-400*x[1]+2.0);
        h[1][1]=(1+lam)*200.0;
        h[1][0]=(-400*x[0]);
        h[0][1]=h[1][0];
        return h;
    }

   static double fcja(double[]x)//obliczanie wartosci w danym punkcie
    {
        return 100*Math.pow(x[0],4)-200*Math.pow(x[0],2)*x[1]+Math.pow(x[0],2)-2*x[0]+100*Math.pow(x[1],2)+1;
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

    static double []mnoz(double a[],double x)//mnozenie wektora przez skalar
    {
        double r[]= new double[a.length];
        for(int i=0;i<a.length;i++)
        {
            r[i]=x*a[i];
        }
        return r;
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


   static double[] luSolve(double [][]mat,double []fx, int n) {//rozwiazywanie ukladu rownan faktoryzacja Lu
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
   static double f_kierun(double[]x,double[]p,double a)//obliczanie wartosci funkcji w danym kierunku i z danym parametrem a
    {
        double q= 100*(Math.pow(x[0],4)+4*Math.pow(x[0],3)*a*p[0]+6*x[0]*x[0]*a*a*p[0]*p[0]+4*x[0]*Math.pow(a*p[0],3)+Math.pow(a*p[0],4));
        double w= 200*(x[0]*x[0]+2*x[0]*p[0]*a+Math.pow(a*p[0],2))*(x[1]+a*p[1]);
        double e= Math.pow((x[0]+a*p[0]),2);
        double r=(-2*(x[0]+a*p[0]));
        double t=100*(Math.pow((x[1]+a*p[1]),2));
        return q-w+e+r+t+1.0;
    }

    static double[] abc(double[]x,double[]p,double a)//znalezenie 3 punktow spełniajacych a<b<c i f(a)>f(b) i f(c)>f(b) , poczatkowy punkt a
    {
        double[] abc=new double[3];
        double[] wynik=new double[3];
        int k=0;
        double e=0.1;
        wynik[0]=f_kierun(x,p,a);
        abc[0]=a;
        wynik[1]=f_kierun(x,p,a+e);
        abc[1]=a+e;
        if(wynik[1]>wynik[0])
        {
            e=(-e);
            abc[2]=abc[1];
            abc[1]=abc[0];
            while(f_kierun(x,p,abc[1]+e)<=wynik[0])
            {
                e=(2*e);
                k++;
                if(k>1e5)
                {
                    System.out.println("Monotoniczna galaz !");
                    break;
                }
            }
            abc[0]=abc[1]+e;
        }
        else
        {
            while(f_kierun(x,p,a+e)<=wynik[0])
            {
                e=(2*e);
                k++;
                if(k>1e5)
                {
                    System.out.println("Monotoniczna galaz !");
                    break;
                }
            }
            abc[0]=a+e;
        }
        return abc;
    }

    static double brent(double[] a,double[] x,double []p )//funkcja obliczajaca minimum kierunkowe funkcji w kierunku p metodą brenta
    {
        double []ab=new double[3];
        double d=0.0;
        double []dlugosc=new double[3];
        double[]temp=new double[3];
        dlugosc [0]=a[2]-a[0];
        ab=bisekcja(x,p,a,d);
        dlugosc[1]=ab[2]-ab[0];
        double b=0.0;
        double tmp;
        int i=2;
        dlugosc[i]=1;
        int k=0;
        while(Math.abs(ab[2]-ab[0])>1e-15*(Math.abs(b)+Math.abs(d)))
        {
            b=ab[1];
            d=d_bren(ab,p,x);
            if(d==420.69)
                break;
            temp = warun(ab, d, x, p);
            tmp = temp[2] - temp[0];
            if (((ab[0] < d) && (ab[2] > d)) && (tmp < (dlugosc[k] / 2)))
            {
                System.arraycopy(temp, 0, ab, 0, ab.length);
                dlugosc[i]=tmp;
            }
            else {
                temp = bisekcja(x, p, ab,d);
                dlugosc[i] = temp[2] - temp[0];
                System.arraycopy(temp, 0, ab, 0, ab.length);
            }
            i=(i+1)%3;
            k=(k+1)%3;
        }
        return ab[1];
    }

    static double[] warun(double[]abc,double d,double[]x,double[]p)//warunki a<b<c i f(a)>f(b) i f(c)>f(b)
    {
        double[]aa=new double[3];
        if(f_kierun(x,p,d)<f_kierun(x,p,abc[1]))
        {
            if(d<abc[1])
            {
                aa[2]=abc[1];
                aa[1]=d;
                aa[0]=abc[0];
            }
            else
            {
                aa[0]=abc[1];
                aa[1]=d;
                aa[2]=abc[2];
            }
        }
        else
        {
            if(d<abc[1])
            {
                aa[0]=d;
                aa[1]=abc[1];
                aa[2]=abc[2];
            }
            else
            {
                aa[2] = d;
                aa[1]=abc[1];
                aa[0]=abc[0];
            }
        }
        return aa;
    }

    static double[] bisekcja(double[]x,double[]p,double[]abc,double q)
    {
        double d=abc[0]+(abc[2]-abc[0])/2;
        if(d==abc[1])
            d=d+(abc[2]-abc[0])/100.0;
        q=d;
        double []ab=warun(abc,d,x,p);
        return ab;
    }

    static double []sprz(double[]xk,double[] gr,int licz)//metoda gradientow sprzeżonych
    {
        double[]x=new double[xk.length];
        double bet;
        System.arraycopy(xk,0,x,0,x.length);
        double []r1=mnoz(gr,-1);
        double []p=mnoz(gr,-1);
        double []r2;
        double s=brent(abc(x,p,0.0),x,p);
        double []x2=dod(x, mnoz(p, s), 2);

        r2 = mnoz(grad(x2), -1);
        bet = skal(r2, dod(r2,r1,0)) / skal(r1, r1);
        p = dod(r2, mnoz(p, bet), 2);

        while(Math.abs(Math.sqrt(skal(x2,x2))-Math.sqrt(skal(x,x)))>1e-16)
        {
            licz++;
            System.out.println(licz+"\n");
            System.out.println("x "+x2[0]);
            System.out.println("y "+x2[1]);
            System.arraycopy(x2,0,x,0,x.length);
            s=brent(abc(x2,p,0.0),x2,p);
            x2 = dod(x, mnoz(p, s), 2);
            System.arraycopy(r2,0,r1,0,r1.length);
            r2 = mnoz(grad(x2), -1);
            bet = skal(r2, dod(r2,r1,0)) / skal(r1, r1);
            p = dod(r2, mnoz(p, bet), 2);
        }

        return x2;
    }

    static double d_bren(double[]ab,double[]p,double[]x)//obliczanie minimum paraboli w metodzie brenta
    {
        double d;
        double gora=Math.pow(ab[0],2)*(f_kierun(x, p, ab[2]) - f_kierun(x, p, ab[1])) + Math.pow(ab[1],2)* (f_kierun(x, p, ab[0]) - f_kierun(x, p, ab[2])) + Math.pow(ab[2],2)*(f_kierun(x, p, ab[1]) - f_kierun(x, p, ab[0]));
        double dol=ab[0]*(f_kierun(x, p, ab[2]) - f_kierun(x, p, ab[1])) + ab[1]*(f_kierun(x, p, ab[0]) - f_kierun(x, p, ab[2])) + ab[2]*(f_kierun(x, p, ab[1]) - f_kierun(x, p, ab[0]));
        if(gora==0&&dol==0)
            return 420.69;
        else if(ab[1]>1e16||ab[2]>1e16||ab[0]>1e16)
            return 1e26;
        return 0.5*gora/dol;
    }
}
