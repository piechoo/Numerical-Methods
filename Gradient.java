//Autor:Piotr Piechowicz
public class Gradient {
    public double a[][] = new double[128][128];
    public double e[]=new double[128];
    public double p[]=new double[128];
    public double r[]=new double[128];
    public double x[]=new double[128];
    public Gradient(){
        for(int i=0; i<128;i++){
            e[i]=1.0;
            x[i]=0.0;
            a[i][i]=4.0;
            if(i>3)
                a[i][i-4]=1.0;
            if(i>0)
                a[i][i-1]=1.0;
            if(i<124)
                a[i][i+4]=1.0;
            if(i<127)
                a[i][i+1]=1.0;
        }
    }
    double [] grad(){
        double alfa;
        double beta;
        double ct=1;
        double roz[]=new double [r.length];
        double roz1[]= new double[r.length];
        System.arraycopy(p, 0, x, 0, p.length);
        double ap[]=new double[p.length];
        ap=mac(a,x);
        int k=0;
        r=dod(e,ap,0);//obliczenie r wstepnego = b-Ax gdzie x to wstepne rozwiazanie
        System.arraycopy(r, 0, p, 0, r.length);
        double bp[]=new double[p.length];
        double r1[]=new double [p.length];
        //while(skal(r,r)>1E-170)
        while(ct>0)//wykonuje dopoki ||x(k+1)-x(k)|| nie rowna sie 0
        {
  
            alfa = skal(r, r)/tran(p,a);//obliczanie alfa
            ap=mac(a,p);
            k++;
            ap=mnoz(ap,alfa);
			
            r1=dod(r,ap,0);//obliczanie wektora r
            
			beta=skal(r1,r1)/skal(r,r);//obliczanie bety
            
			System.arraycopy(r1, 0, r, 0, r.length);
            bp=mnoz(p,beta);
            ap=mnoz(p,alfa);
            p=dod(r1,bp,2);
            
			roz1=dod(roz,ap,2);//obliczenie rozwiazania x(k+1)
            roz1=dod(roz1,roz,0);//roznica x(k+1)-x(k)
            ct=skal(roz1,roz1);//iloczyn skalarny x(k+1)-x(k)
            
			System.out.println(ct);
            
			roz1=dod(roz,ap,2);//obliczenie rozwiazania x(k+1)
            System.arraycopy(roz1, 0, roz, 0, r.length);
        }
        return roz;
    }
    double skal(double a[],double b[])//obliczanie iloczynu skalarnego wektorów
    {
        double sum=0.0;
        for(int i=0;i<a.length;i++)
        {
            sum+=a[i]*b[i];
        }
        return sum;
    }
    double tran(double p[],double a[][])//wektor t razy macierz razy wektor
    {
        double[] il;
        il=mac(a,p);
        return skal(il,p);
    }
    double[] mac(double a[][],double p[])//macierz razy wektro
    {
        double e[] = new double [p.length];
        e[0]=a[0][1]*p[1]+a[0][4]*p[4]+a[0][0]*p[0];
        e[1]=a[1][0]*p[0]+a[1][2]*p[2]+a[1][5]*p[5]+a[1][1]*p[1];
        e[2]=a[2][1]*p[1]+a[2][3]*p[3]+a[2][6]*p[6]+a[2][2]*p[2];
        e[3]=a[3][2]*p[2]+a[3][4]*p[4]+a[3][7]*p[6]+a[3][3]*p[3];
        for(int j=4;j<124;j++){
            e[j]=a[j][j-4]*p[j-4]+a[j][j-1]*p[j-1]+a[j][j+1]*p[j+1]+a[j][j+4]*p[j+4]+a[j][j]*p[j];
        }
        e[124]=a[124][120]*p[120]+a[124][123]*p[123]+a[124][125]*p[125]+a[124][124]*p[124];
        e[125]=a[125][121]*p[121]+a[125][124]*p[124]+a[125][126]*p[126]+a[125][125]*p[125];
        e[126]=a[126][122]*p[122]+a[126][125]*p[125]+a[126][127]*p[127]+a[126][126]*p[126];
        e[127]=a[127][123]*p[123]+a[127][126]*p[126]+a[127][127]*p[127];
        return e;
    }
    double []mnoz(double a[],double x)//mnozenie wektora przez skalar
    {
        double r[]= new double[a.length];
        for(int i=0;i<a.length;i++)
        {
            r[i]=x*a[i];
        }
        return r;
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
    public static void main(String[] args)
    {
        double roz[]=new double [128];
        Gradient g1= new Gradient();
        roz=g1.grad();
        //wypisanie rozwiazan
        for(int i=0;i<128;i++){
            System.out.println("x nr:" + i +" "+ roz[i]);
        }

        //sprawdzenie
        for(int j=4;j<124;j++){
            System.out.println((4*roz[j]+roz[j-4]+roz[j-1]+roz[j+1]+roz[j+4]));
        }
        //System.out.println("Wykonano "+g1.cont+" krokow");

    }
}
