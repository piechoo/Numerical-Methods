//Autor: Piotr Piechowicz
import org.apache.commons.math3.complex.Complex;
public class Laguerre_Method{
    public static void main(String[] args) {

        double []x1 = { 16, -72, -28, 558, -990, 783, -486, 243  };//wspolczynniki pierwszego wielomianu
        double []x2 = { -4, -4, -12, -8, -11,-3, -1, 2, 3, 1, 1  };//wspolczynniki drugiego wielomianu

        Complex f_zesp[]=new Complex[5];//wspolczynniki 3 wielomianu
        f_zesp[0]=new Complex(1,0);
        f_zesp[1]=new Complex(0,-1);
        f_zesp[2]=new Complex(-1,0);
        f_zesp[3]=new Complex(0,1);
        f_zesp[4]=new Complex(1,0);
        Complex fcja []=new Complex[8];
        Complex fcja2[]=new Complex[11];
        for(int i=0;i<8;i++)
            fcja[i]=new Complex(x1[i]);
        Polynomial p1=new Polynomial(fcja);
        for(int i=0;i<11;i++)
            fcja2[i]=new Complex(x2[i]);

        Polynomial p2=new Polynomial(fcja2);

        Polynomial p3=new Polynomial(f_zesp);
        System.out.println("f1: 243z^7 − 486z^6 + 783z^5 − 990z^4 + 558z^3 − 28z^2 − 72z + 16 ");
        p1.findRoots();
        System.out.println();

        System.out.println("f2: z^10 + z^9 + 3z^8 + 2z^7 − z^6 − 3z^5 − 11z^4 − 8z^3 − 12z^2 − 4z − 4 ");
        p2.findRoots();
        System.out.println();
        System.out.println("f3: z^4 + iz^3 − z^2 − iz + 1 ");
        p3.findRoots();
    }


}
class Polynomial {
    // współczynniki danego wielomianu P
    private Complex[] polynomialCoeff;


    public Polynomial(Complex[] polynomialCoeff)
    {
        this.polynomialCoeff = polynomialCoeff;
    }

    // funkcja startująca algorytm
    public Complex[] findRoots(){
        Complex start = new Complex(0,0);//poczatkowe przyblizenie
        Complex[] z = new Complex[polynomialCoeff.length -1];//wektor miejsc zerowych
        int i = 0;

        Complex[] polLess = new Complex[polynomialCoeff.length];
        z[0] = LaguerreMethod(start, polynomialCoeff);

        Complex[] tmp = polynomialCoeff;

        for(i = 1; i < polynomialCoeff.length -1; i++)
        {
            polLess = newPol(z[i-1],tmp);//deflacja
            z[i] = LaguerreMethod(start, polLess);//oblienie miejsca zerowego
            z[i] = LaguerreMethod(z[i], polynomialCoeff);//wygladzanie miejsca zerowego
            tmp = polLess;
        }


        for(i = 0; i < polynomialCoeff.length -1; i++){
            System.out.println(z[i]);//wypisanie rozwiazan
        }

        return z;
    }

    // Metoda Laguerre'a
//  argumenty:
//  start - punkt początkowy
//  pol - tablica przechowująca współczynniki wielomianu, którego szukamy mz.
//  zmienne:
//  fun - P(z)*n
//  der - P'(z), secDer - P''(z)
//  denominator - mianownik
//  zwraca:
//  z - miejsce zerowe
    private Complex LaguerreMethod(Complex start, Complex[] pol){
        Complex z = start;
        Complex fun,der,secDer, denominator;
        Complex tmp = new Complex(100,0);

// iterujemy tak długo aż różnica między aktualnym a poprzednio
// wyliczonym miejscem zerowym będzie nie mniejsca niż
// zadana precyzja(w tym przypadku 1e-10)

        while((horner(z, pol).subtract(horner(tmp, pol))).abs() > 1e-10){
            tmp  = z;
            fun =    horner(z, pol).multiply(pol.length-1);
            der = horner(z, derivative(pol));
            secDer = horner(z, derivative(derivative(pol)));
            denominator = ( (  der.multiply(der).multiply(((pol.length-1) - 1))
                    .subtract(fun.multiply(secDer)) ).multiply((pol.length-1) - 1) ).sqrt();
            denominator = der.subtract(denominator).abs() > der.add(denominator).abs()
                    ? der.subtract(denominator) : der.add(denominator);
            z = z.subtract((fun.divide(denominator)));
        }
        return z;
    }

    // wartość wielomianu w z wyliczona za pomocą Hornera
    private Complex horner(Complex z, Complex[] coefficients){
        Complex P = coefficients[coefficients.length-1];
        int k = coefficients.length-1;
        while(k >0){
            k--;
            P = P.multiply(z).add(coefficients[k]);
        }
        return P;
    }


    // deflacja
    public Complex[] newPol(Complex root, Complex[] oldPol){
        Complex[] newPol = new Complex[oldPol.length-1];
        newPol[oldPol.length -2] = oldPol[oldPol.length -1];
        for(int i = oldPol.length - 3; i >= 0; i--){
            newPol[i] = oldPol[i+1].add(newPol[i+1].multiply(root));
        }
        return newPol;
    }
    // funkcja licząca pochodną
    public Complex[] derivative(Complex[] pol){
        if(pol.length-1 == 0) return new Complex[]{new Complex(0,0)};
        Complex[] der = new Complex[pol.length-1];

        for(int i = pol.length -2; i >= 0; i-- ){
            der[i] = pol[i+1].multiply(i+1);
        }

        return der;
    }
}
