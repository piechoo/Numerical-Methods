//Autor : Piotr Piechowicz
import java.lang.*;
import java.util.*;

import static java.lang.Math.abs;

public class Romberg_Integral {

    static double wart(double x)//obliczanie wartosci funkcji podcalkowej
    {
        return Math.sin(Math.PI*(1+Math.sqrt(x))/(1+x*x))*Math.exp(-x);
    }

    static double trapez(double a,double b,int k)//funkcja obliczajaca calke metoda trapezow( a poczatek przedzialu, b koniec, 2^k ilosc przedzialow)
    {
        double przedz=Math.pow(2,k);
        double sum=0;
        double h=(b-a)/przedz;

        for(int i=0;i<przedz+1;i++)
        {
            if((i==0)||(i==przedz))
                sum+=wart(a+(i/przedz)*(b-a))/2;
            else {
                sum += wart(a + (i / przedz) * (b - a));
            }
        }
        return h*sum;
    }

    static double romberg(double a,double b)
    {
        int z=0;//numer wiersza w tablicy
        int n=0;//numer kolumny
        double temp;//zmienna tymczasowa
        double temp_next;//zmienna tymczasowa
        double comp=0;//zmienna w ktorej zapisywany poprzedni diagonalny element (do porownania)

        ArrayList<Double> arr = new ArrayList();
        //pierwsze dwa wiersze obliczone przed petla aby latwiej napisac warunek
        arr.add(trapez(a,b,z));
        n++;
        z++;
        arr.add((1/(Math.pow(4,n)-1))*(Math.pow(4,n)*trapez(a,b,z)-arr.get(0)));
        z++;
        arr.set(0,trapez(a,b,1));
        n=0;
        while(abs(arr.get(arr.size()-1)-comp)>1E-7)
        {
                temp=trapez(a,b,z);
                while(n<arr.size())
                {
                    temp_next=(1/(Math.pow(4,(n+1))-1)*(Math.pow(4,(n+1))*temp-arr.get(n)));
                    comp=arr.get(n);
                    arr.set(n,temp);
                    temp=temp_next;
                    n++;
                }
                //System.out.println(arr.get(0));
                arr.add(temp);
                n=0;
                z++;
        }
        return arr.get(arr.size()-1);
    }

    public static void main(String[] args)
    {
        System.out.println("Wynik: "+(romberg(0,17)+Math.exp(-17)));//calka od 0 do 17 plus przyblizenie calki od 17 do +inf(tak naprawde nie jest konieczne bo dokladnosc <10-7)
    }
}
