#include "TF1.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TFile.h"
#include "math.h"
#define A 2
#define B 1
#define C 2
#define D 2

double Myfunction(double x)
{
	return A*pow(x,4)*exp(-C*x*x) + 0.5*B*pow(x/2 + 3/2, 4)*exp(-D*pow(x/2 + 3/2,2));
}
void integrazione(int nInter=10)
{
	double a=1.;
	double b=10.;
	TF1* PmyF=new TF1("MyF","Myfunction(x)",a,b);		
	PmyF->Draw();
	double a1 = 1.;
	double b1 = 7.;
	
	//integrazione con la regola dei trapezi
	double h= b1 - a1;
	double IntTrapezi = (Myfunction(a1)+ Myfunction(b1))*h/2.;
	//integrazione con la regola di Simpson
	double h1 = h/2;
	double IntSimpson = h1*(Myfunction(a1)+4*Myfunction((a1+b1)/2)+Myfunction(b1))/3;
	cout<<" Integrale con la regola dei trapezi tra "<< a1 <<" e "<< b1<<" e' "<<IntTrapezi<<endl;
	cout<<" Integrale con la regola di Simpson tra "<< a1 <<" e "<< b1<<" e' "<<IntSimpson<<endl;
	//integrazione con la regola dei trapezi composita
	double IntTrapeziComposito = 0.;
	double H = (b1 - a1)/nInter;
	for(int i=1; i < nInter; i++)
	{
		double x = a1+H*i;
		IntTrapeziComposito = IntTrapeziComposito + Myfunction(x);
	}
	IntTrapeziComposito = IntTrapeziComposito + Myfunction(a1)/2 + Myfunction(b1)/2;
	IntTrapeziComposito = IntTrapeziComposito*H;
	cout<<" Integrale con la regola dei trapezi composita tra "<< a1 <<" e "<< b1<<" e' "<<IntTrapeziComposito<<endl;
	
	//integrazione con la regola di Simpson composita
	double IntSimpsonComposito = 0.;
	double H = (b1 - a1)/(2*nInter);
	for(int i=1; i < 2*nInter; i++)
	{
		double x = a1+H*i;
		if (i%2 == 0)
		{
	    	IntSimpsonComposito = IntSimpsonComposito + 2*Myfunction(x);
		}
		else
		{
			IntSimpsonComposito = IntSimpsonComposito + 4*Myfunction(x);
		}
	}
	IntSimpsonComposito = IntSimpsonComposito + Myfunction(a1) + Myfunction(b1);
	IntSimpsonComposito = IntSimpsonComposito*H/3;
	cout<<" Integrale con la regola di Simpson composita tra "<< a1 <<" e "<< b1<<" e' "<<IntSimpsonComposito<<endl;
}
