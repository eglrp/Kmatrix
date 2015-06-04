#ifndef JAC_EQS_H
#define JAC_EQS_H
#include <math.h>
#include <stdio.h>

template <typename T>
T Power(T base, T pow){return std::pow(base, pow);}

template <typename T>
T Sqrt(T exp){return std::sqrt(exp);}

// jacobian formulas for H^{-1}
template <typename T>
void jac_inv_equations(const vector<T> &P, const matrix<T> &s, T u, T v, T &dEdh1, T &dEdh2, T &dEdh3, T &dEdh4, T &dEdh5,
	T &dEdh6, T &dEdh7, T &dEdh8, T &dEdh9)
{
	T h1 = P[0], h2 = P[1], h3 = P[2];
	T h4 = P[3], h5 = P[4], h6 = P[5];
	T h7 = P[6], h8 = P[7], h9 = P[8];
	T a = s(0,2), b = s(1,2), f = s(2,2);
	dEdh1 = (2*(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3 + a*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + 
		((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(-((h2 + a*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - 
		(h2 + a*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + (2*h1 + 2*a*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*
		(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*
		(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) +(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*
		(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2))*(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*
		(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) +(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*
		(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + 2*(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3 + a*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + 
		((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(-((h2 + a*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - 
		(h2 + a*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + (2*h1 + 2*a*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*
		(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*
		(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + (h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*
		(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - ((-h2 - a*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

	dEdh2 = (2*(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*h2 + 2*a*h8) - 
		(h1 + a*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - (h1 + a*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))
		)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - 
		((2*h2 + 2*a*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + 
		2*(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*
		((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*h2 + 2*a*h8) - 
		(h1 + a*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - (h1 + a*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))
		)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - 
		((-h1 - a*h7)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

	dEdh3 = ((-2*(h1 + a*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))) - 
		(2*(h1 + a*h7)*(-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

	dEdh4 = (2*(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h6 + b*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + 
		((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		(-((h5 + b*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - 
		(h5 + b*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + 
		(2*h4 + 2*b*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2))*
		(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + 
		2*(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h6 + b*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + 
		((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*
		(-((h5 + b*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - 
		(h5 + b*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + 
		(2*h4 + 2*b*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - 
		((-h5 - b*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

	dEdh5 = (2*(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*h5 + 2*b*h8) - 
		(h4 + b*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - (h4 + b*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))
		)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - 
		((2*h5 + 2*b*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + 
		2*(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*
		((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*h5 + 2*b*h8) - 
		(h4 + b*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - (h4 + b*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))
		)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - 
		((-h4 - b*h7)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

	dEdh6 = ((-2*(h4 + b*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))) - 
		(2*(h4 + b*h7)*(-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

	dEdh7 = (2*(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(a*h3 + b*h6 + f*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + 
		((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		(-((a*h2 + b*h5 + f*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - 
		(a*h2 + b*h5 + f*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + 
		(2*a*h1 + 2*b*h4 + 2*f*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*
		(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2))*
		(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + 
		2*(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(a*h3 + b*h6 + f*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) + 
		((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*
		(-((a*h2 + b*h5 + f*h8)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)) - 
		(a*h2 + b*h5 + f*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)) + 
		(2*a*h1 + 2*b*h4 + 2*f*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))*
		(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - 
		((-(a*h2) - b*h5 - f*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

	dEdh8 = (2*(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*a*h2 + 2*b*h5 + 2*f*h8) - 
		(a*h1 + b*h4 + f*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - 
		(a*h1 + b*h4 + f*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - 
		((2*a*h2 + 2*b*h5 + 2*f*h8)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u) + 
		2*(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*
		((h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(2*a*h2 + 2*b*h5 + 2*f*h8) - 
		(a*h1 + b*h4 + f*h7)*(h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8) - 
		(a*h1 + b*h4 + f*h7)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8)))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		Power<T>(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)),2) - 
		((-(a*h1) - b*h4 - f*h7)*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));

	dEdh9 = ((-2*(a*h1 + b*h4 + f*h7)*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*
		(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))) - 
		(2*(a*h1 + b*h4 + f*h7)*(-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*
		(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))))/
		(2.*Sqrt<T>(Power<T>(-(((h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - u,2) + 
		Power<T>(-(((-(h1*(h2 + a*h8)) - h4*(h5 + b*h8) - h7*(a*h2 + b*h5 + f*h8))*(h3*(h1 + a*h7) + h6*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h9))/
		(-((h2*(h1 + a*h7) + h5*(h4 + b*h7) + (a*h1 + b*h4 + f*h7)*h8)*(h1*(h2 + a*h8) + h4*(h5 + b*h8) + h7*(a*h2 + b*h5 + f*h8))) + 
		(h1*(h1 + a*h7) + h4*(h4 + b*h7) + h7*(a*h1 + b*h4 + f*h7))*(h2*(h2 + a*h8) + h5*(h5 + b*h8) + h8*(a*h2 + b*h5 + f*h8)))) - v,2)));
}

#endif