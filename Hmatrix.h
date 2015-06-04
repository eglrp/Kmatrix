#ifndef HMATRIX_H
#define HMATRIX_H

#include "LM.h"
#include "matrix.h"
#include "SVD.h"
#include "jac_eqs.h"
#include <vector>
#include <cmath>
using namespace libNumerics;

// HCS == H0CS
template <typename T>
T evaluateHomography(matrix<T> &H, std::vector<matrix<T> > &S, matrix<T> &H0CS)
{
	int n = S.size();
	vector<T> error(n);
	for (int i = 0; i < n; i++)
	{
		T x1 = H0CS(0,i), y1 = H0CS(1,i);
		matrix<T> s = S[i];
		T csx, csy;
		getEllipseCenter(s, csx, csy);
		matrix<T> cs = matrix<T>::zeros(3,1); cs(0,0) = csx; cs(1,0) = csy; cs(2,0) = 1;
		matrix<T> hcs = H*cs;
		T x2 = hcs(0,0)/hcs(2,0), y2 = hcs(1,0)/hcs(2,0);
		error[i] = std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	}
	return std::sqrt(error.qnorm()/n);
}

//template <typename T>
//void PreconditionerFromPoints(int width, int height, matrix<T> *N) 
//{
//	// Build the normalization matrix
//	double dNorm = 1.0 / sqrt( static_cast<double>(width*height) ); 
//
//	(*N) = matrix<T>::eye(3);
//	(*N)(0,0) = (*N)(1,1) = dNorm;
//	(*N)(0,2) = -.5f*width*dNorm;
//	(*N)(1,2) = -.5*height*dNorm;
//}
//
//template <typename T>
//void ApplyTransformationToPoints(const matrix<T> &points, const matrix<T> &N, matrix<T> *transformed_points) {
//	const int n = points.ncol();
//	(*transformed_points) = matrix<T>(2,n);
//	for (int i = 0; i < n; ++i) {
//		matrix<T> in(3,1), out(3,1);
//		in(0,0) = points(0, i); in(1,0) = points(1, i); in(2,0) = 1.0;
//		out = N * in;
//		(*transformed_points)(0, i) = out(0)/out(2);
//		(*transformed_points)(1, i) = out(1)/out(2);
//	}
//}
//
//template <typename T>
//void NormalizePoints(const matrix<T> &points, matrix<T> *normalized_points, matrix<T> *N, int width, int height) 
//{
//	PreconditionerFromPoints(width, height, N);
//	ApplyTransformationToPoints(points, *N, normalized_points);
//}
//
//template <typename T>
//void Unnormalize(const matrix<T> &T1, const matrix<T> &T2, matrix<T> *H)  {
//	*H = T2.t() * (*H) * T1;
//}

template <typename T>
matrix<T> solveHomography(matrix<T>& x1_, int w1, int h1, matrix<T>& x2_, int w2, int h2) {
	
	//matrix<T> x1_(x1.nrow(), x1.ncol()); 
	//matrix<T> x2_(x2.nrow(), x2.ncol()); 
	//matrix<T> N1_(3,3); 
	//matrix<T> N2_(3,3); 
	//assert(2 == x1_.nrow()); assert(x1_.nrow() == x2_.nrow()); assert(x1_.ncol() == x2_.ncol());
	//int w=std::max(w1,w2), h=std::max(h1,h2);
	//NormalizePoints(x1, &x1_, &N1_, w, h);
	//NormalizePoints(x2, &x2_, &N2_, w, h);
	
	int n = x1_.ncol();
	if (n < 4)
		exit(1);
	matrix<T> A = matrix<T>::zeros(n*2,9);
	for (int i = 0; i < n; ++i) {
		int j = 2*i;
		A(j,0) = x1_(0, i);
		A(j,1) = x1_(1, i);
		A(j,2) = 1.0;
		A(j,6) = -x2_(0, i) * x1_(0, i);
		A(j,7) = -x2_(0, i) * x1_(1, i);
		A(j,8) = -x2_(0, i);

		++j;
		A(j,3) = x1_(0, i);
		A(j,4) = x1_(1, i);
		A(j,5) = 1.0;
		A(j,6) = -x2_(1, i) * x1_(0, i);
		A(j,7) = -x2_(1, i) * x1_(1, i);
		A(j,8) = -x2_(1, i);
	}
	vector<T> vecNullspace(9);
	matrix<T> M(3,3); M.fill(0);
	if( SVD<T>::Nullspace(A, &vecNullspace) ) {
		M.read(vecNullspace);
		if(M.det() < 0)
			M = -M;
		M /= M(2,2);
		//Unnormalize<T>(N1_, N2_, &M);
	}
	return M;
}

// angles are given in radians
template <typename T>
matrix<T> groundH(matrix<T>& K, T thetax, T thetay, T thetaz, T tx, T ty, T tz)
{
	matrix<T> Rx = matrix<T>::zeros(3,3);
	matrix<T> Ry = matrix<T>::zeros(3,3);
	matrix<T> Rz = matrix<T>::zeros(3,3);

	Rx(0,0) = 1; 
	Rx(1,1) = std::cos(thetax); Rx(1,2) = -std::sin(thetax);
	Rx(2,1) = std::sin(thetax); Rx(2,2) = std::cos(thetax);

	Ry(0,0) = std::cos(thetay); Ry(0,2) = std::sin(thetay);
	Ry(1,1) = 1;
	Ry(2,0) = -std::sin(thetay); Ry(2,2) = std::cos(thetay);

	Rz(0,0) = std::cos(thetaz); Rz(0,1) = -std::sin(thetaz);
	Rz(1,0) = std::sin(thetaz); Rz(1,1) = std::cos(thetaz);
	Rz(2,2) = 1;

	matrix<T> R = Rx*Ry*Rz;
	matrix<T> H = matrix<T>::zeros(3,3);
	H.paste(0,0,R.col(0));
	H.paste(0,1,R.col(1));
	H(0,2) = tx; H(1,2) = ty; H(2,2) = tz;
	
	return K*H;
}

template <typename T>
T errorShape(const matrix<T>& H, const matrix<T>& S, T u, T v) {
	matrix<T> hs = H.t()*S*H;// because H is already inversed, we cant use rotateCircle(S,H);
	T x,y;
	getEllipseCenter(hs, x,y);
	return std::sqrt( (x-u)*(x-u) + (y-v)*(y-v) );
}

template <typename T>
void jacdif(const vector<T>& P, const std::vector<matrix<T> >& S, matrix<T>& J, matrix<T> UV)
{
	T EPS = 0.0001;
	matrix<T> H(3,3);
	T h1 = H(0,0) = P[0]; T h2 = H(0,1) = P[1]; T h3 = H(0,2) = P[2];
	T h4 = H(1,0) = P[3]; T h5 = H(1,1) = P[4]; T h6 = H(1,2) = P[5];
	T h7 = H(2,0) = P[6]; T h8 = H(2,1) = P[7]; T h9 = H(2,2) = P[8];
	int ncircle = J.nrow();
	for (int i = 0; i < 9; i++)
	{
		T dx = 0.25*EPS;
		vector<T> Pd(P);
		Pd[i] += dx;
		matrix<T> Hd(3,3);
		Hd(0,0) = Pd[0]; Hd(0,1) = Pd[1]; Hd(0,2) = Pd[2];
		Hd(1,0) = Pd[3]; Hd(1,1) = Pd[4]; Hd(1,2) = Pd[5];
		Hd(2,0) = Pd[6]; Hd(2,1) = Pd[7]; Hd(2,2) = Pd[8];
		for (int j = 0; j < ncircle; j++) {
			T r  = errorShape<T>(H, S[j], UV(0,j), UV(1,j));
			T rd = errorShape<T>(Hd, S[j], UV(0,j), UV(1,j));
			J(j,i) = rd-r/dx;
		}
	}
}

template <typename T>
class LMhomography : public MinLM<T> {
public:	
    LMhomography(std::vector<matrix<T> >& s, matrix<T>& uv) {
		S = s;
		UV = uv;
	}
private:
	matrix<T> Hini;
    std::vector<matrix<T> > S;
	matrix<T> UV;
public:
	virtual void modelData(const vector<T>& P, vector<T>& ymodel) const
	{
		int ncircle = ymodel.size();
		matrix<T> H(3,3);
		H(0,0) = P[0]; H(0,1) = P[1]; H(0,2) = P[2];
		H(1,0) = P[3]; H(1,1) = P[4]; H(1,2) = P[5];
		H(2,0) = P[6]; H(2,1) = P[7]; H(2,2) = P[8];
		for (int i = 0; i < ncircle; i++)
		{
			T err = errorShape(H, S[i], UV(0,i), UV(1,i));
			ymodel[i] = err;
		}
	}

	virtual void modelJacobian(const vector<T>& P, matrix<T>& J) const
	{
		//jacdif(P,S,J,UV);
		int ncircle = J.nrow();
		for (int i = 0; i < ncircle; i++){
			jac_inv_equations(P, S[i], UV(0,i), UV(1,i), J(i,0), J(i,1), J(i,2),
				J(i,3), J(i,4), J(i,5), J(i,6), J(i,7), J(i,8));
		}
	}

};

#endif
