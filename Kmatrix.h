#ifndef KMATRIX_H
#define KMATRIX_H
#include "LM.h"
#include "matrix.h"
#include "SVD.h"
#include <vector>
using namespace libNumerics;

template <typename T> 
void errorKparams(const matrix<T> &K0,const  matrix<T> &K, T &ea, T &eb, T &eu, T &ev, T &eg) {
	ea = 100*std::abs(K0(0,0)-K(0,0))/K0(0,0);
	eb = 100*std::abs(K0(1,1)-K(1,1))/K0(1,1);
	eu = std::abs(K0(0,2)-K(0,2));
	ev = std::abs(K0(1,2)-K(1,2));
	eg = std::abs(K0(0,1)-K(0,1));
}

template <typename T>
matrix<T> composevij(matrix<T>& H, int i, int j) {
	matrix<T> vij(1,6);
	vij(0,0) = H(0,i)*H(0,j);
	vij(0,1) = H(0,i)*H(1,j)+H(1,i)*H(0,j);
	vij(0,2) = H(1,i)*H(1,j);
	vij(0,3) = H(2,i)*H(0,j)+H(0,i)*H(2,j);
	vij(0,4) = H(2,i)*H(1,j)+H(1,i)*H(2,j);
	vij(0,5) = H(2,i)*H(2,j); 
	return vij;
}


template <typename T>
matrix<T> extractKfromH(std::vector<matrix<T> > H)
{
	int n = H.size();
	if (n < 3)
		exit(1);
	matrix<T> V(2*n,6);
    //Leman: i<3 -> i<n
    for (int i = 0; i < n; i++) {
		matrix<T> h = H[i];
		matrix<T> v12 = composevij<T>(h,0,1);
		matrix<T> v11 = composevij<T>(h,0,0);
		matrix<T> v22 = composevij<T>(h,1,1);
		V.paste(i*2, 0, v12);
		V.paste(i*2+1, 0, v11-v22);
	}
	V = V.t()*V;
	vector<T> b(6);
	SVD<T>::Nullspace(V, &b);
	T v0 = (b[1]*b[3]-b[0]*b[4]) / (b[0]*b[2]-b[1]*b[1]);
	T lambda = b[5] - (b[3]*b[3] + v0*(b[1]*b[3]-b[0]*b[4])) / b[0];
	T alpha = std::sqrt(lambda/b[0]);
	T beta = std::sqrt(lambda*b[0] / (b[0]*b[2]-b[1]*b[1]));
	T gamma = -b[1]*alpha*alpha*beta/lambda;
	T u0 = gamma*v0/beta - b[3]*alpha*alpha/lambda;
	matrix<T> K(3,3); K.fill(0);
	K(0,0) = alpha; K(0,1) = gamma; K(0,2) = u0;
	K(1,1) = beta; K(1,2) = v0;
	K(2,2) = 1;
	return K;
}

template <typename T>
matrix<T> groundK(T alpha, T beta, T gamma, T u0, T v0) {
	matrix<T> K = matrix<T>::zeros(3,3);
	K(0,0) = alpha;
	K(0,1) = gamma;
	K(0,2) = u0;
	K(1,1) = beta;
	K(1,2) = v0;
	K(2,2) = 1;
	return K;
}

#endif
