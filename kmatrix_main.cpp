#include "ellipse_operations.h"
#include "matrix.h"
#include "Hmatrix.h"
#include "Kmatrix.h"
#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
using namespace libNumerics;

#define PI 3.14159265358979323846
//Leman:
#define RAD (PI/180)
#define DEG (180/PI)

//Leman random:
template <typename T>
T random(T a,T b){
    T r=T(rand())/RAND_MAX;
    return a+(b-a)*r;
}

template <typename T>
matrix<T> vec2mat(vector<T> v, int w, int h) {
    matrix<T> m(h,w);
    int n = v.size();
    assert(n == w*h);
    for (int i=0; i<w; i++) {
        for (int j=0; j<h; j++) {
            int idx = j*w+i;
            m(j,i) = v[idx];
        }
    }
    return m;
}

template <typename T>
vector<T> mat2vec(const matrix<T>& m, int l) {
    vector<T> v(l);
    int w = m.ncol();
    int h = m.nrow();
    assert(l == w*h);
    for (int i=0; i<w; i++) {
        for (int j=0; j<h; j++) {
            int idx = j*w+i;
            v[idx] = m(j,i);
        }
    }
    return v;
}

template <typename T>
std::vector<matrix<T> > virtualImage(T ncirclex, T ncircley, T x_cm, T y_cm,
                                     T cm1, T scale)
{
    int WI = x_cm*scale*cm1;
    int HE = y_cm*scale*cm1;
    int stepX = ceil((T)WI/(T)ncirclex);
    int stepY = ceil((T)HE/(T)ncircley);
    T radi = 0.3*std::min(stepX, stepY);
    int begin = 1, DELTA = 0;
    T XC = 0.5*(stepX*stepX-begin)/(stepX-begin) + DELTA;
    T YC = 0.5*(stepY*stepY-begin)/(stepY-begin) + DELTA;
    int ncircle = ncirclex*ncircley;
    std::vector<matrix<T> > S(ncircle);
    int cnt = 0;
    for (int j = 0; j < HE; j=j+stepY) {
        for (int i = 0; i < WI; i=i+stepX) {
            T XCi = XC + i;
            T YCi = YC + j;
            S[cnt] = getCircleMatrix(XCi, YCi, radi);
            cnt++;
        }
    }
    return S;
}

template <typename T>
std::vector<matrix<T> > buildCameras(matrix<T> &K0, T scale, T cm1)
{
    int nimage = 10;
    // [Matrix K]
    T alpha = 1250*scale, beta = 900*scale;
    T gamma = 1.09083;
    T u0 = 255*scale, v0 = 255*scale;
    // [Matrix K]

    matrix<T> K = groundK(alpha,beta,gamma,u0,v0);
    K0.paste(0,0,K);
    std::vector<matrix<T> > H(nimage);
    //Leman:
    printf("\nImage\tthetaX\tthetaY\tthetaZ\ttx\tty\ttz\n");
    for(int i=0;i<nimage;i++){
        T thetaX, thetaY, thetaZ, tx, ty, tz;
        if(true){//random data
            thetaX=random(-40*RAD,40*RAD);
            thetaY=random(-40*RAD,40*RAD);
            thetaZ=random(-40*RAD,40*RAD);
            tx=0.0*cm1*scale;
            ty=0.0*cm1*scale;
            tz=50.0*cm1*scale;
        }else{//given data
            thetaX=i*5*RAD;
            thetaY=0.0;
            thetaZ=0.0;
            tx=0.0*cm1*scale;
            ty=0.0*cm1*scale;
            tz=50.0*cm1*scale;
        }
        H[i] = groundH(K, thetaX,thetaY,thetaZ,tx,ty,tz);
        printf("%d\t% 3.1f\t% 3.1f\t% 3.1f",i,thetaX/RAD,thetaY/RAD,thetaZ/RAD);
        printf("\t% 3.1f\t% 3.1f\t% 3.1f\n",tx/cm1/scale,ty/cm1/scale,tz/cm1/scale);
    }
    //H[0] = groundH(K, 20*RAD, 0*RAD, 0*RAD, -9*cm1*scale, -12.5*cm1*scale, 50*cm1*scale);
    //H[1] = groundH(K, 0*RAD, 20*RAD, 0*RAD, -9*cm1*scale, -12.5*cm1*scale, 51*cm1*scale);
    //T ctmp = 1/std::sqrt((double)5);
    //H[2] = groundH(K, -30*RAD*ctmp, -30*RAD*ctmp, -15*RAD*ctmp,-10.5*cm1*scale, -12.5*cm1*scale, 52.5*cm1*scale);
    return H;
}

template <typename T>
matrix<T> extractK(std::vector<matrix<T> >& S, T scale, T cm1, T x_cm, T y_cm, T tolFun)
{
    int WI = x_cm*scale*cm1, HE = y_cm*scale*cm1;
    int wi = 512*scale, he = 512*scale;
    matrix<T> K0(3,3);
    std::vector<matrix<T> > H = buildCameras(K0, scale, cm1);
    int nimage = H.size();
    int ncircle = S.size();
    matrix<T> XY(2,ncircle);
    for (int i = 0; i < ncircle; i++) {
        T x,y;
        getEllipseCenter(S[i], x,y);
        XY(0,i) = x;
        XY(1,i) = y;
    }

    std::vector<matrix<T> > H_ellip(nimage);
    std::vector<matrix<T> > H_point(nimage);

    // calculate homographies through minimization of d|CHS-CH0S|
    printf("\nImage\tIter\tH_err_ini\tvs\tH_err_fin:\t\tDifference\n");
    for (int i = 0; i < nimage; i++) {
        matrix<T> CH0S = matrix<T>::zeros(2,ncircle);
        matrix<T> H0CS = matrix<T>::zeros(2,ncircle);
        rotateVirtualImg(S, H[i], CH0S, H0CS);
        matrix<T> UV(CH0S);

        matrix<T> H_ini = solveHomography(XY, WI, HE, UV, wi, he);
        H_point[i] = H_ini;
        LMhomography<T> HeLMA(S, UV);
        vector<T> trgData = vector<T>::zeros(ncircle);
        vector<T> h_ = mat2vec(H_ini.inv(),9);
        T rmse_ellip = HeLMA.minimize(h_, trgData, tolFun);
        matrix<T> H_fin = (vec2mat(h_,3,3)).inv();
        H_ellip[i] = H_fin;
        T rmse_ini = evaluateHomography<T>(H_ini, S, H0CS);
        T rmse_fin = evaluateHomography<T>(H_fin, S, H0CS);
        printf("%i\t%d",i,HeLMA.iterations);
        printf("\t%16.14f\t%16.14f\t%16.14f%%\n", rmse_ini, rmse_fin,rmse_fin/rmse_ini*100);
    }
    matrix<T> K_ellip = extractKfromH(H_ellip);
    matrix<T> K_point = extractKfromH(H_point);
    T eae, ebe, eue, eve, ege;
    errorKparams(K0, K_ellip, eae, ebe, eue, eve, ege);
    T eap, ebp, eup, evp, egp;
    errorKparams(K0, K_point, eap, ebp, eup, evp, egp);
    printf("\npoint vs ellip alfa:	%f	%f\n", eap,eae);
    printf("point vs ellip beta:	%f	%f\n", ebp,ebe);
    printf("point vs ellip    u:	%f	%f\n", eup,eue);
    printf("point vs ellip    v:	%f	%f\n", evp,eve);
    printf("point vs ellip gama:	%f	%f\n", egp,ege);

    return K_ellip;
}

int main(int argc, char ** argv)
{
    srand(time(NULL));for(int i=0;i<1000;i++)rand();
    int x_cm = 18, y_cm = 25;
    int cm1 = 50;
    int scale = 4;
    int ncirclex = 10, ncircley = 14;
    std::vector<matrix<double> > S = virtualImage<double>(ncirclex, ncircley, x_cm, y_cm, cm1, scale);
    matrix<double> K = extractK<double>(S, scale, cm1, x_cm, y_cm, 0.000001);

    return 0;
}


