#include <math.h> // required for sqrt(), fabs();
#include <stdio.h>  // required for exit
#include <stdlib.h> // required for malloc(), free();
#include <Rinternals.h>  // required for SEXP et.al.;
#include <float.h>  // required for DBL_EPSILON
#include <string.h> // required for memcpy()
#include "var_hdh.h"

/*void mh2_rrn_test(double *sigma2, double *h2, double *x1, double *y, double *l, double eta_k, int n, int p, int q, int *debias){
	int i, j, k;

	double An1_tr = 0.0, YAn1Y, bias = 0.0, vy0, lYAn1Yl, Ylnorm, sigma2tmp, yl0, vyl0;
    double **xtx 	= requirematrix(n, n, 0);
    double *yl 	= (double*)malloc(sizeof(double)*n);
	double *lam 	= (double*)malloc(sizeof(double)*n);		
	double **x 		= (double**)malloc(sizeof(double*)*n);
	double *vy 		= (double*)malloc(sizeof(double)*n*q);
    double *vyl 	= (double*)malloc(sizeof(double)*n);

	for (j = 0; j < n; ++j){
	    yl0=0.0;
	    for (k = 0; k < q; k++) yl0 += y[k*n+j]*l[k];
	    yl[j] = yl0;
	}

	for (i = 0; i < n; i++) x[i] = &x1[i*p];

	// matrixprod(y, l, yl, n, q, 1, 0 );
	matrixprod(x, x, xtx, n, p, n, 2);
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) xtx[i][j] /= n;
	}
	eigencmp(xtx, lam, n, 1);

	for (j = 0; j < n; ++j)	An1_tr += 1/(eta_k+lam[j]);

	for (k = 0; k < q; k++){
		for (j = 0; j < n; ++j){		
			vy0 = 0.0;
			for (i = 0; i < n; ++i) vy0 += xtx[i][j]*y[k*n+i];
			vy[k*n+j] = vy0;
		}
	}

    for (j = 0; j < n; ++j){
        vyl0=0.0;
        for (k = 0; k < q; k++) vyl0 += vy[k*n+j]*l[k];
        vyl[j] = vyl0;
    }


	for (k = 0; k < q; k++){
		for (i = k; i < q; i++){
			YAn1Y	= 0.0;
			for (j = 0; j < n; ++j) {
				YAn1Y += vy[k*n+j]*vy[i*n+j]/(eta_k+lam[j]);
			}
			if (debias[0] == 1) bias = eta_k*YAn1Y/An1_tr;			
			sigma2[k*q+i] = YAn1Y/An1_tr - bias;
		}
	}
    lYAn1Yl = 0.0;
	for (j = 0; j < n; ++j){
        lYAn1Yl += vyl[j]*vyl[j]*lam[j]/(eta_k+lam[j]);
    }
	Ylnorm = 0.0;
    for (j = 0; j < n; ++j){
        Ylnorm += yl[j]*yl[j];
    }
    Ylnorm /=n;

	for (k = 1; k < q; k++){
		for (i = 0; i < k; i++){		
			sigma2[k*q+i] = sigma2[i*q+k];
		}
	}
	sigma2tmp = 0.0;
    sigma2tmp = (Ylnorm*n - lYAn1Yl)/(n-An1_tr);
	h2[0] = 1-sigma2tmp/Ylnorm;
    freematrix(xtx, n);
	free(lam);
	free(x);
	free(vyl);
	free(yl);
	
}*/

void mh2_rrn_test(double *sigma2, double *h2, double *x1, double *y, double *l,double eta_k, int n, int p, int q, int *debias){
    int i, j, k;

    double An1_tr = 0.0, YAn1Y, bias = 0.0, vy0, yl0, lsigmal, Ylnorm;
    double **xtx 	= requirematrix(n, n, 0);
    double *lam 	= (double*)malloc(sizeof(double)*n);
    double **x 		= (double**)malloc(sizeof(double*)*n);
    double *vy 		= (double*)malloc(sizeof(double)*n*q);
    double *yl 	= (double*)malloc(sizeof(double)*n);

    for (j = 0; j < n; ++j){
        yl0=0.0;
        for (k = 0; k < q; k++) yl0 += y[k*n+j]*l[k];
        yl[j] = yl0;
    }
    Ylnorm = 0.0;
    for (j = 0; j < n; ++j) {
        Ylnorm += yl[j] * yl[j];
    }
    Ylnorm /=n;

    for (i = 0; i < n; i++) x[i] = &x1[i*p];

    matrixprod(x, x, xtx, n, p, n, 2);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) xtx[i][j] /= n;
    }

    eigencmp(xtx, lam, n, 1);

    for (j = 0; j < n; ++j)	An1_tr += 1/(eta_k+lam[j]);
    // trA1[0]	= eta_k*An1_tr;  // n - tr(A1n)

    for (k = 0; k < q; k++){
        for (j = 0; j < n; ++j){
            vy0 = 0.0;
            for (i = 0; i < n; ++i) vy0 += xtx[i][j]*y[k*n+i];
            vy[k*n+j] = vy0;
        }
    }

    for (k = 0; k < q; k++){
        for (i = k; i < q; i++){
            YAn1Y	= 0.0;
            for (j = 0; j < n; ++j) {
                YAn1Y += vy[k*n+j]*vy[i*n+j]/(eta_k+lam[j]);
            }
            if (debias[0] == 1) bias = eta_k*YAn1Y/An1_tr;
            sigma2[k*q+i] = YAn1Y/An1_tr - bias;
        }
    }
    for (k = 1; k < q; k++){
        for (i = 0; i < k; i++){
            sigma2[k*q+i] = sigma2[i*q+k];
        }
    }
    lsigmal = 0.0;
    for (i = 0; i < q; i++){
        for (j = 0 ; j < q; j++){
            lsigmal += l[i]*sigma2[i*q+j]*l[j];
        }
    }

    h2[0] = 1-(lsigmal/Ylnorm);
    freematrix(xtx, n);
    free(lam);
    free(x);
    free(vy);

}
/*void mh2_rrn_test(double *sigma2, double *trA1, double *h2, double *normly, double *x1, double *y, double *l, double eta_k, int n, int p, int q, int *debias){
    int i, j, k;
    double An1_tr = 0.0,  YAn1Y, fenzi = 0.0, bias = 0.0, vy0,tmp;
    double **xtx 	= requirematrix(n, n, 0);  // **代表指针的指针，这句相当于在声明这个变量
    double *lam 	= (double*)malloc(sizeof(double)*n);
    double **x 		= (double**)malloc(sizeof(double*)*n);
    double *vy 		= (double*)malloc(sizeof(double)*n*q);


    for (i = 0; i < n; i++) x[i] = &x1[i*p];

    for (j=0;j<n;j++){
        tmp  =  0.0 ;
        for (k=0;k<q;k++) tmp += y[k*n+j]*l[k];
        normly[0] += tmp*tmp;
    }

    matrixprod(x, x, xtx, n, p, n, 2); // transpose=0:X*Y;transpose=1:t(X)*Y;transpose=2：X*t(Y)
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) xtx[i][j] /= n;
    }

    eigencmp(xtx, lam, n, 1);

    for (j = 0; j < n; ++j)	An1_tr += 1/(eta_k+lam[j]);
    trA1[0]	= eta_k*An1_tr;  // n - tr(A1n)

    for (k = 0; k < q; k++){
        for (j = 0; j < n; ++j){
            vy0 = 0.0;
            for (i = 0; i < n; ++i) vy0 += xtx[i][j]*y[k*n+i];
            // for (i = 0; i < n; ++i) vy0 += xtx[i][j]*ly[i];
            vy[k*n+j] = vy0;
        }
    }

    for (k = 0; k < q; k++){
        for (i = k; i < q; i++){
            YAn1Y	= 0.0;
            for (j = 0; j < n; ++j) {
                YAn1Y += vy[k*n+j]*vy[i*n+j]/(eta_k+lam[j]);
            }
            if (debias[0] == 1) bias = eta_k*YAn1Y/An1_tr;
            sigma2[k*q+i] = YAn1Y/An1_tr - bias;
        }
    }
    for (k = 1; k < q; k++){
        for (i = 0; i < k; i++){
            sigma2[k*q+i] = sigma2[i*q+k];
        }
    }


    for (k = 0; k < q; k++){
        for (j = 0; j< q; j++){
            fenzi +=l[k]*sigma2[k*q+j]*l[j];
        }
    }
    h2[0] = 1 - fenzi/(normly[0]/n);
    freematrix(xtx, n);
    free(lam);
    free(x);
    free(vy);

}*/

void mh2_rrp_test(double *sigma2, double *h2, double *x1, double *y, double *l, double eta_k, int n, int p, int q, int is_eta){
	int i, j, k;

	double An1_tr = 0.0, tmp, normqy, normrqy, yl0, lsigmal, Ylnorm;
	double *x, *Q, *R, *qy, *invR, *rqy;
    double *yl 	= (double*)malloc(sizeof(double)*n);
	x 		= (double*)malloc(sizeof(double)*n*p); 
	Q 		= (double*)malloc(sizeof(double)*n*p); 
	R 		= (double*)malloc(sizeof(double)*p*p); 
	invR 	= (double*)malloc(sizeof(double)*p*p);
	qy 		= (double*)malloc(sizeof(double)*p*q);
	rqy 	= (double*)malloc(sizeof(double)*p*q);


    for (j = 0; j < n; ++j){
        yl0=0.0;
        for (k = 0; k < q; k++) yl0 += y[k*n+j]*l[k];
        yl[j] = yl0;
    }
    Ylnorm = 0.0;
    for (j = 0; j < n; ++j) {
        Ylnorm += yl[j] * yl[j];
    }
    Ylnorm /=n;

	for(j=0;j<p;j++){
		for(i=0;i<n;i++)	x[j*n+i]	= x1[i*p+j];
	}
	QRDecompN(Q, R, x, n, p);
	LowTriangularInv(invR, p, R);	

	An1_tr = n - p;
	if(is_eta){
		tmp = 0.0;
		for(i=0;i<p;i++)
			for(j=0;j<=i;j++)
				tmp	+= invR[i*p+j]*invR[i*p+j];		
		An1_tr	+= n*eta_k*tmp;	
	}

	
	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			tmp = 0.0; 
			for(i=0;i<n;i++) tmp += Q[j*n+i]*y[k*n+i];
			qy[k*p+j] = tmp;					
		}
	}

	if(is_eta){
		for(k=0;k<q;k++){
			for(j=0;j<p;j++){
				tmp = 0.0;
				for(i=j;i<p;i++)	tmp	+= invR[i*p+j]*qy[k*q+i];
				rqy[k*p+j]	= tmp;	
			}	
		}		
	}


	for(k=0;k<q;k++){
		for(j=k;j<q;j++){
			tmp = normqy = 0.0; 
			for(i=0;i<n;i++) tmp 	+= y[j*n+i]*y[k*n+i];
			for(i=0;i<p;i++) normqy += qy[k*p+i]*qy[j*p+i];
			sigma2[k*q+j]	= (tmp - normqy)/An1_tr;

			if(is_eta){
				normrqy = 0.0;
				for(i=0;i<p;i++)	normrqy += rqy[k*p+i]*rqy[j*p+i];			
				sigma2[k*q+j]	+=  n*eta_k*normrqy/An1_tr;
			}
		}
	}
	
	for (k = 1; k < q; k++){
		for (i = 0; i < k; i++){		
			sigma2[k*q+i] = sigma2[i*q+k];
		}
	}

	lsigmal = 0.0;
	for (i = 0; i < q; i++){
        for (j = 0 ; j < q; j++){
            lsigmal += l[i]*sigma2[i*q+j]*l[j];
        }
	}

    h2[0] = 1-(lsigmal/Ylnorm);

	free(x);
	free(Q);
	free(R);
	free(invR);
	free(qy);	
	free(rqy);
	free(yl);

}

void h2_rrn_test(double *sigma2, double *h2, double *x1, double *y, double eta_k, int n, int p, int q, int *debias){
    int i, j, k;

    double An1_tr = 0.0, YAn1Y, bias = 0.0, vy, normy;
    double **xtx 	= requirematrix(n, n, 0);
    double *lam 	= (double*)malloc(sizeof(double)*n);
    double **x 		= (double**)malloc(sizeof(double*)*n);

    for (i = 0; i < n; i++) x[i] = &x1[i*p];

    matrixprod(x, x, xtx, n, p, n, 2);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) xtx[i][j] /= n;
    }

    eigencmp(xtx, lam, n, 1);

    for (j = 0; j < n; ++j)	An1_tr += 1/(eta_k+lam[j]);
    // trA1[0]	= eta_k*An1_tr;


    for (k = 0; k < q; k++){
        YAn1Y	= 0.0;
        normy = 0.0;
        for(i=0;i<n;i++)	normy += y[k*n+i]*y[k*n+i];
        for (j = 0; j < n; ++j) {
            vy = 0.0;
            for (i = 0; i < n; ++i) vy += xtx[i][j] * y[k * n + i];
            YAn1Y += vy*vy/(eta_k+lam[j]);
        }
        if (debias[0] == 1) bias = eta_k*YAn1Y/An1_tr;
        sigma2[k] = YAn1Y/An1_tr - bias;
        normy /= n;
        h2[k] = 1 - sigma2[k]/normy;
    }

    freematrix(xtx, n);
    free(lam);
    free(x);
}

void h2_rrp_test(double *sigma2, double *h2, double *x1, double *y, double eta_k, int n, int p, int q, int is_eta){
    int i, j, k;

    double An1_tr = 0.0, tmp, normy, normqy, normrqy;
    double *x, *Q, *R, *qy, *invR;
    x 		= (double*)malloc(sizeof(double)*n*p);
    Q 		= (double*)malloc(sizeof(double)*n*p);
    R 		= (double*)malloc(sizeof(double)*p*p);
    qy 		= (double*)malloc(sizeof(double)*p);
    invR 	= (double*)malloc(sizeof(double)*p*p);

    for(j=0;j<p;j++){
        for(i=0;i<n;i++)	x[j*n+i]	= x1[i*p+j];
    }
    QRDecompN(Q, R, x, n, p);
    LowTriangularInv(invR, p, R);

    An1_tr = n - p;
    if(is_eta){
        tmp = 0.0;
        for(i=0;i<p;i++)
            for(j=0;j<=i;j++)
                tmp	+= invR[i*p+j]*invR[i*p+j];

        An1_tr	+= n*eta_k*tmp;
    }
    // trA1[0]	= An1_tr;
    for(k=0;k<q;k++){
        normqy = 0.0;
        for(j=0;j<p;j++){
            tmp = 0.0;
            for(i=0;i<n;i++) tmp += Q[j*n+i]*y[k*n+i];
            qy[j] = tmp;   //
            normqy	+= tmp*tmp;
        }
        normy = 0.0;
        for(i=0;i<n;i++)	normy += y[k*n+i]*y[k*n+i];
        sigma2[k]	= (normy - normqy)/An1_tr;

        if(is_eta){
            normrqy = 0.0;
            for(j=0;j<p;j++){
                tmp = 0.0;
                for(i=0;i<=j;i++)	tmp	+= invR[j*p+i]*qy[i];
                normrqy	+= tmp*tmp;
            }
            sigma2[k]	+=  n*eta_k*normrqy/An1_tr;
        }
        normy /= n;
        h2[k] = 1 - sigma2[k]/normy;
    }
    free(x);
    free(Q);
    free(R);
    free(invR);
    free(qy);
}

SEXP h2_RR_test(SEXP X_, SEXP Y_, SEXP ETAK_, SEXP PARA_, SEXP DEBIAS_, SEXP ALPHA_)
{
    // dimensions
    int *para = INTEGER(PARA_);
    int n     = para[0];
    int p     = para[1];
    int q     = para[2];
    int is_eta= para[3];

    // Pointers
    double *y	= REAL(Y_);
    double *x	= REAL(X_);
    double eta	= 0.0;

    // Outcome
    SEXP _output, _sigma2,  _h2, _r_names;
    PROTECT(_output 	= allocVector(VECSXP, 	2));
    PROTECT(_r_names   	= allocVector(STRSXP, 	2));
    // PROTECT(_trA1 		= allocVector(REALSXP, 	1));
    PROTECT(_sigma2 	= allocVector(REALSXP, 	q));
    PROTECT(_h2 	    = allocVector(REALSXP, 	q));

    if(p<n){
        if(is_eta)	eta = REAL(ETAK_)[0];
        h2_rrp_test(REAL(_sigma2), REAL(_h2), x, y, eta, n, p, q, is_eta);
    }
    else{
        if(is_eta)	eta = REAL(ETAK_)[0];
        else	eta = setup_eta(y, x, REAL(ALPHA_), para);
        h2_rrn_test(REAL(_sigma2), REAL(_h2), x, y, eta, n, p, q, INTEGER(DEBIAS_));
    }


    SET_STRING_ELT(_r_names, 	0,	mkChar("sigma2"));
    SET_STRING_ELT(_r_names, 	1,	mkChar("h2"));
    SET_VECTOR_ELT(_output,		0, 	_sigma2);
    SET_VECTOR_ELT(_output,		1, 	_h2);
    setAttrib(_output, 			R_NamesSymbol, _r_names);

    UNPROTECT(4);
    return _output;
}

SEXP mh2_RR_test(SEXP X_, SEXP Y_, SEXP L_, SEXP ETAK_, SEXP PARA_, SEXP DEBIAS_, SEXP ALPHA_)
{
	// dimensions
	int *para = INTEGER(PARA_);
	int n     = para[0];
	int p     = para[1];
	int q     = para[2];
	int is_eta= para[3];

	// Pointers
	double *y	= REAL(Y_);
	double *x	= REAL(X_);
	double *l   = REAL(L_);
	double eta	= 0.0;

	// Outcome
	SEXP _output, _sigma2, _h2, _r_names;
	PROTECT(_output 	= allocVector(VECSXP, 	2));
	PROTECT(_r_names   	= allocVector(STRSXP, 	2));
	PROTECT(_h2 		= allocVector(REALSXP, 	1));
	PROTECT(_sigma2 	= allocVector(REALSXP, 	q*q));

	if(p<n){	
		if(is_eta)	eta = REAL(ETAK_)[0];
		mh2_rrp_test(REAL(_sigma2), REAL(_h2), x, y, l, eta, n, p, q, is_eta);
	}
	else{
		if(is_eta)	eta = REAL(ETAK_)[0];
		else	eta = setup_eta_cov(y, x, REAL(ALPHA_), para);
		mh2_rrn_test(REAL(_sigma2), REAL(_h2), x, y, l, eta, n, p, q, INTEGER(DEBIAS_));
	}
	

	SET_STRING_ELT(_r_names, 	0,	mkChar("sigma2"));
	SET_STRING_ELT(_r_names, 	1,	mkChar("h2"));
	SET_VECTOR_ELT(_output,		0, 	_sigma2);
	SET_VECTOR_ELT(_output,		1, 	_h2);
	setAttrib(_output, 			R_NamesSymbol, _r_names); 

	UNPROTECT(4);
	return _output;
}

SEXP h2_MM_test(SEXP Y_, SEXP X_, SEXP SIGMA_, SEXP PARA_, SEXP IDENTITY_, SEXP ISNULL_)
{
    // dimensions
    int *para = INTEGER(PARA_);
    int n     = para[0];
    int p     = para[1];

    // Pointers
    double **x = (double**)malloc(sizeof(double*)*n);
    double *y  = REAL(Y_);

    // Intermediate quantities
    int i, j;
    double ynorm = 0.0, xynorm = 0.0;
    double *temp = (double*)calloc(p, sizeof(double));

    // Outcome
    SEXP _output, _sigma2, _h2, _r_names;
    PROTECT(_output = allocVector(VECSXP, 2));
    PROTECT(_sigma2 = allocVector(REALSXP, 1));
    PROTECT(_h2 = allocVector(REALSXP, 1));
    PROTECT(_r_names   = allocVector(STRSXP, 2));

    for (i = 0; i < n; ++i) {
        x[i] = &REAL(X_)[i*p];
        ynorm += pow(y[i], 2);
    }
    for (j = 0; j < p; ++j) {
        for (i = 0; i < n; ++i) temp[j] += x[i][j]*y[i];
    }
    ynorm /= n;

    if (INTEGER(ISNULL_)[0] == 0) {
        double temp1;
        for (i = 0; i < p; ++i) {
            for (temp1 = j = 0; j < p; ++j) temp1 += REAL(SIGMA_)[j*p+i]*temp[j];
            xynorm += temp1*temp[i];
        }
        REAL(_sigma2)[0] = ((p+n+1)*ynorm - xynorm/n)/(n+1);
        REAL(_h2)[0] = 1 - REAL(_sigma2)[0]/ynorm;
    } else {
        if (INTEGER(IDENTITY_)[0] != 0) {
            for (j = 0; j < p; ++j) xynorm += temp[j]*temp[j];
            REAL(_sigma2)[0] = ((p+n+1)*ynorm - xynorm/n)/(n+1);
        } else {
            double m1 = 0.0, m2 = 0.0;
            for (j = 0; j < p; ++j) xynorm += pow(temp[j], 2);
            double **xtx = requirematrix(n, n, 0);
            double *lam = (double*)malloc(sizeof(double)*n);

            matrixprod(x, x, xtx, n, p, n, 2);
            eigencmp(xtx, lam, n, 0);
            for (i = 0; i < n; ++i) {
                lam[i] /= n;
                m1 += lam[i];
                m2 += pow(lam[i], 2);
            }
            m2 = m2/p - pow(m1, 2)/p/n;
            m1 /= p;
            REAL(_sigma2)[0] = (1+p*pow(m1, 2)/m2/(n+1))*ynorm - m1*xynorm/n/m2/(n+1);
            REAL(_h2)[0] = 1 - REAL(_sigma2)[0]/ynorm;
            freematrix(xtx, n);
            free(lam);
        }
    }

    SET_STRING_ELT(_r_names, 0,  mkChar("sigma2"));
    SET_STRING_ELT(_r_names, 1,  mkChar("h2"));
    SET_VECTOR_ELT(_output, 0, _sigma2);
    SET_VECTOR_ELT(_output, 1, _h2);
    setAttrib(_output, R_NamesSymbol, _r_names);

    free(x);
    free(temp);
    UNPROTECT(4);
    return _output;
}

SEXP h2_MLE_test(SEXP Y_, SEXP X_, SEXP PARA_, SEXP MAXITER_, SEXP TOL_)
{
    // dimensions
    int *para = INTEGER(PARA_);
    int n     = para[0];
    int p     = para[1];

    // Pointers
    double **x = (double**)malloc(sizeof(double*)*n);
    double *y  = REAL(Y_);

    // Intermediate quantities
    int i, j, step = 1;
    double sigma2, eta2, likelih0, likelih;
    double yUy, trUxtx, yUxtxUy, trUxtxUxtx, yUxtxUxtxUy, theta1;
    double **H = requirematrix(2, 2, 0);
    double **xtx = requirematrix(n, n, 0);
    double *lam = (double*)malloc(sizeof(double)*n);
    double *vy = (double*)calloc(n, sizeof(double));
    double *vy2 = (double*)malloc(sizeof(double)*n);
    double *temp = (double*)malloc(sizeof(double)*n);

    double *theta0 = (double*)malloc(sizeof(double)*2);
    double *theta = (double*)malloc(sizeof(double)*2);
    double *grad = (double*)malloc(sizeof(double)*2);


    for (i = 0; i < n; ++i) x[i] = &REAL(X_)[i*p];
    theta0[0] = log(var(y, n, 1));
    theta0[1] = 0.01;
    sigma2 = theta0[0];
    eta2 = theta0[1];

    matrixprod(x, x, xtx, n, p, n, 2);
    eigencmp(xtx, lam, n, 1);
    for (i = 0; i < n; ++i) {
        lam[i] /= p;
        for (j = 0; j < n; ++j) vy[i] += xtx[j][i]*y[j];
        vy2[i] = pow(vy[i], 2);
        temp[i] = eta2*lam[i]+1;
    }
    for (likelih0 = i = 0; i < n; ++i) likelih0 += vy2[i]/(temp[i]);
    likelih0 *= exp(-1.0*sigma2);
    likelih0 += n*sigma2;
    for (i = 0; i < n; ++i) likelih0 += log(temp[i]);

    while(step < INTEGER(MAXITER_)[0]) {
        step += 1;
        yUy = trUxtx = yUxtxUy = trUxtxUxtx = yUxtxUxtxUy = 0.0;
        for (i = 0; i < n; ++i) {
            yUy += vy2[i]/temp[i];
            trUxtx += lam[i]/temp[i];
            yUxtxUy += vy2[i]*lam[i]/pow(temp[i], 2);
            trUxtxUxtx += pow(lam[i]/temp[i], 2);
            yUxtxUxtxUy += pow(vy[i]*lam[i], 2)/pow(temp[i], 3);
        }
        theta1 = exp(-1.0*sigma2);
        grad[0] = n-theta1*yUy;
        grad[1] = trUxtx-theta1*yUxtxUy;
        H[0][0] = theta1*yUy;
        H[0][1] = theta1*yUxtxUy;
        H[1][1] = -1.0*trUxtxUxtx + 2*theta1*yUxtxUxtxUy;
        grad[1] -= grad[0]*yUxtxUy/yUy;
        H[1][1] -= H[0][1]*yUxtxUy/yUy;
        theta[1] = theta0[1] - grad[1]/H[1][1];
        theta[0] = theta0[0] - (grad[0]-H[0][1]*grad[1]/H[1][1])/H[0][0];
        sigma2 = theta[0];
        eta2 = theta[1];
        for (i = 0; i < n; ++i) temp[i] = eta2*lam[i]+1;
        for (likelih = i = 0; i < n; ++i) likelih += vy2[i]/(temp[i]);
        likelih *= exp(-1.0*sigma2);
        likelih += n*sigma2;
        for (i = 0; i < n; ++i) likelih += log(temp[i]);
        if (likelih > likelih0) {
            sigma2 = theta0[0];
            eta2 = theta0[1];
            break;
        } else {
            if (FMAX(fabs(theta[0]-theta0[0]), fabs(theta[1]-theta0[1])) < REAL(TOL_)[0]) {
                break;
            } else {
                theta0[0] = theta[0];
                theta0[1] = theta[1];
                likelih0 = likelih;
            }
        }
    }

    freematrix(xtx, n);
    freematrix(H, 2);
    free(lam);
    free(vy);
    free(vy2);
    free(temp);
    free(x);
    free(theta);
    free(theta0);
    free(grad);

    // Outcome
    SEXP _output, _sigma2, _h2, _r_names;
    PROTECT(_output = allocVector(VECSXP, 2));
    PROTECT(_sigma2 = allocVector(REALSXP, 1));
    PROTECT(_h2 = allocVector(REALSXP, 1));
    PROTECT(_r_names   = allocVector(STRSXP, 2));
    REAL(_sigma2)[0] = exp(sigma2);
    REAL(_h2)[0] = eta2/(eta2+1);

    SET_STRING_ELT(_r_names, 0,  mkChar("sigma2"));
    SET_STRING_ELT(_r_names, 1,  mkChar("h2"));
    SET_VECTOR_ELT(_output, 0, _sigma2);
    SET_VECTOR_ELT(_output, 1, _h2);
    setAttrib(_output, R_NamesSymbol, _r_names);

    UNPROTECT(4);
    return _output;
}