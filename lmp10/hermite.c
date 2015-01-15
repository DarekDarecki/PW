#include "makespl.h"
#include "piv_ge_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
zmiennej środowiskowej APPROX_BASE_SIZE
*/

double hermite(double x, int n){
	if(n == 0)
		return 1;
	if(n == 1)
		return 2 * x;
	if(n >= 2)
		return 2 * x * hermite(x, n - 1) - 2 * (n - 1) * hermite(x, n - 2);
}

double hermiteD1(double x, int n){
	if(n == 1)	
		return 2;
	if(n >= 2){
		return 2 * hermite(x, n - 1) + 2 * x * hermiteD1(x, n - 1) - 2 * (n - 1) * hermiteD1(x, n - 2);
	} else return 0;
}

double hermiteD2(double x, int n){
	if(n == 2)
		return 8 ;
	if(n >= 3){
		return 4 * hermiteD1(x, n - 1) + 2 * x * hermiteD2(x, n - 1) - 2*(n-1) * hermiteD2(x, n - 2);
	} else return 0;
}

double hermiteD3(double x, int n){
	if(n == 3)
		return 48;
	if(n >= 4){
		return 6 * hermiteD2(x, n - 1) + 2 * x * hermiteD3(x, n - 1) - 2 * (n - 1) * hermiteD3(x, n - 2);
	} else return 0;
}

void wypelnijMacierze(matrix_t *eqs, points_t *pts, int n){
    double *x = pts->x;
    double *y = pts->y;
    int i, j, k;
    double suma;
    for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++){
                    for (k = 0; k < pts->n; k++)
                            suma += hermite(x[k], i) * hermite(x[k], j);
                add_to_entry_matrix(eqs, i, j, suma);
                suma = 0;
            }
            for (k = 0; k < pts->n; k++)
                suma += y[k] * hermite(x[k], i);
            add_to_entry_matrix(eqs, i, n, suma);
            suma = 0;
        }
}	

void make_spl(points_t * pts, spline_t * spl){
	matrix_t *eqs= NULL;
	double *x = pts->x;
	double *y = pts->y;
	int	i, j, k;
	int	n = pts->n > 8 ? 8 : pts->n;
	double		a = x[0];
	double		b = x[pts->n - 1];
	
	char *nEnv= getenv( "APPROX_BASE_SIZE" );
	if( nEnv != NULL && atoi( nEnv ) > 0 )
		n = atoi( nEnv );
	
	eqs = make_matrix(n, n+1);
	
	wypelnijMacierze(eqs, pts, n);
	
	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
	
	if (alloc_spl(spl, n) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i * (b - a)/(spl->n - 1);
			xx+= 10.0 * DBL_EPSILON; // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < n; k++) {
				double	ck = get_entry_matrix(eqs, k, n);
				spl->f[i] += ck * hermite(xx, k);
				spl->f1[i] += ck * hermiteD1(xx, k);
				spl->f2[i] += ck * hermiteD2(xx, k);
				spl->f3[i] += ck * hermiteD3(xx, k);
			}
		}
	}
}