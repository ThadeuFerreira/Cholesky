/**
 * @file main.c
 * @author  Heitor Freitas <heitor.freitas@gmail.com>
 * @author  Mariana Rodrigues <rodrigues.mariana@gmail.com>
 * @author  Thadeu A. F. Melo <thadeu.afm@gmail.com>
 * @version 1.0
 *
 * @license GPL
 *
 * @brief Matrix inversion through Cholesky Decomposition.
 */


/**
 ******************
 * Sources:
 ******************
 *
 * Cholesky:
 * http://www.icmc.usp.br/~andretta/ensino/aulas/sme0301-1-10/SistemasLinearesCholesky.pdf
 * http://wwwp.fc.unesp.br/~arbalbo/Iniciacao_Cientifica/sistemaslineares/teoria/3_Met_Cholesky.pdf
 * https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
 * http://rosettacode.org/wiki/Cholesky_decomposition#C
 * http://files.calcnum.webnode.com/200000021-605b36154f/Apostila%20-%20Cuminato.pdf
 *
 * Matrix:
 * http://rosettacode.org/wiki/Matrix_multiplication#C
 * http://www.programming-techniques.com/2011/09/numerical-methods-inverse-of-nxn-matrix.html
 *
 * GNU GCC Profiling Tool:
 * http://www.thegeekstuff.com/2012/08/gprof-tutorial/
 *
 * GCC optimizations:
 * ftp://gcc.gnu.org/pub/gcc/summit/2003/Optimizing%20for%20space.pdf
 */

#define PRINT_MATRICES 1 /*! @brief Configures if matrices are printed in screen. */
#define TIMED 0 /*! @brief Inserts timing checks in the program */
#define ERROR_ANALYSIS 2  /*!
                       	* @brief Computes mean error. Values are:
                       	* \n 0 No error is computed
                       	* \n 1 Only error in relation to Identity matrix is computed
                       	* \n 2 Both error in relation to Identity and original
                       	* matrices are computed.
                       	*/

int order = 0; /*! @brief Square matrix order               	*/
int n_int = 0; /*! @brief Number of iterations in routine sqrt_ */


#include <stdio.h>
#if(TIMED == 1)
#include <time.h>
#endif

typedef double fpoint_t; /*! @typedef fpoint_t @brief Variable type for easy
                          	* interchange between \c float and \c double
                          	* primitive types.
                          	*/
/**
 * Macro for debug purposes.
 */
#define DBG(x) printf("[%s:%d]%s():%s = %f\n",__FILE__,__LINE__,__FUNCTION__,#x,x);
int allocate(fpoint_t *A)
 {
        int i;
     for(i = 0; i < 25; i ++)
        A[i] = 0;

}
/*!
 * @brief Print matrix \b A in screen
 */
void show_matrix(fpoint_t * A) {
	int i, j;
	for (i = 0; i < order; i++) {
    	for (j = 0; j < order; j++) {
        	printf("%.9f ", A[i * order + j]);
    	}
    	printf("\n");
	}
}

/*!
 * @brief Calculates the square root of a floating-point number \b x.
 */
fpoint_t sqrt_(fpoint_t x) {
	fpoint_t val = 1.0;
	int times = n_int;
	while (times-- > 0) {
    	val = val / 2 + x / (val * 2);
	}
	return val;
}

/*!
 * @brief Computes in \b L the lower triangular matrix generated from Cholesky
 * decomposition of matrix \b A.
 * \n Modified from \url{http://rosettacode.org/wiki/Cholesky_decomposition#C}
 */
void cholesky(fpoint_t * A, fpoint_t * L) {
	int i, j, k;
	fpoint_t s;

	for (i = 0; i < order; i++) {
    	for (j = 0; j < (i + 1); j++) {
        	s = 0;
        	for (k = 0; k < j; k++) {
            	s += L[i * order + k] * L[j * order + k];
        	}
        	L[i * order + j] = (i == j) ? sqrt_(A[i * order + i] - s) : (1.0 / L[j * order + j] * (A[i * order + j] - s));
    	}
	}
}

/*!
 * @brief Computes the dot product between a row (pointed by \c fpoint_t \c *
 * \c row ) and a column (pointed by \c fpoint_t \c * \c col) of two square
 * matrix of order #order.
 * \n Modified from \url{http://rosettacode.org/wiki/Matrix_multiplication#C}
 */
fpoint_t dot(fpoint_t * row, fpoint_t * col) {
	fpoint_t r = 0;
	int len = order;
	while (len--) {
    	r += *row++ * *col;
    	col += order;
	}
	return r;
}

/*!
 * @brief Computes in \b C the multiplication of square matrices \b A and \b B (
 * \c C \c = \c A \c * \c B ).
 * \n Modified from \url{http://rosettacode.org/wiki/Matrix_multiplication#C}
 */
void multiply(fpoint_t * A, fpoint_t * B, fpoint_t * C) {
	fpoint_t * pa;
	int i, j;

	for (pa = A, i = 0; i < order; i++, pa += order) {
    	for (j = 0; j < order; j++) {
        	*C++ = dot(pa, B + j);
    	}
	}
}

/*!
 * @brief Stores in matrix \b B the multiplication between n-order square
 * matrix \b A and its transpose.
 * \n Modified from \url{http://rosettacode.org/wiki/Matrix_multiplication#C}
 */
void multiplyByTranspose(fpoint_t * A, fpoint_t * B) {
	int i, j, k;
	for (i = 0; i < order; i++) {
    	for (j = 0; j < order; j++) {
        	B[i * order + j] = (fpoint_t) 0.00;
        	for (k = 0; k < order; k++) {
            	B[i * order + j] += (A[i * order + k] * A[j * order + k]);
        	}
    	}
	}
}

/*!
 * @brief Stores in matrix \b AI the multiplication between the transpose of
 * square matrix \b G and the original \b G matrix.
 * \n Code based on \url{http://rosettacode.org/wiki/Matrix_multiplication#C}
 */
void inverseFromInvertedL(fpoint_t * AI, fpoint_t * G) {
	int i, j, k;
	for (i = 0; i < order; i++) {
    	for (j = 0; j < order; j++) {
        	AI[i * order + j] = 0;
        	for (k = 0; k < order; k++) {
            	AI[i * order + j] += (G[k * order + i] * G[k * order + j]);
        	}
    	}
	}
}

/**
 * @brief Computes in \b G the inverse of a lower triangular matrix \b L.
 * \n From book \e Digital \e Media \e Processing: \e DPS \e Algorithms \e Using
 * \e C. (H. Malepati, 2010) Pcode 9.6.
 */
void inverseTriangular(fpoint_t * L, fpoint_t * G) {
	int i, j, k;
	fpoint_t sum;

	for (i = 0; i < order; i++) {
    	G[i * (order + 1)] = ((fpoint_t) 1) / L[i * (order + 1)];
    	for (j = 0; j < i; j++) {
        	sum = 0;
        	for (k = j; k < i; k++) {
            	sum = sum - L[i * order + k] * G[k * order + j];
        	}
        	G[i * order + j] = sum / L[i * (order + 1)];
    	}
	}

}

/**
 * @brief Computes the mean error of matrix \b R in relation to matrix \b A.
 * \note It is up to the user to ensure that matrices \b A and \b R have the
 * same order.
 */
fpoint_t computeError(fpoint_t * A, fpoint_t * R) {
	int i, j;
	fpoint_t error = 0, sum = 0;

	for (i = 0; i < order; i++) {
    	for (j = 0; j < order; j++) {
        	error = A[i * order + j] - R[i * order + j];
        	sum += error/(order * order);
    	}
	}
	return sum;
}

/**
 * @brief Calculates a transpose matrix.
 * @param[in] A square matrix of ORDER #ORDER to be transposed.
 * @param[out] AT transpose matrix of "A".
 */
void transpose(fpoint_t * A, fpoint_t * AT) {
	int i, j;

	for (i = 0; i < order; i++) {
    	for (j = 0; j < order; j++) {
        	AT[i * order + j] = A[j * order + i];
    	}
	}
}

/**
 * @brief Computes the mean error of n-square matrix \b A in relation to Identity
 * matrix or order \c n.
 */
fpoint_t computePrecision(fpoint_t * A) {
	int i, j;
	fpoint_t error = 0, sum = 0;

	for (i = 0; i < order; i++) {
    	for (j = 0; j < order; j++) {
        	error = (((i == j) ? (fpoint_t) 1 : (fpoint_t) 0) - A[i * order + j]);
        	sum += error/(order * order);
    	}
	}
	return sum;
}

/*!
 * @brief Main function. For automated tests, two parameters must be given:
 * \n 1) Name of file in which the matrix is stored;
 * \n 2) Number of iteractions used in sqrt_ method.
 */
int main(int argc, char* argv[]) {

#if (TIMED == 1)
	clock_t start, end;
#endif
	int i;
	fpoint_t averageError;
	fpoint_t A[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; /*! \brief Matrix to be inverted */
	fpoint_t L[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; /*! \brief Lower triangular matrix obtained from Cholesky
               	* decomposition
               	*/
	fpoint_t G[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; /*! \brief Inverse of lower triangular matrix */

	// If tests are automated
	if (argc > 1) {
    	FILE *fid = fopen(argv[1], "r");
    	fscanf(fid, "%d", &order);
    	for (i = 0; i < (order * order); i++) {
        	fscanf(fid, "%lf", &A[i]);
    	}
    	fclose(fid);
    	n_int = atoi(argv[2]);
	}//otherwise, read from prompt.
	else {
    	FILE *fid = fopen("Matriz_5x5.txt", "r");
    	fscanf(fid, "%d", &order);
    	for (i = 0; i < (order * order); i++) {
        	fscanf(fid, "%lf", &A[i]);
    	}
    	fclose(fid);
    	n_int = 10;
	}


#if(PRINT_MATRICES == 1)
	// Print original matrix in screen
	printf("\nOriginal Matrix:\n");
	show_matrix(A);
#endif

#if (TIMED == 1)
	start = clock();
#endif
	// Computes lower trangular matrix from Cholesky Decomposition
	cholesky(A, L);

#if(PRINT_MATRICES == 1)
	// Print lower triangular matrix in screen
	printf("\nLower Triangular Matrix from Cholesky Decomposition:\n");
	show_matrix(L);
#endif

#if(ERROR_ANALYSIS > 1)
	// Computes mean error in recovering original matrix
	multiplyByTranspose(L, G);
	averageError = computeError(A, G);
	printf("\nAverage error in original matrix: %.9f\n", averageError);
	allocate(G);
#endif

	// Computes the inverse of lower triangular matrix and find inverse of A
	inverseTriangular(L, G); 	// G = L ^ (-1)
	allocate(L);
	inverseFromInvertedL(L, G); // L = G' * G

#if (TIMED == 1)
	end = clock();
	printf("\nTime: %f ms.\n", (((double) (end - start)) * 1000 / (CLOCKS_PER_SEC)));
#endif

#if(PRINT_MATRICES == 1)
	// Print inverse matrix in screen
	printf("\nInverse Matrix:\n");
	show_matrix(L);
#endif

	// Computes the Identity matrix from original matrix and its inverse
	allocate(G);
	multiply(A, L, G);

#if(PRINT_MATRICES == 1)
	// Print identity matrix in screen
	printf("\nIdentity Matrix:\n");
	show_matrix(G);
#endif

#if(ERROR_ANALYSIS > 0)
	// Computes mean error of multiplication
	averageError = computePrecision(G);
	printf("\nAverage error in identity matrix: %.9f\n", averageError);
#endif
	/* Deallocate matrices. */

	return 0;
}
