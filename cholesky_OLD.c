/**
 * @file main.c
 * @author  Heitor Freitas <heitor.freitas@gmail.com>
 * @author  Mariana Rodrigues <rodrigues.mariana@gmail.com>
 * @author  Thadeu A. F. Melo <thadeu.afm@gmail.com>
 * @version 1.0
 *
 * @license GPL
 *
 * @brief Matrix inversion throught Cholesky Decomposition.
 */


/**
 ******************
 * References:
 ******************
 *
 * Cholesky Decomposition:
 * http://www.icmc.usp.br/~andretta/ensino/aulas/sme0301-1-10/SistemasLinearesCholesky.pdf
 * http://wwwp.fc.unesp.br/~arbalbo/Iniciacao_Cientifica/sistemaslineares/teoria/3_Met_Cholesky.pdf
 * https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
 * http://rosettacode.org/wiki/Cholesky_decomposition#C
 * http://files.calcnum.webnode.com/200000021-605b36154f/Apostila%20-%20Cuminato.pdf
 *
 * Matrices:
 * http://rosettacode.org/wiki/Matrix_multiplication#C
 * http://www.programming-techniques.com/2011/09/numerical-methods-inverse-of-nxn-matrix.html
 *
 * GNU GCC Profiling Tool:
 * http://www.thegeekstuff.com/2012/08/gprof-tutorial/
 *
 * GCC optimization:
 * ftp://gcc.gnu.org/pub/gcc/summit/2003/Optimizing%20for%20space.pdf
 */



#include <stdio.h>
#include <stdlib.h>

typedef float fpoint_t;  /*! @typedef fpoint_t @brief Floating-point variable type. (typedef created in order to easily interchange between \c float and \c double 
types.*/

/*! Debug messages macro  */
#define DBG(x) printf("[%s:%d]%s():%s = %f\n",__FILE__,__LINE__,__FUNCTION__,#x,x);

int n_int;
/**
 * @brief Print n-order square matrix \b \c A in screen
 */
void show_matrix(fpoint_t * A, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%.9f ", A[i * n + j]);
        }
        printf("\n");
    }
}

/**
 * @brief Compute the square root of floating-point \b \c x using \b \c n_int iterations;
 */
fpoint_t sqrt_(fpoint_t x) {
    fpoint_t val = 1.0;
    while (n_int-- > 0) {
        val = val / 2 + x / (val * 2);
    }
    return val;
}

/**
 * @brief Computes lower trianguar matrix \b \c L component in Cholesky Decomposition of n-order square matrix \b \c A. Code based on \url{http://rosettacode.org/wiki/Cholesky_decomposition#C}
 */
void cholesky(fpoint_t * A, fpoint_t * L, int n) {
    int i, j, k;
    fpoint_t s;

    for (i = 0; i < n; i++) {
        for (j = 0; j < (i + 1); j++) {
            s = 0;
            for (k = 0; k < j; k++) {
                s += L[i * n + k] * L[j * n + k];
            }
            L[i * n + j] = (i == j) ? sqrt_(A[i * n + i] - s) : (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }
    }
}

/**
 * @brief Computes the dot product between a row of n-order square matrix \b \c A and a column of n-order square matrix \b \c B. Code based on \url{http://rosettacode.org/wiki/Matrix_multiplication#C}
 */
fpoint_t dot(fpoint_t *A, fpoint_t *B, int n) {
    fpoint_t r = 0;
    int len = n;
    while (len--) {
        r += *A++ * *B;
        B += n;
    }
    return r;
}


/**
 * @brief Stores in matrix \b \c C the result of multiplication between n-order square matrixes \b \c A and \b \c B. Code based on \url
{http://rosettacode.org/wiki/Matrix_multiplication#C}
 */
void multiply(fpoint_t * A, fpoint_t * B, fpoint_t * C, int n) {
    fpoint_t * pa;
    int i, j;

    for (pa = A, i = 0; i < n; i++, pa += n) {
        for (j = 0; j < n; j++) {
            *C++ = dot(pa,B + j,n);
        }
    }
}

/**
 * @brief Stores in matrix \b \c R the multiplication between n-order square matrix \b \c A and its transpose. Code based on \url
{http://rosettacode.org/wiki/Matrix_multiplication#C}
 */
void multiplyByTranspose(fpoint_t * A, fpoint_t * B, int n) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        for(j = 0; j < n; j++){
            B[i * n + j] = 0;
            for(k = 0; k < n; k++){
                B[i*n + j] += (A[i*n+k] * A[j*n+k]);
            }
        }
    }
}

/**
 * @brief Computes in \b \c LI the inverse of n-order lower triangular matrix \b \c L. Code based on Pcode 9.6 from book \e Digital \e Media \e Processing: \e DPS \e 
Algorithms \e Using \e C. (H. Malepati, 2010).
 */
void inverseTriangular(fpoint_t * L, fpoint_t * LI,int n) {
    int i, j, k;
    fpoint_t sum;

    for (i = 0; i < n; i++) {
        LI[i * (n + 1)] = ((fpoint_t) 1) / L[i * (n + 1)];
        for (j = 0; j < i; j++) {
            sum = 0;
            for (k = j; k < i; k++) {
                sum = sum - L[i * n + k] * LI[k * n + j];
            }
            LI[i * n + j] = sum / L[i * (n + 1)];
        }
    }

}

/**
 * @brief Computes the average error of n-square matrix \b \c R in relation to n-order square matrix \b \c A. 
 */
fpoint_t computeError(fpoint_t * A, fpoint_t * R, int n){
    int i,j;
    fpoint_t error = 0;

    for(i = 0; i< n; i++){
        for(j=0;j<n;j++){
            error += A[i*n + j] - R[i*n +j];
        }
    }
    return (erro/(n*n));
}

/**
 * @brief Computes in n-order square matrix \b \c AT the transpose of n-order square matrix \b \c A.
 */
void transpose(fpoint_t * A, fpoint_t * AT, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            AT[i * n + j] = A[j * n + i];
        }
    }
}

/**
 * @brief Computes the average error of n-square matrix \b \c R in relation to n-order Identity matrix.
 */
fpoint_t computePrecision(fpoint_t * R, int n){
    int i,j;
    fpoint_t error = 0;

    for(i = 0; i< n; i++){
        for(j=0;j<n;j++){
            error += ((i==j)? (fpoint_t) 1: (fpoint_t) 0)- R[i*n +j];
        }
    }
    return (error/(n*n));
}

int main(void) {

    int i, j;
    fpoint_t averageError;
    n_int = 5;		// \TODO: Insert as parameter
    int order = 3;	// \TODO: Insert as parameter
    
    fpoint_t A[] = {25, 15, -5,
             15, 18,  0,
             -5,  0, 11};
    fpoint_t * L = (fpoint_t*) calloc(order * order, sizeof (fpoint_t)); //Matriz triangular inferior resultante da decomposição de Cholesky (L)
    fpoint_t *LI = (fpoint_t*) calloc(order * order, sizeof (fpoint_t)); //Matriz inversa da matriz triangular inferior (L^-1)
    fpoint_t *LIT = (fpoint_t*) calloc(order * order, sizeof (fpoint_t)); //Matriz transposta da matriz inversa da matriz triangular inferior ((L^-1)^T)
    fpoint_t *AI = (fpoint_t*) calloc(order * order, sizeof (fpoint_t)); //Matriz inversa de A (A^-1)
    fpoint_t *I = (fpoint_t*) calloc(order * order, sizeof (fpoint_t)); //Matriz identidade (I)


    //Mostra a matriz original
    printf("\nMatriz Original (A):\n");
    show_matrix(A,order);

    printf("\nMatriz triangular inferior por Cholesky (L):\n");
    cholesky(A, L,order);
    show_matrix(L,order);

    //Calcula Transposta
    transpose(L, LIT,order);

    printf("\nRecupera Matriz A:\n");
    multiply(L, LIT, I,order);
    show_matrix(I,order);

    printf("\nErro de recuperacao:\n");
    averageError = calculaErro(A,I,AI,order);
    show_matrix(AI,order);
    printf("\nErro medio: %.9f\n",averageError);

    //Calcula a matriz inversa da matriz triangular inferior
    printf("\nLI:\n");
    inverseTriangular(L, LI,order);
    show_matrix(LI,order);

    //Calcula a matriz transposta da matriz inversa da matriz triangular inferior.
    /**
     * Código LaTeX:
     * \begin{equation}
     *      (L^{-1})^{T}
     * \end{equation}
     */
    printf("\n(LI)T:\n");
    transpose(LI, LIT,order);
    show_matrix(LIT,order);

    //Calcula a matriz inversa, ie: AI = (LI)T * (LI) */
    /**
     * Código LaTeX:
     * \begin{equation}
     *      A^{-1} = (L^{-1})^T(L^{-1})
     * \end{equation}
     */
    printf("\nMatriz inversa (A^-1):\n");
    multiply(LIT, LI, AI,order);
    show_matrix(AI,order);

    printf("\nMatriz inversa nova rotina (A^-1):\n");
    multiplyByTranspose(LIT,AI,order);
    show_matrix(AI,order);
    

    //Calcula a matriz identidade, ie: I = A*AI
    /**
     * Código LaTeX:
     * \begin{equation}
     *      I = AA^{-1}
     * \end{equation}
     */

    printf("\nMatriz identidade (A*A^-1):\n");
    multiply(A, AI, I,order);
    show_matrix(I,order);

    printf("\nPrecision Error:\n");
    averageError = computePrecision(I,order);
    printf("\nAverage: %.9f\n",averageError);

    /* Desaloca as matrizes da RAM. */
    free(L);
    free(LI);
    free(LIT);
    free(AI);
    free(I);

    return 0;
}
