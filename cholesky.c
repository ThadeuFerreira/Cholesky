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


int ORDER=0;                 /*!< Order of square matrix. */
#define SQRT_INT 5 /*! @def SQRT_INT @brief Iteration number of routine sqrt_ */



/******************************************************************************/
/*                               LIBRARIES                                     */
/******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
typedef float fpoint_t; /*! @typedef fpoint_t @brief Type of variable for float point \c float and \c fpoint_t */
/**
 * Macro for debug purposes.
 */
#define DBG(x) printf("[%s:%d]%s():%s = %f\n",__FILE__,__LINE__,__FUNCTION__,#x,x);

/**
 * @brief Show the matrix \b A 
 * @param[in] A matrix to be showed
 */
void show_matrix(fpoint_t * A) {
    int i, j;
    for (i = 0; i < ORDER; i++) {
        for (j = 0; j < ORDER; j++) {
            printf("%.9f ", A[i * ORDER + j]);
        }
        printf("\n");
    }
}

/**
 * @brief Calculates the square root of a number "x".
 * @param[in] x number to be calculated square root.
 * @return Square root of a number "x".
 */
fpoint_t sqrt_(fpoint_t x) {
    fpoint_t val = 1.0;
    int times = SQRT_INT;
    while (times-- > 0) {
        val = val / 2 + x / (val * 2);
    }
    return val;
}

/**
 * @brief Do the Cholesky decomposition. Modified of http://rosettacode.org/wiki/Cholesky_decomposition#C
 * @param[in] A matrix of order #ORDER.
 * @param[out] L lower triangular matrix from of Cholesky decomposition.
 */
void cholesky(fpoint_t * A, fpoint_t * L) {
    int i, j, k;
    fpoint_t s;

    for (i = 0; i < ORDER; i++) {
        for (j = 0; j < (i + 1); j++) {
            s = 0;
            for (k = 0; k < j; k++) {
                s += L[i * ORDER + k] * L[j * ORDER + k];
            }
            L[i * ORDER + j] = (i == j) ? sqrt_(A[i * ORDER + i] - s) : (1.0 / L[j * ORDER + j] * (A[i * ORDER + j] - s));
        }
    }
}

/**
 * @brief Do a matrix elements multiplication.
 * Modified of http://rosettacode.org/wiki/Matrix_multiplication#C
 * @param[in] A pointer to matrix.
 * @param[in] B pointer to matrix.
 * @return multiplication of two matrices elements.
 */
fpoint_t dot(fpoint_t *A, fpoint_t *B) {
    fpoint_t r = 0;
    int len = ORDER;
    while (len--) {
        r += *A++ * *B;
        B += ORDER;
    }
    return r;
}

/**
 * @brief Do a matrix elements multiplication.
 * Modified of http://rosettacode.org/wiki/Matrix_multiplication#C
 * @param[in] A pointer to matrix.
 * @param[in] B pointer to matrix.
 * @return multiplication of two matrices elements.
 */
fpoint_t dot2(fpoint_t *A, fpoint_t *B) {
    fpoint_t r = 0;
    int len = ORDER;
    
    while (len--) {
        r += *A * *B;
        A += ORDER;
        B += ORDER;
    }
    return r;
}

/**
 * @brief Do a multiplication of two matrices. Modified of http://rosettacode.org/wiki/Matrix_multiplication#C
 * @param[in] A first matrix.
 * @param[in] B second matrix.
 * @param[out] C result of the multiplication.
 */
void multiply(fpoint_t * A, fpoint_t * B, fpoint_t * C) {
    fpoint_t * pa;
    int i, j;

    for (pa = A, i = 0; i < ORDER; i++, pa += ORDER) {
        for (j = 0; j < ORDER; j++) {
            *C++ = dot(pa, B + j);
        }
    }
}

/**
 * @brief Do a multiplication of a lower triangular with its transpose. 
 * Note they are the same elements, so only one input matrix is need 
 * Adapted from http://rosettacode.org/wiki/Matrix_multiplication#C
 * @param[in] A matrix.
 * @param[out] C result of the multiplication.
 */
void multiplyByTranspose(fpoint_t * A, fpoint_t * C) {
    fpoint_t *pa;
    int i, j;

    for (pa = A, i = 0; i < ORDER; i++, pa += ORDER) {
        for (j = 0; j < ORDER; j++) {
            *C++ = dot2(A+i, A + j);
        }
    }
}

/**
 * @brief Computes the inverse of a lower triangular matrix L.
 * From book \e Digital \e Media \e Processing: \e DPS \e Algorithms \e Using \e C. (H. Malepati, 2010) Pcode 9.6.
 * @param[in] L lower inverse triangular matrix of order #ORDER to be inverted.
 * @param[out] G inverse matrix \b L
 */
void inverseTriangular(fpoint_t * A, fpoint_t * AI) {
    int i, j, k;
    fpoint_t sum;

    for (i = 0; i < ORDER; i++) {
        AI[i * (ORDER + 1)] = ((fpoint_t) 1) / A[i * (ORDER + 1)];
        for (j = 0; j < i; j++) {
            sum = 0;
            for (k = j; k < i; k++) {
                sum = sum - A[i * ORDER + k] * AI[k * ORDER + j];
            }
            AI[i * ORDER + j] = sum / A[i * (ORDER + 1)];
        }
    }

}

/**
 * @brief Calculates the error matrix.
 * @param[in] A square matrix of #ORDER original order.
 * @param[in] R square matrix of #ORDER recovered order.
 * @param[out] E error matrix with the difference between the elements.
 * @return mean square error.
 */
fpoint_t calculaErro(fpoint_t * A, fpoint_t * R, fpoint_t * E){
    int i,j;
    fpoint_t erro = 0;

    for(i = 0; i< ORDER; i++){
        for(j=0;j<ORDER;j++){
            E[i * ORDER + j] = A[i * ORDER + j] - R[i * ORDER + j];
            erro += E[i*ORDER + j] * E[i*ORDER +j];
        }
    }
    return (sqrt_(erro)/(ORDER*ORDER));
}

/**
 * @brief Calculates a transpose matrix.
 * @param[in] A square matrix of ORDER #ORDER to be transposed.
 * @param[out] AT transpose matrix of "A".
 */
void transpose(fpoint_t * A, fpoint_t * AT) {
    int i, j;

    for (i = 0; i < ORDER; i++) {
        for (j = i; j < ORDER; j++) {
            AT[i * ORDER + j] = A[j * ORDER + i];
        }
    }
}

/**
 * @brief Computes the precision of ...
 * @para[in] A matrix ...
 * @return the precision of...
 */
fpoint_t computePrecision(fpoint_t * A){
    int i,j;
    fpoint_t erro = 0, soma = 0;

    for(i = 0; i< ORDER; i++){
        for(j=0;j<ORDER;j++){
            erro = ((i==j)? (fpoint_t) 1: (fpoint_t) 0)- A[i*ORDER +j];
            soma = soma + erro * erro;
        }
    }
    return (sqrt_(soma)/(ORDER*ORDER));
}

/**
 * @brief
 * @return
 */
int main(int argc, char* argv[]) {

    /*
     * Faz a alocação dinâmica das matrizes.
     * Algumas pode ser "reaproveitadas" para economizar memória.
     */
    int i;
    fpoint_t averageError;
    
	// Matriz a ser invertida
    fpoint_t *A;
    
    //Caso dados estejam em arquivo
    if(argc > 1){
	    FILE *fid = fopen(argv[1], "r");
	    fscanf(fid, "%d", &ORDER);
	    A = (fpoint_t*) calloc(ORDER*ORDER, sizeof(fpoint_t));
	    for (i=0; i < (ORDER * ORDER); i++){
    		fscanf(fid, "%f", &A[i]);
    	}
    	fclose(fid);
    }
    //Caso os dados sejam digitados via linha de comando
    else {
    	scanf("%d", &ORDER);
    	A = (fpoint_t*) calloc(ORDER*ORDER, sizeof(fpoint_t));
    	for (i=0; i < (ORDER * ORDER); i++){
    		scanf("%f", &A[i]);
    	}
    }
    
             
    
    fpoint_t * L = (fpoint_t*) calloc(ORDER * ORDER, sizeof (fpoint_t)); //Matriz triangular inferior resultante da decomposição de Cholesky (L)
    fpoint_t *LI = (fpoint_t*) calloc(ORDER * ORDER, sizeof (fpoint_t)); //Matriz inversa da matriz triangular inferior (L^-1)
    fpoint_t *AUX = (fpoint_t*) calloc(ORDER * ORDER, sizeof (fpoint_t)); //Matriz auxiliar
    fpoint_t *AI = (fpoint_t*) calloc(ORDER * ORDER, sizeof (fpoint_t)); //Matriz inversa de A (A^-1)
    fpoint_t *I = (fpoint_t*) calloc(ORDER * ORDER, sizeof (fpoint_t)); //Matriz identidade (I)


    //Mostra a matriz original
    printf("\nMatriz Original (A):\n");
    show_matrix(A);

    //Calcula a matriz triangular inferior pela decomposição de Cholesky
    printf("\nMatriz triangular inferior por Cholesky (L):\n");
    cholesky(A, L);
    show_matrix(L);

    //Calcula Transposta
    transpose(L, AUX);

    printf("\nRecupera Matriz A:\n");
    multiply(L, AUX, I);
    show_matrix(I);

    printf("\nErro de recuperacao:\n");
    averageError = calculaErro(A,I,AI);
    show_matrix(AI);
    printf("\nErro medio: %.9f\n",averageError);

    //Calcula a matriz inversa da matriz triangular inferior
    printf("\nLI:\n");
    inverseTriangular(L, LI);
    show_matrix(LI);

    //Calcula a matriz inversa, ie: AI = (LI)T * (LI) */
    /**
     * Código LaTeX:
     * \begin{equation}
     *      A^{-1} = (L^{-1})^T(L^{-1})
     * \end{equation}
     */
    printf("\nMatriz inversa nova rotina (A^-1):\n");
    multiplyByTranspose(LI,AI);
    show_matrix(AI);

    //Calcula a matriz identidade, ie: I = A*AI
    /**
     * Código LaTeX:
     * \begin{equation}
     *      I = AA^{-1}
     * \end{equation}
     */

    printf("\nMatriz identidade (A*A^-1):\n");
    multiply(A, AI, I);
    show_matrix(I);

    printf("\nPrecision Error:\n");
    averageError = computePrecision(I);
    printf("\nAverage: %.9f\n",averageError);

    /* Desaloca as matrizes da RAM. */
    free(L);
    free(LI);
    free(AUX);
    free(AI);
    free(I);

    return 0;
}
