/*
 *
 * matrix.cpp
 * Hetero2
 *
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2014, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * All rights reserved. CSIRO is willing to grant you a license to this Hetero Version 2 on the following terms,
 * except where otherwise indicated for third party material.
 * Redistribution and use of this software in source and binary forms, with or without modification, are permitted
 * provided that the following conditions are met:
 * •	Redistributions of source code must retain the above copyright notice, this list of conditions and the
 *      following disclaimer.
 * •	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
 *      following disclaimer in the documentation and/or other materials provided with the distribution.
 * •	Neither the name of CSIRO nor the names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission of CSIRO.
 * EXCEPT AS EXPRESSLY STATED IN THIS AGREEMENT AND TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, THE SOFTWARE
 * IS PROVIDED "AS-IS". CSIRO MAKES NO REPRESENTATIONS, WARRANTIES OR CONDITIONS OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO ANY REPRESENTATIONS, WARRANTIES OR CONDITIONS REGARDING THE CONTENTS OR ACCURACY
 * OF THE SOFTWARE, OR OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, THE ABSENCE
 * OF LATENT OR OTHER DEFECTS, OR THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT DISCOVERABLE.
 * TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT SHALL CSIRO BE LIABLE ON ANY LEGAL THEORY (INCLUDING,
 * WITHOUT LIMITATION, IN AN ACTION FOR BREACH OF CONTRACT, NEGLIGENCE OR OTHERWISE) FOR ANY CLAIM, LOSS, DAMAGES OR
 * OTHER LIABILITY HOWSOEVER INCURRED.  WITHOUT LIMITING THE SCOPE OF THE PREVIOUS SENTENCE THE EXCLUSION OF LIABILITY
 * SHALL INCLUDE: LOSS OF PRODUCTION OR OPERATION TIME, LOSS, DAMAGE OR CORRUPTION OF DATA OR RECORDS; OR LOSS OF
 * ANTICIPATED SAVINGS, OPPORTUNITY, REVENUE, PROFIT OR GOODWILL, OR OTHER ECONOMIC LOSS; OR ANY SPECIAL, INCIDENTAL,
 * INDIRECT, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES, ARISING OUT OF OR IN CONNECTION WITH THIS AGREEMENT, ACCESS
 * OF THE SOFTWARE OR ANY OTHER DEALINGS WITH THE SOFTWARE, EVEN IF CSIRO HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
 * CLAIM, LOSS, DAMAGES OR OTHER LIABILITY.
 * APPLICABLE LEGISLATION SUCH AS THE AUSTRALIAN CONSUMER LAW MAY APPLY REPRESENTATIONS, WARRANTIES, OR CONDITIONS, OR
 * IMPOSES OBLIGATIONS OR LIABILITY ON CSIRO THAT CANNOT BE EXCLUDED, RESTRICTED OR MODIFIED TO THE FULL EXTENT SET OUT
 * IN THE EXPRESS TERMS OF THIS CLAUSE ABOVE "CONSUMER GUARANTEES".  TO THE EXTENT THAT SUCH CONSUMER GUARANTEES
 * CONTINUE TO APPLY, THEN TO THE FULL EXTENT PERMITTED BY THE APPLICABLE LEGISLATION, THE LIABILITY OF CSIRO UNDER
 * THE RELEVANT CONSUMER GUARANTEE IS LIMITED (WHERE PERMITTED AT CSIRO’S OPTION) TO ONE OF FOLLOWING REMEDIES OR
 * SUBSTANTIALLY EQUIVALENT REMEDIES:
 * (a) THE REPLACEMENT OF THE SOFTWARE, THE SUPPLY OF EQUIVALENT SOFTWARE, OR SUPPLYING RELEVANT SERVICES AGAIN;
 * (b) THE REPAIR OF THE SOFTWARE;
 * (c) THE PAYMENT OF THE COST OF REPLACING THE SOFTWARE, OF ACQUIRING EQUIVALENT SOFTWARE, HAVING THE RELEVANT
 *     SERVICES SUPPLIED AGAIN, OR HAVING THE SOFTWARE REPAIRED.
 * IN THIS CLAUSE, CSIRO INCLUDES ANY THIRD PARTY AUTHOR OR OWNER OF ANY PART OF THE SOFTWARE OR MATERIAL DISTRIBUTED
 * WITH IT.  CSIRO MAY ENFORCE ANY RIGHTS ON BEHALF OF THE RELEVANT THIRD PARTY.
 * Third Party Components
 * The following third party components are distributed with the Software.  You agree to comply with the license terms
 * for these components as part of accessing the Software.  Other third party software may also be identified in
 * separate files distributed with the Software.
 * ___________________________________________________________________
 * 
 * JACOBI_EIGENVALUE.C (http://people.sc.fsu.edu/~jburkardt/c_src/jacobi_eigenvalue/jacobi_eigenvalue.c)
 * Copyright (C) 2003-2013 John Burkardt
 * This software is licensed under GNU LGPL (http://www.gnu.org/licenses/lgpl.html)
 * ___________________________________________________________________
 */

#include "matrix.h"

// construct a matrix
// the values inside are not resetted
double* matrix(int row, int col) {
    return (double*) malloc (sizeof(double)*row*col);
}

// reset the matrix to zero
void resetMat(double* mat, int row, int col) {
    memset(mat,0,sizeof(double)*row*col);
}

// reset a square matrix to an identity matrix
void ident(double* sqMat, int dim) {
    memset(sqMat,0,sizeof(double)*dim*dim);
    for (int i=0; i<dim; i++)
        sqMat[i*dim+i] = 1.0;
}

// set the values of the particular row of the matrix same as inArray
void setMatrix(double* mat, int col, int row, double inArray[]) {
    double* fr_ptr = &(inArray[0]);
    double* to_ptr = &(mat[row*col]);
    memcpy(to_ptr,fr_ptr,col*sizeof(double));
}

// set all entries inside the matrix to a specific value
void setMatrixToVal(double* mat, int row, int col, double value) {
    for (int i=0; i<row*col; i++) {
        mat[i] = value;
    }
}

// to set the diagional of the square matrix same as inArray
void setDiag(double* sqMatrix, double* inArray, int dim) {
    for (int i=0; i<dim; i++)
        sqMatrix[i*dim+i]=inArray[i];
}

// output a diagional matrix from an array
double* toDiagMat(double* inArray, int dim) {
    double* result = matrix(dim,dim);
    resetMat(result, dim, dim);
    setDiag(result, inArray, dim);
    return result;
}

// to collect the diagional of a square matrix
void diag(double* outArray, double* inSqMat, int dim) {
    for (int i=0; i<dim; i++)
        outArray[i] = inSqMat[i*dim+i];
}
double* diag(double* inSqMat, int dim) {
    double* result = matrix(1, dim);
    diag(result,inSqMat,dim);
    return result;
}

// duplicate
void duplicate(double* toMat, double* frMat, int row, int col) {
    memcpy(toMat,frMat,row*col*sizeof(double));
}
double* duplicate(double* frMat, int row, int col) {
    double* toMat = matrix(row, col);
    memcpy(toMat,frMat,row*col*sizeof(double));
    return toMat;
}

// print out a matrix
void printMatrix(double* matrix, int row, int col){
    for (int i=0; i<row; i++) {
        for (int j=0; j<col; j++) {
            cout << matrix[i*col+j] << "\t";
        }
        cout << endl;
    }
}

void jacobi(double* sqMatrix, int dim, double* eigenVals, double* eigenVects, int& numRound)
{
    int rot_num;
    jacobi_eigenvalue ( dim, sqMatrix, MAX_ROUND, eigenVects, eigenVals, &numRound, &rot_num );
    
    // transpose eigenVects
    double tmp;
    int i, j;
    for (i=0; i<dim; i++) {
        for (j=i+1; j<dim; j++) {
            tmp = eigenVects[i*dim+j];
            eigenVects[i*dim+j] = eigenVects[j*dim+i];
            eigenVects[j*dim+i] = tmp;
        }
    }
}

// to compute the transpose of a matrix with dimension row x col
void transpose(double* outMat, double* inMat, int row, int col) {
    for (int i=0; i<row; i++) {
        for (int j=0; j<col; j++) {
            outMat[j*row+i] = inMat[i*col+j];
        }
    }
}
double* transpose(double* inMat, int row, int col) {
    double* result = matrix(col,row);
    transpose(result, inMat, row, col);
    return result;
}

// to compute the multiplication of two matrice
// dimension of matrix 1: row x common_len
// dimension of matrix 2: common_len x col
void multiply(double *outMat, double* inMat1, double* inMat2, int row, int common_len, int col){
    // reset the outSqMat to zero
    memset(outMat,0,sizeof(double)*row*col);
    for (int i=0; i<row; i++)
        for (int j=0; j<col; j++)
            for (int k=0; k<common_len; k++)
                outMat[i*col+j] += inMat1[i*common_len+k]*inMat2[k*col+j];
}
double* multiply(double * inMat1, double* inMat2, int row, int common_len, int col) {
    double* result = matrix(row,col);
    multiply(result,inMat1,inMat2,row,common_len,col);
    return result;
}
