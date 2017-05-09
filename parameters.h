/*
 *
 * parameters.h
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

#ifndef __RAL_RAS__parameters__
#define __RAL_RAS__parameters__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "matrix.h"
#include "tool_box.h"
#include "definitions.h"

class VariableSet {
public:
    double beta;
    double* alpha;
    double* probXGivenInv;
    int num_alpha;
    
    double* rootNodeFreq;
    
    VariableSet(int num_alpha); // constructor
    ~VariableSet(); //destructor
    
    void readSiteInfoFile(char* fileName, vector<string>* siteCatNames);
    // read the values of the variables from the site info file
    
    // copy the whole content (except rateIDs) from another instance of ParameterSet
    void copyFrom(VariableSet& vs);
    
    void showContent();
    void showContent(string& outStr);
    // show the content
    void showContent(ofstream& fout);
};

class ParameterSet {
    
public:
    double* w;  // variables in the S vector for all edges
    double* pi; // stationary probabilities for all edges
    double* t;  // time for all edges
    int num_w; // number of variables in S vector
    int num_edges; // number of edges
    
    int num_rate_matrices; // number of rate matrices
    vector<vector<int>* > rateIDs; // the ID of rate matrices for each edge
    
    double* allEigenMat;  // the eigen matrix for all edges
    // i.e. sqrt(pi) * R * inverse of sqrt(pi), where R = S * Pi
    double* allCondProb; // conditional probability matrices for all edges
    
    ParameterSet(); // constructor
    ~ParameterSet(); // destructor
    
    // read the values of w and pi from the parameter file
    void readParamFile(char* fileName, vector<string>& nodeList, vector<double>& edgeLens,
                        int catID, int num_edges, int edgeRepresent);

    // load the file of rate groups
    // prerequisite: the value of "num_edges" has to be set
    // return the number of rate groups
    void loadRateMat(char* rateGrpFile);
    
    void initialize(int num_w, int num_edges); // initialize the parameters
    
    // update the parameters such that
	// outFormat = 1 : Edge length is set to the rate of substitution; OR
	//             2 : S6 in the rate matrix is set to 1
	void updateContent(int outFormat);

    void showContent(string &outStr, int* topMatrix, vector<string>* leafList);
    void showContent();
    void showContent(ofstream& fout);
    
    // compute the eigen-matrix for all edges
    // output: eigen_matrix
    void computeAllEigenMatrix();
    
    // OBJECTIVE:             compute the conditional probabilities for all the edges.
    // OUTPUT:                allCondProb
    // SIMILAR FUNCTION IN R: getCondProb
    // PREREQUISITE:          computeAllEigenMatrix()
    void computeAllCondProb(int isReversible);
    
    void printAllCondProb();
    void printAllCondProb2();
    
    void print_t();
    
    // output the rate matrix in single array form
    // 1-based
    void outputRateMat(vector<int>& outRateMatArray);
    
    // update the rate matrix
    void updateRateMat(vector<int>& rateMatArray, int numRateGrp);
    
    void printRateMatID();
    
    // copy the whole content (except rateIDs) from another instance of ParameterSet
    void copyFrom(ParameterSet& ps);
    
    
private:
    // Temporary memory allocation
    // for internal computation use only
    double* s;
    double* eigenVals;
    double* eigenVects;
    double* tmpMat1;
    double* tmpMat2;
    double* sqrtPi;
    
    // get the initialize ID of rate matrices for each edge
    // before: the values of we have been read from the file
    void initializeRateMatID();
};


class AllParameterSet {
public:
    vector<ParameterSet*> ps;
    
    vector<double*> allCondProbSet;
    
    int numRateCat;
    
    // constructor
    AllParameterSet(int numRateCat);
    
    // destructor
    ~AllParameterSet();
    
    // read the values of w and pi from the parameter file list
    void readParamFileList(char* fileName, vector<string>& categoryList, vector<string>& nodeList,
                           vector<double>& edgeLens, int num_edges, int edgeRepresent);
                                        
    // load the file of rate groups
    // return the number of rate groups
    int loadRateMat(char* rateGrpFile);

    // update the parameters such that
	// outFormat = 1 : Edge length is set to the rate of substitution; OR
	//             2 : S6 in the rate matrix is set to 1
	void updateContent(int outFormat);
    
    void showContent();
    void showContent(string& outStr, int* topMatrix, vector<string>* leafList);

    // compute the eigen-matrix for all edges for all parameter sets
    // output: eigen_matrix
    void computeAllEigenMatrix();
    
    // OBJECTIVE:             compute the conditional probabilities for all the edges for all parameter sets
    // OUTPUT:                allCondProb
    // SIMILAR FUNCTION IN R: getCondProb
    // PREREQUISITE:          computeAllEigenMatrix()
    void computeAllCondProb(int isReversible);

    void printAllCondProb();
    void printAllCondProb2();
    
    int size();
};


// obtain the s-matrix
void computeSMatrix(double* s_matrix, double* w, double* pi);

// compute the eigen-matrix
void computeEigenMatrix(double* eigenMat, double* w, double* pi, double* s);

// print all conditional probabilities for all rate categories
void printAllCatCondProb(vector<ParameterSet*>& ps_set);

#endif /* defined(__RAL_RAS__parameters__) */
