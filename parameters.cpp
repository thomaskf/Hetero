/*
 *
 * parameters.cpp
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

#include "parameters.h"

VariableSet::VariableSet(int num_alpha) {
    // constructor
    
    this->num_alpha = num_alpha;
    beta = 0.0;
    alpha = (double*) malloc (sizeof(double)*num_alpha);
    memset(alpha,0,sizeof(double)*num_alpha);
    probXGivenInv = (double*) malloc (4 * sizeof(double));
    memset(probXGivenInv,0,4 * sizeof(double));
    rootNodeFreq = (double*) malloc (4 * sizeof(double)*num_alpha);
    memset(rootNodeFreq,0,4 * sizeof(double)*num_alpha);
}

VariableSet::~VariableSet() {
    // destructor
    if (alpha != NULL)
        free(alpha);
    if (probXGivenInv != NULL)
        free(probXGivenInv);
    if (rootNodeFreq != NULL)
        free(rootNodeFreq);
}


void VariableSet::readSiteInfoFile(char* fileName, vector<string>* siteCatNames) {
    // read the values of the variables from the site info file
    
    ifstream fin;
    fin.open(fileName);
    if (!fin.is_open()) {
        cerr << "Error opening the variable file :" << fileName << endl;
        exit(1);
    }
    string aline;
    vector<string> token;
    string catName;
    int isVariant;
    int catID;
    int* catLoaded = new int[num_alpha+1]; // catLoaded[num_alpha] is for the invariant site
    memset(catLoaded, 0, sizeof(int)*(num_alpha+1));
    while (getline(fin,aline)) {
        if ((aline.length() > 0) && (aline[0] != '#')) {
            tokenizer(aline," \t", &token);
            if (token.size() < 7)
                continue;
            catName = token[0];
            if (token[1] == "invariant") {
                isVariant = 0;
            } else if (token[1] == "variant") {
                isVariant = 1;
                // check for the category name
                for (catID=0; catID<siteCatNames->size(); catID++)
                    if (catName == siteCatNames->at(catID))
                        break;
                if (catID >= siteCatNames->size()) {
                    cerr << "Error! The variant category " << catName << " does not appear in the tree file" << endl;
                    exit(1);
                }
            } else {
                cerr << "Error! No indication on whether the category " << catName << " is variant or invariant in the \"site info file\"" << endl;
                exit(1);
            }
            if (isVariant) {
                if (catLoaded[catID]) {
                    cerr << "Error! The variant category " << catName << " appears twice in the \"site info file\"" << endl;
                    exit(1);
                }
                alpha[catID] = atof(token[2].c_str()); // proportion
                rootNodeFreq[catID*4  ] = atof(token[3].c_str()); // freq(A)
                rootNodeFreq[catID*4+1] = atof(token[4].c_str()); // freq(C)
                rootNodeFreq[catID*4+2] = atof(token[5].c_str()); // freq(G)
                rootNodeFreq[catID*4+3] = atof(token[6].c_str()); // freq(T)
                catLoaded[catID] = 1;
            } else {
                if (catLoaded[num_alpha]) {
                    cerr << "Error! There are more than on invariant category in the \"site info file\"" << endl;
                    exit(1);
                }
                beta = atof(token[2].c_str()); // proportion
                probXGivenInv[0] = atof(token[3].c_str()); // freq(A)
                probXGivenInv[1] = atof(token[4].c_str()); // freq(C)
                probXGivenInv[2] = atof(token[5].c_str()); // freq(G)
                probXGivenInv[3] = atof(token[6].c_str()); // freq(T)
                catLoaded[num_alpha] = 1;
            }
        }
    }
    fin.close();
    // if there is no constant site
    if (!catLoaded[num_alpha]) {
        beta = 0; // proportion
        probXGivenInv[0] = 0; // freq(A)
        probXGivenInv[1] = 0; // freq(C)
        probXGivenInv[2] = 0; // freq(G)
        probXGivenInv[3] = 0; // freq(T)
    }
}

/*
void VariableSet::resetAllVariables(Alignment& alignment) {
    // reset the variables according to the constant sites
    
    // set the value of beta to 0.75 * proportion of invariable sites
    beta = 0.75 * (double) alignment.numConstSites / (double) alignment.numSites;
    
    for (int i=0; i<num_alpha; i++) {
        alpha[i] = (1.0 - beta) / num_alpha;
        
        for (int k=0; k<4; k++)
            rootNodeFreq[i*4 + k] = 1.0 / 4;
    }
    
    for (int i=0; i<4; i++) {
        probXGivenInv[i] = 1.0 / 4;
    }
    
}
 */

// copy the whole content (except rateIDs) from another instance of ParameterSet
void VariableSet::copyFrom(VariableSet& vs) {
    beta = vs.beta;
    num_alpha = vs.num_alpha;
    if (alpha==NULL)
        alpha = (double*) malloc (sizeof(double)*num_alpha);
    memcpy(alpha,vs.alpha,sizeof(double)*num_alpha);
    if (probXGivenInv==NULL)
        probXGivenInv = (double*) malloc (4 * sizeof(double));
    memcpy(probXGivenInv,vs.probXGivenInv,4 * sizeof(double));
    if (rootNodeFreq==NULL)
        rootNodeFreq = (double*) malloc (4 * sizeof(double)*num_alpha);
    memcpy(rootNodeFreq,vs.rootNodeFreq,4 * sizeof(double)*num_alpha);
}


void VariableSet::showContent() {
    // show the content
    
    cout << "beta : " << beta << endl;
    cout << "alpha : ";
    for (int i=0; i<num_alpha; i++)
        cout << alpha[i] << ",";
    cout << endl;
    cout << "probXGivenInv : ";
    for (int i=0; i<4; i++)
        cout << probXGivenInv[i] << ",";
    cout << endl;
    cout << "rootNodeFreq : ";
    for (int i=0; i<4*num_alpha; i++)
        cout << rootNodeFreq[i] << ",";
    cout << endl;
}

void VariableSet::showContent(string& outStr) {
    // show the content
    
		outStr.append("[Nucleotide frequencies in root]\n");
		outStr.append("               \tfreq(A)\tfreq(C)\tfreq(G)\tfreq(T)\n");
		outStr.append("Invariable_site");
		for (int i=0; i<4; i++) {
            outStr.append("\t");
            if (beta==0.0)
                outStr.append("N/A");
            else
                outStr.append(doubleToStr(probXGivenInv[i],PARAM_DECI));
		}
		for (int i=0; i<4*num_alpha; i++) {
			if (i%4==0)
				outStr.append("\nSite_category_" + intToStr((int)i/4+1));
			outStr.append("\t" + doubleToStr(rootNodeFreq[i],PARAM_DECI));
		}
		outStr.append("\n\n");
		outStr.append("[Proportion of different sites]\n");
		outStr.append("Invariable_site\t" + doubleToStr(beta,PARAM_DECI) + "\n");
		for (int i=0; i<num_alpha; i++) {
			outStr.append("Site_category_" + intToStr(i+1) + "\t" + doubleToStr(alpha[i],PARAM_DECI) + "\n");
		}
}

void VariableSet::showContent(ofstream& fout) {
    // show the content
    string str = "";
    showContent(str);
    fout << str << endl;
}

void ParameterSet::readParamFile(char* fileName, vector<string>& nodeList, vector<double>& edgeLens, 
                                                        int catID, int num_edges, int edgeRepresent) {
        
    // read the values of w and pi from the parameter file
    // For reversible, there are six w values (i.e. num_w should be 6)

    // edgeRepresent: representation of the edge length
    //                1 - Average number of substitutions per site
    //                2 - Time
    
    ifstream fin;
    fin.open(fileName);
    if (!fin.is_open()) {
        cerr << "Error opening the file :" << fileName << endl;
        exit(1);
    }
    
    // this function is designed for num_w = 6
    int num_w = 6;
    
    // allocate the space to the arrays
    this->num_w = num_w;
    this->num_edges = num_edges;
    if (w == NULL)
        w = matrix(num_edges, num_w);
    if (pi == NULL)
        pi = matrix(num_edges, 4);
    if (t == NULL)
        t = matrix(num_edges, 1);
    
    // read the file
    vector<string> token;
    string aline;
    string nodeName;
    int nodeID;
    double base;
    double* curr_w;
    double* curr_pi;
    int* accessed = new int[num_edges];
    memset(accessed, 0, num_edges*sizeof(int));
    while(getline(fin,aline)) {
        if (aline.length() > 0 && aline[0]!='#') {
            tokenizer(aline, " \t, ", &token);
            if (token.size() <= num_w+4) {
                cerr << "Error! The number of items are too few in the following line inside the file " << fileName << " :" << endl;
                cerr << aline << endl;
                exit(1);
            }
            // get the node name and the node ID
            nodeName = token[0];
            for (nodeID=0; nodeID<num_edges; nodeID++) {
                if (nodeName == nodeList[nodeID])
                    break;
            }
            if (nodeID >= num_edges) {
                cerr << "Error! The node name : " << nodeName << " does not appear in the tree file" << endl;
                exit(1);
            }
            if (accessed[nodeID]) {
                cerr << "Error! The node name : " << nodeName << " appears twice in the file " << fileName << endl;
                exit(1);
            }
            accessed[nodeID] = 1;
            curr_w = &(w[nodeID*num_w]);
            curr_pi = &(pi[nodeID*4]);
            for (int i=0; i<num_w; i++) {
                curr_w[i] = atof(token[i+1].c_str());
            }
            for (int i=0; i<4; i++) {
                curr_pi[i] = atof(token[num_w+i+1].c_str());
            }
            if (edgeRepresent == 1) {
                // edge length represents the average number of substitutions per site
                base = (curr_pi[0]*curr_pi[1]*curr_w[0] + curr_pi[0]*curr_pi[2]*curr_w[1] + curr_pi[0]*curr_pi[3]*curr_w[2] + curr_pi[1]*curr_pi[2]*curr_w[3] + curr_pi[1]*curr_pi[3]*curr_w[4] + curr_pi[2]*curr_pi[3]*curr_w[5])*2.0;
                t[nodeID] = edgeLens[catID*num_edges + nodeID] / base;
            } else {
                // edge length represents the time
                t[nodeID] = edgeLens[catID*num_edges + nodeID];
            }
            // normalize
            /*
            for (int i=0; i<num_w; i++) {
                curr_w[i] = curr_w[i] / base;
            }
            t[nodeID] = edgeLens[catID*num_edges + nodeID];
             */
        }
    }
    fin.close();
    // check whether there is missing node name
    for (nodeID=0; nodeID<num_edges; nodeID++) {
        if (accessed[nodeID] == 0) {
            cerr << "Error! The node name : " << nodeList[nodeID] << " does not appear in the file " << fileName << endl;
            exit(1);
        }
    }
    delete[] accessed;
}

void ParameterSet::initialize(int num_w, int num_edges) {
    // initialize the parameters
    
    // allocate the space to the arrays
    this->num_w = num_w;
    this->num_edges = num_edges;
    
    if (w == NULL) {
        w = matrix(num_edges, num_w);
        // set all values of w to 1
        setMatrixToVal(w, num_w, num_edges, 1.0);
    }
    if (pi == NULL) {
        pi = matrix(num_edges, 4);
        // set all values of pi to 0.25
        setMatrixToVal(pi, 4, num_edges, 0.25);
    }
    if (t == NULL) {
        t = matrix(num_edges, 1);
        // set all values of t to 1
        setMatrixToVal(t, 1, num_edges, 1.0);
    }
    // initalize the rate matrix ID
    initializeRateMatID();
}

ParameterSet::ParameterSet() {
    // constructor
    w = NULL;
    pi = NULL;
    t = NULL;
    allEigenMat=NULL;
    allCondProb=NULL;
    
    // for temporary usage
    s = NULL;
    eigenVals = NULL;
    eigenVects = NULL;
    tmpMat1 = NULL;
    tmpMat2 = NULL;
    sqrtPi = NULL;
}

ParameterSet::~ParameterSet() {
    // destructor
    if (w!=NULL)
        free(w);
    if (pi!=NULL)
        free(pi);
    if (t!=NULL)
        free(t);
    if (allEigenMat!=NULL)
        free(allEigenMat);
    if (allCondProb!=NULL)
        free(allCondProb);
    for (int i=0; i<rateIDs.size(); i++) {
        delete(rateIDs[i]); // checked
    }
    
    // for temporary usage
    if (s!=NULL)
        free(s);
    if (eigenVals!=NULL)
        free(eigenVals);
    if (eigenVects!=NULL)
        free(eigenVects);
    if (tmpMat1!=NULL)
        free(tmpMat1);
    if (tmpMat2!=NULL)
        free(tmpMat2);
    
}

// update the parameters such that
// outFormat = 1 : Edge length is set to the rate of substitution; OR
//             2 : S6 in the rate matrix is set to 1
void ParameterSet::updateContent(int outFormat) {
    if (w!=NULL && pi!=NULL && t!=NULL) {
        for (int i=0; i<num_edges; i++) {
            if (num_w == 6) {
                double* curr_w = &(w[i*6]);
                double* curr_pi = &(pi[i*4]);
                double base = 1.0;
				if (outFormat==1)
					base = (curr_pi[0]*curr_pi[1]*curr_w[0] + curr_pi[0]*curr_pi[2]*curr_w[1] + curr_pi[0]*curr_pi[3]*curr_w[2] + curr_pi[1]*curr_pi[2]*curr_w[3] + curr_pi[1]*curr_pi[3]*curr_w[4] + curr_pi[2]*curr_pi[3]*curr_w[5])*2.0;
				else if (outFormat==2)
					base = curr_w[num_w-1];
                for (int j=0; j<6; j++) {
					curr_w[j] = curr_w[j]/base;
                }
				t[i] = t[i]*base;
			}
		}
	}
}

void ParameterSet::showContent(string &outStr, int* topMatrix, vector<string>* leafList) {
	
	if (topMatrix!=NULL) {
		if (w!=NULL) {
			// show the information of the rate matrix
			outStr.append("[GTR rate parameters of edge leading to each node]\n");
			outStr.append("      \t");
			outStr.append("S1\tS2\tS3\tS4\tS5\tS6\n");
			for (int i=0; i<num_edges; i++) {
				if (topMatrix!=NULL) {
					if (topMatrix[i] < 0) {
						// a leaf
						outStr.append(leafList->at(-topMatrix[i]-1));
					} else {
						// an internal node
						outStr.append(intToStr(topMatrix[i]));
					}
				}
				for (int j=0; j<num_w; j++) {
					outStr.append("\t" + doubleToStr(w[i*num_w+j],PARAM_DECI));
				}
				outStr.append("\n");
			}
		}
		outStr.append("\n");
		if (pi!=NULL) {
			// show the information of Nucleotide distribution
			outStr.append("[Nucleotide distribution of edge leading to each node]\n");
			if (topMatrix!=NULL)
				outStr.append("      \t");
			outStr.append("Pi_1\tPi_2\tPi_3\tPi_4\n");
			for (int i=0; i<num_edges; i++) {
				if (topMatrix!=NULL) {
					if (topMatrix[i] < 0) {
						// a leaf
						outStr.append(leafList->at(-topMatrix[i]-1));
					} else {
						// an internal node
						outStr.append(intToStr(topMatrix[i]));
					}
				}
				for (int j=0; j<4; j++) {
					outStr.append("\t" + doubleToStr(pi[i*4+j],PARAM_DECI));
				}
				outStr.append("\n");
			}
		}
		outStr.append("\n");
	} else {
		if (w!=NULL && pi!=NULL && t!=NULL) {
			for (int i=0; i<num_edges; i++) {
				for (int j=0; j<num_w; j++) {
					outStr.append(doubleToStr(w[i*num_w+j],PARAM_DECI) + "\t");
				}
				for (int j=0; j<4; j++) {
					outStr.append(doubleToStr(pi[i*4+j],PARAM_DECI) + "\t");
				}
				outStr.append(doubleToStr(t[i],PARAM_DECI) + "\n");
			}
		}
	}
}


void ParameterSet::showContent() {
    if (w!=NULL && pi!=NULL && t!=NULL) {
        for (int i=0; i<num_edges; i++) {
            for (int j=0; j<num_w; j++) {
                cout << w[i*num_w+j] << "\t";
            }
            for (int j=0; j<4; j++) {
                cout << pi[i*4+j] << "\t";
            }
            cout << t[i] << endl;
        }
    }
}

void ParameterSet::showContent(ofstream& fout) {
    if (w!=NULL && pi!=NULL && t!=NULL) {
        for (int i=0; i<num_edges; i++) {
            for (int j=0; j<num_w; j++) {
                fout << w[i*num_w+j] << "\t";
            }
            for (int j=0; j<4; j++) {
                fout << pi[i*4+j] << "\t";
            }
            fout << t[i] << endl;
        }
    }
}

struct paramComp {
    double* w;
    double* pi;
    int num_w;
    bool operator() (const paramComp& edge1, const paramComp& edge2) const {
        int i=0;
        while (i<num_w && edge1.w[i]==edge2.w[i])
            i++;
        if (i<num_w)
            return (edge1.w[i]<edge2.w[i]);
        int k=0;
        while (k<4 && edge1.pi[k]==edge2.pi[k])
            k++;
        if (k<4)
            return (edge1.pi[i]<edge2.pi[k]);
        else
            return false;
    }
};


// get the initialize ID of rate matrices for each edge
// before: the values of the parameters (i.e. w and pi) have been read from the file
void ParameterSet::initializeRateMatID() {
    map<paramComp,int,paramComp> paramList;
    map<paramComp,int,paramComp>::iterator itr;
    for (int i=0; i<num_edges; i++) {
        paramComp currParam;
        currParam.num_w = num_w;
        currParam.w = &(w[num_w * i]);
        currParam.pi = &(pi[4 * i]);
        itr = paramList.find(currParam);
        if (itr!=paramList.end()) {
            // found
            rateIDs[itr->second]->push_back(i);
        } else {
            // new
            int newID = (int) paramList.size();
            rateIDs.push_back(new vector<int>);
            rateIDs[newID]->push_back(i);
            paramList.insert(pair<paramComp,int>(currParam,newID));
        }
    }
    num_rate_matrices = (int) paramList.size();
}

void ParameterSet::printRateMatID() {
    for (int i=0; i<rateIDs.size(); i++) {
        cout << "Rate matrix " << i << " :";
        for (int j=0; j<rateIDs[i]->size(); j++) {
            cout << " " << rateIDs[i]->at(j);
        }
        cout << endl;
    }
}

// output the rate matrice in single array form
void ParameterSet::outputRateMat(vector<int>& outRateMatArray) {
    outRateMatArray.assign(num_edges,0);
    for (int i=0; i<rateIDs.size(); i++) {
        for (int j=0; j<rateIDs[i]->size(); j++) {
            outRateMatArray[rateIDs[i]->at(j)]=i+1;
        }
    }
}

// update the rate matrix
void ParameterSet::updateRateMat(vector<int>& rateMatArray, int numRateGrp) {
    // clear the rateIDs
    for (int i=0; i<rateIDs.size(); i++) {
        rateIDs[i]->clear();
        delete(rateIDs[i]); // checked
    }
    rateIDs.clear();
    // initialize the rateIDs
    for (int i=0; i<numRateGrp; i++) {
        vector<int> * new_vector = new vector<int>;
        rateIDs.push_back(new_vector);
    }
    // set the rateIDs
    for (int i=0; i<rateMatArray.size(); i++) {
        rateIDs[rateMatArray[i]-1]->push_back(i);
    }
    num_rate_matrices=numRateGrp;
}

// load the rate matrix from the file
// prerequisite: the value of "num_edges" has to be set
void ParameterSet::loadRateMat(char* rateGrpFile) {
    ifstream fin;
    fin.open(rateGrpFile);
    if (!fin.is_open()) {
        cerr << "Error opening the rate matrix file :" << rateGrpFile << endl;
        exit(1);
    }
    string rateMatrixStr;
    // get the rate matrix from the first line of the file
    if (!getline(fin, rateMatrixStr)) {
        cerr << "Error! The rate matrix file " << rateGrpFile << " is empty" << endl;
        exit(1);
    }
    fin.close(); // close the file
    
    trim(rateMatrixStr);
    vector<string> token;
    tokenizer(rateMatrixStr, ",\t ", &token);
    if (token.size() < num_edges) {
        cerr << "Error! The number of items in the first line of the rate matrix file " << rateGrpFile << " is less than the number of edges : " << num_edges << endl;
        exit(1);
    }
    
    int maxRateGrp = 0;
    vector<int> rateMatArray;
    for (int i=0; i<num_edges; i++) {
        int rateGrp = atoi(token[i].c_str());
        rateMatArray.push_back(rateGrp);
        if (maxRateGrp < rateGrp)
            maxRateGrp = rateGrp;
    }
    
    updateRateMat(rateMatArray, maxRateGrp);
}


// compute the eigen-matrix for all edges
// output: allEigenMat
void ParameterSet::computeAllEigenMatrix() {
    
    // allocate memory to eigen-matrix
    if (allEigenMat==NULL) {
        allEigenMat = matrix(num_edges,16);
    }
    
    if (s==NULL)
        s = matrix(4,4);
    
    for (int edge=0; edge<num_edges; edge++) {
        double* curr_w = &(w[edge*num_w]);
        double* curr_pi = &(pi[edge*4]);
        double* curr_eigenMat = &(allEigenMat[edge*16]);
        computeSMatrix(s, curr_w, curr_pi);
        computeEigenMatrix(curr_eigenMat, curr_w, curr_pi, s);
    }
    
}

// OBJECTIVE:             compute the conditional probabilities for all the edges.
// OUTPUT:                allCondProb
// SIMILAR FUNCTION IN R: getCondProb
// PREREQUISITE:          computeAllEigenMatrix()
void ParameterSet::computeAllCondProb(int isReversible) {
    if (!isReversible) {
        cerr << "Error! the function ""computeCondProb"" cannot support not isreversible condition yet." << endl;
    }
    if (allCondProb==NULL) {
        allCondProb=matrix(num_edges, 16);
    }
    
    if (eigenVals==NULL)
        eigenVals = matrix(1,4);
    if (eigenVects==NULL)
        eigenVects = matrix(4,4);
    if (tmpMat1==NULL)
        tmpMat1 = matrix(4,4);
    if (tmpMat2==NULL)
        tmpMat2 = matrix(4,4);
    if (sqrtPi==NULL)
        sqrtPi = matrix(1,4);
    
    int numRound;
    
    for (int edge=0; edge<num_edges; edge++) {
        double* condProb = &(allCondProb[edge*16]);
        double* eigenMat = &(allEigenMat[edge*16]);
        double* currPi = &(pi[edge*4]);
        double time = t[edge];
        
        // compute the eigenvalues and eigenvectors
        jacobi(eigenMat, 4, eigenVals, eigenVects, numRound);
        
        // tmpMat1 = eigen-vector * diag(exp(eigen-values * time))
        for (int j=0; j<4; j++) {
            double f = exp(eigenVals[j] * time);
            for (int i=0; i<4; i++) {
                tmpMat1[i*4+j] = f * eigenVects[i*4+j];
            }
        }
        
        // tmpMat2 = transpose of eigen-vector
        transpose(tmpMat2, eigenVects, 4, 4);
        
        // condProb = tmpMat1 * tmpMat2
        multiply(condProb, tmpMat1, tmpMat2, 4, 4, 4);
        
        // sqrtPi = sqrt(Pi)
        for (int i=0; i<4; i++)
            sqrtPi[i] = sqrt(currPi[i]);
        
        // condProb = diag(invserse(sqrtPi)) * condProb * diag(sqrtPi)
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                if (i==j)
                    continue;
                condProb[i*4+j] = condProb[i*4+j] * sqrtPi[j] / sqrtPi[i];
            }
        }
    }
}

void ParameterSet::printAllCondProb() {
    if (allCondProb != NULL) {
        // print out in R format
        for (int k=0; k<4; k++) {
            cout << endl << ", , " << k+1 << endl;
            cout << "\t[,1]\t[,2]\t[,3]\t[,4]" << endl;
            for (int i=0; i<num_edges; i++) {
                cout << "[" << i+1 << ",]";
                for (int j=0; j<4; j++) {
                    cout << "\t" << allCondProb[i*16+j*4+k];
                }
                cout << endl;
            }
        }
    }
}

void ParameterSet::printAllCondProb2() {
    if (allCondProb != NULL) {
        
        for (int i=0; i<num_edges; i++) {
            cout << "Edge " << i << endl;
            for (int j=0; j<4; j++) {
                for (int k=0; k<4; k++) {
                    cout << "\t" << allCondProb[i*16+j*4+k];
                }
                cout << endl;
            }
            cout << endl;
        }
    }
}

void ParameterSet::print_t() {
    for (int i=0; i<num_edges; i++) {
        cout << t[i] << " ";
    }
    cout << endl;
}


// copy the whole content (except rateIDs) from another instance of ParameterSet
void ParameterSet::copyFrom(ParameterSet& ps) {
    
    num_w = ps.num_w;
    num_edges = ps.num_edges;

    if (w == NULL) {
        w = matrix(num_edges, num_w);
    }
    memcpy(w,ps.w,num_w*num_edges*sizeof(double));
           
    if (pi == NULL) {
        pi = matrix(num_edges, 4);
    }
    memcpy(pi,ps.pi,num_edges*4*sizeof(double));
    
    if (t == NULL) {
        t = matrix(num_edges, 1);
    }
    memcpy(t,ps.t,num_edges*sizeof(double));
    
    if (allEigenMat==NULL) {
        allEigenMat = matrix(num_edges,16);
    }
    memcpy(allEigenMat,ps.allEigenMat,num_edges*16*sizeof(double));

    if (allCondProb==NULL) {
        allCondProb=matrix(num_edges, 16);
    }
    memcpy(allCondProb,ps.allCondProb,num_edges*16*sizeof(double));
}



// constructor
AllParameterSet::AllParameterSet(int numRateCat) {
    this->numRateCat = numRateCat;
    // then create the ParameterSet objects
    for (int i=0; i<numRateCat; i++) {
        ParameterSet* a_ps = new ParameterSet();
        ps.push_back(a_ps);
    }
}

// destructor
AllParameterSet::~AllParameterSet() {
    for (int i=0; i<numRateCat; i++) {
        delete (ps[i]); // checked
    }
    ps.clear();
}

// read the values of w and pi from the parameter file list
void AllParameterSet::readParamFileList(char* fileName, vector<string>& categoryList, vector<string>& nodeList, 
                                        vector<double>& edgeLens, int num_edges, int edgeRepresent) {
    int i = 0;
    string aline;
    ifstream fin;
    fin.open(fileName);
    if (!fin.is_open()) {
        cerr << "Error! Cannot open the file : " << fileName << endl;
        exit(1);
    }
    vector<string> token;
    int catID;
    int catNum = categoryList.size();
    int* accessed = new int[catNum];
    memset(accessed, 0, catNum*sizeof(int));
    while (getline(fin, aline)) {
        if (aline.length()==0 || aline[0]=='#')
            continue;
        tokenizer(aline, "\t ", &token);
        if (token.size() >= 2) {
            for (catID=0; catID<catNum; catID++) {
                if (token[0] == categoryList[catID])
                    break;
            }
            if (catID>=catNum) {
                cerr << "Error! The category name : " << token[0] << " does not appear in the tree file" << endl;
                exit(1);
            }
            if (accessed[catID]) {
                cerr << "Error! The category name : " << token[0] << " appears twice in the file: " << fileName << endl;
                exit(1);
            }
            accessed[catID] = 1;
            // read the values of w and pi from the parameter file list
            ps[catID]->readParamFile((char*) token[1].c_str(), nodeList, edgeLens, catID, num_edges, edgeRepresent);
        }
    }
    // check whether there is category which does not appear in the parameter file list
    for (catID=0; catID<catNum; catID++) {
        if (accessed[catID]==0) {
            cerr << "Error! The cateogry name : " << categoryList[catID] << " does not appear in the file: " << fileName << endl;
            exit(1);
        }
    }
    fin.close();
    delete[] accessed;
}

// load the file of rate groups
// return the number of rate groups
int AllParameterSet::loadRateMat(char* rateGrpFile) {
    if (numRateCat == 0) {
        cerr << "[allParameterSet::loadRateMat] Error! The number of rate category is 0!" << endl;
    }
    for (int i=0; i<numRateCat; i++) {
        ps[i]->loadRateMat(rateGrpFile);
    }
    return ps[0]->num_rate_matrices;;
}


void AllParameterSet::showContent(string& outStr, int* topMatrix, vector<string>* leafList) {
    for (int i=0; i<numRateCat; i++) {
		outStr.append("--------------------------------\n");
        outStr.append("Site category " + intToStr(i+1) + "\n");
		outStr.append("--------------------------------\n");
        ps[i]->showContent(outStr, topMatrix,  leafList);
    }
}

void AllParameterSet::showContent() {
    for (int i=0; i<numRateCat; i++) {
        cout << "[parameter set " << i << "]" <<  endl;
        ps[i]->showContent();
    }
}

void AllParameterSet::updateContent(int outFormat) {
    for (int i=0; i<numRateCat; i++) {
        ps[i]->updateContent(outFormat);
    }
}

void AllParameterSet::computeAllEigenMatrix() {
    for (int i=0; i<numRateCat; i++) {
        ps[i]->computeAllEigenMatrix();
    }
}

void AllParameterSet::computeAllCondProb(int isReversible) {
    allCondProbSet.clear();
    for (int i=0; i<numRateCat; i++) {
        ps[i]->computeAllCondProb(isReversible);
        allCondProbSet.push_back(ps[i]->allCondProb);
    }
}

void AllParameterSet::printAllCondProb() {
    for (int k=0; k<4; k++) {
        for (int j=0; j<4; j++) {
            cout << endl;
            cout << ", ," << j+1 << "," << k+1 << endl;
            cout << endl;
            for (int s=0; s<2; s++) {
                for (int i=s*7; i<(s+1)*7; i++) {
                    cout << "\t[," << i+1 << "]";
                }
                cout << endl;
                for (int p=0; p<numRateCat; p++) {
                    cout << "[" << p+1 << ",]";
                    for (int i=s*7; i<(s+1)*7; i++) {
                        cout << "\t" << ps[p]->allCondProb[i*16+j*4+k];
                    }
                    cout << endl;
                }
            }
        }
    }
}

void AllParameterSet::printAllCondProb2() {
    
    for (int i=0; i<ps.size(); i++) {
        ps[i]->printAllCondProb2();
        cout << endl;
    }
}

int AllParameterSet::size() {
    return (int) ps.size();
}

// compute the s-matrix
// output: s_matrix
void computeSMatrix(double* s_matrix, double* w, double * pi) {
    // compute the s-matrix (dimension: 4 x 4)
    for (int i=0; i<3; i++)
        s_matrix[0*4+(i+1)]=s_matrix[(i+1)*4+0]=w[i]; //w1,w2,w3
    s_matrix[1*4+2]=s_matrix[2*4+1]=w[3]; //w4
    s_matrix[1*4+3]=s_matrix[3*4+1]=w[4]; //w5
    s_matrix[2*4+3]=s_matrix[3*4+2]=w[5]; //w6
    for (int i=0; i<4; i++) {
        s_matrix[i*4+i]=0.0;
        s_matrix[i*4+i]=(- sumProduct(&(s_matrix[i*4]), pi, 4)) / pi[i];
    }
}

// compute the eigen-matrix
// output: eigenMat
void computeEigenMatrix(double* eigenMat, double* w, double* pi, double* s) {
    int k = 0;
    for (int i=0; i<3; i++) {
        for (int j=i+1; j<4; j++) {
            // (i,j) = sqrt(pi[i] * pi[j]) * w[k]
            eigenMat[i*4+j] = sqrt(pi[i] * pi[j]) * w[k++];
        }
    }
    for (int i=0; i<4; i++)
        // (i,i) = pi[i] * s[i][j]
        eigenMat[i*4+i] = pi[i] * s[i*4+i];
    for (int i=1; i<4; i++)
        for (int j=0; j<i; j++)
            // (i,j) = (j,i)
            eigenMat[i*4+j] = eigenMat[j*4+i];
}

