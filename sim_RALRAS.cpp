/*
 *
 * sim_RALRAS.cpp
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

#include "sim_RALRAS.h"

// helping function for displaying the topology matrix into a tree format
void subTopMatrixToTreeFormat(int* topMatrix, int line, string& str, vector<string>* leafList, vector<double>& edgeLens, int catID, int numEdges) {
    int node1 = topMatrix[line*2];
    int node2 = topMatrix[line*2+1];
	str.append("(");
    if (node1 < 0) {
        // a leaf
        if (leafList!=NULL)
            str.append(leafList->at(-node1-1));
    } else {
        // an internal node
        subTopMatrixToTreeFormat(topMatrix,node1-1, str, leafList, edgeLens, catID, numEdges);
    }
    str.append(":" + doubleToStr(edgeLens[catID*numEdges+line*2],PARAM_DECI));
    str.append(",");
    if (node2 < 0) {
        // a leaf
        if (leafList!=NULL)
            str.append(leafList->at(-node2-1));
    } else {
        // an internal node
        subTopMatrixToTreeFormat(topMatrix,node2-1, str, leafList, edgeLens, catID, numEdges);
    }
    str.append(":" + doubleToStr(edgeLens[catID*numEdges+line*2+1],PARAM_DECI));
    str.append(")");
}

// display the arrangement of topology matrix into a tree format
string topMatrixToTreeFormat(int* topMatrix, vector<string>* leafList, vector<double>& edgeLens, int numCat) {
    int numLineTopMatrix=leafList->size()-1;
    int numEdges = 2 * numLineTopMatrix;
    string str = "";
	
    for (int catID=0; catID<numCat; catID++) {
        str.append("Tree " + intToStr(catID+1) + ":\n");
        subTopMatrixToTreeFormat(topMatrix, numLineTopMatrix-1, str, leafList, edgeLens, catID, numEdges);
        str.append(";\n");
    }
	return str;	
}


int howManyCharInside(string str, char c) {
    int n=0;
    for (int i=0; i<str.length(); i++) {
        if (str[i]==c)
            n++;
    }
    return n;
}

// Input: topology string (can also include the rate group and the edge lengths)
// Output: 1. number of rows of the matrix, 2. List of leaf names,
//         3. edge length, 4. List of node names
// Return: Topology matrix
int* genTopMatrixFrStr(string newickStr, int& rowNum, vector<string>& leafList,
                       vector<double>& edgeLens, vector<string>& nodeList) {
    
    trim(newickStr);
    if (newickStr.length()==0) {
        cerr << "Error! The topology string is empty" << endl;
        exit(1);
    }
    
    // get the number of rows in the matrix
    rowNum = howManyCharInside(newickStr, ',');
    
    if (rowNum == 0) {
        return NULL;
    }
    
    // clear all the arrays
    leafList.clear();
    edgeLens.clear();
    nodeList.clear();
    
    // allocate memory for the topMatrix
    int* topMatrix = new int[rowNum*2];
    
    // create another array to store the nodes waiting to be processed
    vector<int> waitNodeIDList;
    vector<string> waitNodeNameList;
    vector<double> waitNodeDist;
    
    // the number of new nodes formed
    int newNodeNum = 0;
    
    // current node's name and its edge length
    string currNodeName = "";
    string currDist = "";
    int distMode = 0;
    int isLeaf = 1;
    int i;
    for (i=0; i<(int)newickStr.length(); i++) {
        char c = newickStr.at(i);
        if (c == ';')
            break;
        if (c == '(' || c == ' ') {
            continue; // do nothing
        } else if (c == ',' || c==')') {
            if (currNodeName=="") {
                // no node name is provided
                cerr << "There is missing node name in the tree file" << endl;
                exit(1);
            }
            if (currDist=="") {
                // no edge length information is provided
                cerr << "No edge length information is provided for the node " << currNodeName << endl;
                exit(1);
            }
            // add the current node name to waitNodeNameList
            waitNodeNameList.push_back(currNodeName);
            // add the current edge length to waitNodeDist
            waitNodeDist.push_back((double) atof(currDist.c_str()));
            if (isLeaf) {
                // add the current leaf to leafList
                leafList.push_back(currNodeName);
                // add the current node ID to waitNodeIDList
                waitNodeIDList.push_back(-1 * (int)(leafList.size()));
            } else {
                // add the current node ID to waitNodeIDList
                waitNodeIDList.push_back(++newNodeNum);
            }
            // the following is applied for ')' only
            if (c == ')') {
                // add a new row into the topMatrix
                int lastID = waitNodeIDList.back();
                waitNodeIDList.pop_back();
                topMatrix[newNodeNum*2] = waitNodeIDList.back();
                waitNodeIDList.pop_back();
                topMatrix[newNodeNum*2+1] = lastID;
                // add the corresponding edge lengths into the edgeLens
                double lastDist = waitNodeDist.back();
                waitNodeDist.pop_back();
                edgeLens.push_back(waitNodeDist.back());
                waitNodeDist.pop_back();
                edgeLens.push_back(lastDist);
                // add the corresponding node name into the nodeList
                string lastName = waitNodeNameList.back();
                waitNodeNameList.pop_back();
                nodeList.push_back(waitNodeNameList.back());
                waitNodeNameList.pop_back();
                nodeList.push_back(lastName);
            }
            // reset the variables
            currNodeName = "";
            currDist = "";
            distMode = 0;
            if (c==')')
                isLeaf = 0;
            else
                isLeaf = 1;
        } else if (c == ':') {
            distMode = 1;
        } else if (distMode==1) {
            currDist.append(1,c);
        } else {
            currNodeName.append(1,c);
        }
    }
    if (waitNodeNameList.size() != 0) {
    	cerr << "Error! The input tree is unrooted." << endl;
    	cerr << "This program only supports a rooted tree" << endl;
    	exit(1);
    }
    return topMatrix;
}

// Input: Trees file
// Output: 1. number of rows of the matrix, 2. List of leaf names,
//         3. edge length (# of items = # of edges * # of site categories)
//         4. List of node names, 5. List of category
// return: Topology matrix
int* genTopMatrix(char* treesFile, int& rowNum, vector<string>& leafList, vector<double>& edgeLens,
                            vector<string>& nodeList, vector<string>& siteCatNames) {
    
    // read the trees file
    ifstream fin;
    fin.open(treesFile); // open the file
    if (!fin.is_open()) {
        cerr << "Error opening the topology file :" << treesFile << endl;
        exit(1);
    }
    string aline;
    int* firstTopMatrix = NULL;
    int* currTopMatrix = NULL;
    int currRowNum;
    vector<string> currLeafList;
    vector<double> currEdgeLens;
    vector<string> currNodeList;
    bool sameTree = true;
    int i;
    vector<string> token;
    
    while (getline(fin, aline) && sameTree) {
        if (aline.length()>0 && aline[0]=='#')
            continue;
        tokenizer(aline, " \t", &token);
        if (token.size() >= 2) {
            if (siteCatNames.size()==0) {
                firstTopMatrix = genTopMatrixFrStr(token[1], rowNum, leafList, edgeLens, nodeList);                
                if (rowNum>0) {
                    siteCatNames.push_back(token[0]);
                }
                /*
                // print out the topMatrix
                cout << "topMatrix:" << endl << flush;
                for (i=0; i<rowNum*2; i+=2) {
                    cout << firstTopMatrix[i] << "," << firstTopMatrix[i+1] << endl;
                }
                // print out the leafList
                cout << "leafList:" << endl;
                for (i=0; i<rowNum+1; i++) {
                    if (i>0)
                        cout << ",";
                    cout << leafList[i];
                }
                cout << endl;
                // print out the edge lengths
                cout << "edgeLens:" << endl;
                for (i=0; i<rowNum*2; i++) {
                    if (i>0)
                        cout << ",";
                    cout << edgeLens[i];
                }
                cout << endl;
                // print out the nodeList
                cout << "nodeList:" << endl;
                for (i=0; i<rowNum*2; i++) {
                    if (i>0)
                        cout << ",";
                    cout << nodeList[i];
                }
                cout << endl;*/
            } else {
                currTopMatrix = genTopMatrixFrStr(token[1], currRowNum, currLeafList, currEdgeLens, currNodeList);
                /*
                // print out the topMatrix
                cout << "currTopMatrix:" << endl;
                for (i=0; i<currRowNum*2; i+=2) {
                    cout << currTopMatrix[i] << "," << currTopMatrix[i+1] << endl;
                }
                // print out the leafList
                cout << "currLeafList:" << endl;
                for (i=0; i<currRowNum+1; i++) {
                    if (i>0)
                        cout << ",";
                    cout << currLeafList[i];
                }
                cout << endl;
                // print out the edge lengths
                cout << "currEdgeLens:" << endl;
                for (i=0; i<currRowNum*2; i++) {
                    if (i>0)
                        cout << ",";
                    cout << currEdgeLens[i];
                }
                cout << endl;
                // print out the nodeList
                cout << "currNodeList:" << endl;
                for (i=0; i<currRowNum*2; i++) {
                    if (i>0)
                        cout << ",";
                    cout << currNodeList[i];
                }
                cout << endl;*/
                if (currRowNum>0) {
                    // check whether all the tree topologies are the same
                    if (currRowNum != rowNum) {
                        sameTree = false;
                        break;
                    }
                    // topMatrix vs currTopMatrix
                    for (i=0; i<currRowNum*2; i++) {
                        if (firstTopMatrix[i] != currTopMatrix[i]) {
                            sameTree = false;
                            break;
                        }
                    }
                    // leafList vs currLeafList
                    if (sameTree) {
                        for (i=0; i<currRowNum+1; i++) {
                            if (leafList[i] != currLeafList[i]) {
                                sameTree = false;
                                break;
                            }
                        }
                    }
                    // nodeList vs currNodeList
                    if (sameTree) {
                        for (i=0; i<currRowNum*2; i++) {
                            if (nodeList[i] != currNodeList[i]) {
                                sameTree = false;
                                break;
                            }
                        }
                    }
                    // if the topology is not the same, then give out error message
                    if (!sameTree) {
                        cerr << "Error! The topoloies are not the same inside the tree file" << endl;
                        exit(1);
                    }
                    // append the new items in currEdgeLens to edgeLens
                    for (i=0; i<currRowNum*2; i++) {
                        edgeLens.push_back(currEdgeLens[i]);
                    }
                    siteCatNames.push_back(token[0]);
                    // reset the variables
                    delete[] currTopMatrix;
                    currTopMatrix = NULL;
                }
            }
        }
    }
    fin.close(); // close the file
    return firstTopMatrix;
}

void genRandSeq(char* seq, int start, int num, double* ratios, int numChars) {
    for (int i=start; i<start+num; i++) {
        int base=1000;
        int rand_num = rand()%base+1;
        int j;
        int acc = 0;
        // cout << "acc : ";
        for (j=0; j<numChars; j++) {
            acc += (int) round(ratios[j]*base);
            // cout << acc << ",";
            if (j==numChars-1) {
                seq[i]=j;
            } else if (rand_num <= acc) {
                seq[i]=j;
                break;
            }
        }
        // cout << endl;
    }
    // show the statistics
    int* stat = new int[numChars];
    for (int i=0; i<numChars; i++)
        stat[i] = 0;
    for (int i=start; i<start+num; i++) {
        stat[seq[i]]++;
    }
    cout << "[" << start+1 << "-" << start+num << "] ";
    for (int i=0; i<numChars; i++) {
        cout << i << ":" << stat[i] << "; ";
    }
    cout << endl;
    
}


int* genRandIntArray(int size) {
     // generate an integer array from 0 to size-1 in random order
     int* randIntArray = new int[size];
     int i;
     for (i=0; i<size; i++)
         randIntArray[i] = i;
     for (i=0; i<size; i++) {
         int randPos = (rand()%(size-i)) + i; // a random number from i to seqLen-1
         if (randPos > i) {
             int tmp = randIntArray[i];
             randIntArray[i] = randIntArray[randPos];
             randIntArray[randPos] = tmp;
         }
     }
     return randIntArray;
}

void distNucl(char* seq, int start, int num, int* distribute) { 
    // get the nucleotide distribution 
    int k;
    for (k=0; k<4; k++)
        distribute[k] = 0;
    for (k=start; k<start+num; k++) {
        distribute[seq[k]]++;
    }
}

void genSeq(char* seq, int start, int num, double* ratios, int numChars) {
    
    int curr_start = start;
    double cum_ratio = 0.0;
    int cum_num = 0;
    int i,j;
    for (j=0; j<numChars; j++) {
        cum_ratio += ratios[j];
        cum_num = (int) round(cum_ratio*num);
        for (i=curr_start; i<cum_num+start; i++) {
            seq[i] = j;
        }
        curr_start = cum_num+start;
    }

    // show the statistics
#ifdef DEBUG_MODE
    int* stat = new int[numChars];
    distNucl(seq, start, num, stat);
    cerr << "[" << start+1 << "-" << start+num << "] (";
    for (int i=0; i<numChars; i++) {
        if (i>0)
            cerr << ", ";
        fprintf(stderr, "%7.5f", (double) stat[i] / num);
    }
    cerr << ")" << endl;
#endif
}


void genSeqFrParent(char* parentSeq, char* seq, int start, int num, double* probMatrix, int numChars) {
    
     // generate an integer array from 0 to num-1 in random order
    int* randIntArray = genRandIntArray(num);
    
    // get the distribution of the parentSeq[start, start+num-1]
    int* stat = new int[numChars];
    distNucl(parentSeq, start, num, stat);
    
    int* distribute = new int[numChars*numChars];
    int cumDist;
    int i,j,pos;
    double cumProb, currNum;
    
    for (i=0; i<numChars; i++) {
        cumProb = probMatrix[i*numChars];
        currNum = cumProb*stat[i];
        distribute[i*numChars] = (int) round(currNum);
        cumDist = distribute[i*numChars];
        
        for (j=1; j<numChars; j++) {
            cumProb += probMatrix[i*numChars+j];
            currNum = cumProb*stat[i]-cumDist;
            if (currNum < 0.0) currNum=0.0;
            distribute[i*numChars+j] = (int) round(currNum);
            cumDist += distribute[i*numChars+j];
        }
    }
    
    // generate the sequence
    char c;
    for (i=0; i<num; i++) {
        pos = start+randIntArray[i];
        c = parentSeq[pos];
        for (j=0; j<numChars; j++) {
            if (distribute[c*numChars+j] > 0) {
                seq[pos] = j;
                distribute[c*numChars+j]--;
                break;
            }
        }
    }
    
    // release the memory
    delete[] randIntArray;
    delete[] stat;
    delete[] distribute;
}


// simulate the sequences according to the parameters resulting from RAL-RAS model
void simSeqFrRASParam(char* treeFile, char* siteInfoFile, char* paramFilelist, int seqLen, int outFormat, string outPrefix, int edgeRepresent) {
    
    char int2char[] = {'A','C','G','T'};
    
    // load the topology
    int* topMatrix;
    int numLineTopMat, numSpecies, numEdges, numInterNodes, numCategories;
    vector<string> leafList;
    vector<double> edgeLens;
    vector<string> nodeList;
    vector<string> siteCatList;
    topMatrix = genTopMatrix(treeFile, numLineTopMat, leafList, edgeLens, nodeList, siteCatList);
    
    numSpecies=(int)leafList.size();
    numEdges=numSpecies*2-2;
    numInterNodes=numSpecies-1; // # of internal nodes (including root)
    numCategories=siteCatList.size();
    
    // load the variables
    VariableSet vs(numCategories);
    vs.readSiteInfoFile(siteInfoFile, &siteCatList);
    
    // load the parameter list
    AllParameterSet allPS(numCategories);
    allPS.readParamFileList(paramFilelist, siteCatList, nodeList, edgeLens, numEdges, edgeRepresent);

    // compute all the eigen matrices
    allPS.computeAllEigenMatrix();
    // compute all the conditional probabilities
    int isReversible = 1;
    allPS.computeAllCondProb(isReversible);

    
#ifdef DEBUG_MODE
    // print out the topology matrix, parameters and variables
    cout << "topology matrix:" << endl;
    for (int i=0; i<numLineTopMat; i++) {
        for (int j=0; j<2; j++) {
            if (j>0)
                cout << ",";
            cout << topMatrix[i*2+j];
        }
        cout << endl;
    }
    allPS.showContent();
    allPS.printAllCondProb2();
    vs.showContent();
    int n;
#endif
    
    
    // initialize the variables
    int numChars = 4;
    char* interSeqs = new char[numInterNodes * seqLen];
    char* leafSeqs = new char[numSpecies * seqLen];
    char* currSeq;
    char* parentSeq;
    int i,j,k;
    
    // initialize the randome nucleotide generators
    RandNuclGenerator* randNuclGen = new RandNuclGenerator[vs.num_alpha * numEdges];
    for (i=0; i<vs.num_alpha; i++) {
        for (j=0; j<numEdges; j++) {
            double* curr_prob = &(((allPS.ps)[i])->allCondProb[j*numChars*numChars]);
            randNuclGen[i*numEdges + j].init(curr_prob, numChars);
#ifdef DEBUG_MODE
            cout << "num_alpha : " << i << " EdgeID: " << j << endl;
            randNuclGen[i*numEdges + j].RandNuclGenerator::showAcc();
#endif
        }
    }

    // site categories
    //  0 : constant site
    // >0 : variable site
    int* siteCat = new int[seqLen];
    int* siteNum = new int[vs.num_alpha+1];
    int start = 0;
    int num = vs.beta*seqLen;
    siteNum[0] = num;
    for (i=start; i<start+num && i<seqLen; i++) {
        siteCat[i] = 0;
    }
#ifdef DEBUG_MODE
    // show the site categories
    cout << "Sequence length: " << seqLen << endl;
    // cout << "[" << start+1 << "-" << start+num <<"] constant site" << endl;
#endif
    for (i=0; i<vs.num_alpha; i++) {
        start += num;
        if (i==vs.num_alpha-1)
            num = seqLen - start; // last variable group
        else
            num = vs.alpha[i]*seqLen;
        siteNum[i+1] = num;
        for (j=start; j<start+num && j<seqLen; j++) {
            siteCat[j] = i+1;
        }
// #ifdef DEBUG_MODE
//        cout << "[" << start+1 << "-" << start+num <<"] variable site " << i+1 << endl;
// #endif
    }


    // generate the corresponding sites of the root sequence
    currSeq = &(interSeqs[(numInterNodes-1)*seqLen]);
    /*
    // for testing only, all are A's
    for (i=0; i<seqLen; i++)
        currSeq[i] = 0;
     */
    start = 0;
    for (i=0; i<vs.num_alpha+1; i++) {
        num = siteNum[i];
        if (i==0) {
            // constant sites
            genSeq(currSeq, start, num, vs.probXGivenInv, numChars);
        } else {
            genSeq(currSeq, start, num, &(vs.rootNodeFreq[(i-1)*numChars]), numChars);
        }
        start += num;
    }
    
    // start simulations of the sequences for each internal and terminal node
    for (i=numLineTopMat-1; i>=0; i--) {
        parentSeq = &(interSeqs[i*seqLen]);
        /*
        // for testing only, all are A's
        parentSeq = &(interSeqs[(numInterNodes-1)*seqLen]);
         */
        
        for (j=0; j<2; j++) {
            // ========================================
            // consider the left node or the right node
            // ========================================
            int currEdgeID = i*2+j;
            int currNode = topMatrix[i*2+j];
            if (currNode < 0) {
                // leaf node
                int leafID = -currNode-1;
#ifdef DEBUG_MODE
                cout << "leafID : " << leafID << endl;
#endif
                currSeq = &(leafSeqs[leafID*seqLen]);
            } else {
                // internode
                int interNodeID = currNode-1;
#ifdef DEBUG_MODE
                cout << "interNodeID : " << interNodeID << endl;
#endif
                currSeq = &(interSeqs[interNodeID*seqLen]);
            }
            
            for (k=0; k<seqLen; k++) {
                if (siteCat[k]==0) {
                    // constant site
                    currSeq[k] = parentSeq[k];
                } else {
                    // variable site
                    int  currSiteCat = siteCat[k]-1;
                    currSeq[k] = randNuclGen[currSiteCat*numEdges + currEdgeID].randomNucl(parentSeq[k]);
                }
            }
            
#ifdef DEBUG_MODE
            // show the distribution of nucl of parent sequence
            cout << "distribution of nucl of parent" << endl;
            int dist[4];
            for (k=0; k<4; k++)
                dist[k] = 0;
            for (k=0; k<seqLen; k++) {
                if (k>0 && siteCat[k-1]!=siteCat[k]) {
                    // print out the distribution
                    int tot = 0;
                    for (int p=0; p<4; p++) {
                        tot+=dist[p];
                    }
                    for (int p=0; p<4; p++) {
                        if (p>0)
                            cout << ",";
                        printf("%7.5f", (double)dist[p]/tot);
                    }
                    cout << " ";
                    // reset the numbers to zeros
                    for (int p=0; p<4; p++) {
                        dist[p]=0;
                    }
                }
                dist[parentSeq[k]]++;
            }
            // print out the distribution
            int tot = 0;
            for (int p=0; p<4; p++) {
                tot+=dist[p];
            }
            for (int p=0; p<4; p++) {
                if (p>0)
                    cout << ",";
                printf("%7.5f", (double)dist[p]/tot);
            }
            cout << endl;

            // show the distribution of nucl of current sequence
            cout << "distribution of nucl of current sequence" << endl;
            for (k=0; k<4; k++)
                dist[k] = 0;
            for (k=0; k<seqLen; k++) {
                if (k>0 && siteCat[k-1]!=siteCat[k]) {
                    // print out the distribution
                    int tot = 0;
                    for (int p=0; p<4; p++) {
                        tot+=dist[p];
                    }
                    for (int p=0; p<4; p++) {
                        if (p>0)
                            cout << ",";
                        printf("%7.5f", (double)dist[p]/tot);
                    }
                    cout << " ";
                    // reset the numbers to zeros
                    for (int p=0; p<4; p++) {
                        dist[p]=0;
                    }
                }
                dist[currSeq[k]]++;
            }
            // print out the distribution
            tot = 0;
            for (int p=0; p<4; p++) {
                tot+=dist[p];
            }
            for (int p=0; p<4; p++) {
                if (p>0)
                    cout << ",";
                printf("%7.5f", (double)dist[p]/tot);
            }
            cout << endl;
            
            
            // show the distribution of the conversion
            int conv[16];
            for (k=0; k<16; k++)
                conv[k]=0;
            for (k=0; k<seqLen; k++) {
                if (k>0 && siteCat[k-1]!=siteCat[k]) {
                    // print out the distribution
                    for (int p=0; p<4; p++) {
                        int tot = 0;
                        for (int q=0; q<4; q++) {
                            tot+=conv[p*4+q];
                        }
                        for (int q=0; q<4; q++) {
                            if (q>0)
                                cout << ",";
                            printf("%7.5f", (double)conv[p*4+q]/tot);
                        }
                        cout << endl;
                    }
                    cout << endl;
                    // reset the numbers to zeros
                    for (int p=0; p<4; p++) {
                        for (int q=0; q<4; q++) {
                            conv[p*4+q] = 0;
                        }
                    }
                }
                conv[parentSeq[k]*4+currSeq[k]]++;
            }
            // print out the distribution
            for (int p=0; p<4; p++) {
                int tot = 0;
                for (int q=0; q<4; q++) {
                    tot+=conv[p*4+q];
                }
                for (int q=0; q<4; q++) {
                    if (q>0)
                        cout << ",";
                    printf("%7.5f", (double)conv[p*4+q]/tot);
                }
                cout << endl;
            }
            cout << endl;
#endif
        }
    }
    
    /*
    // show the resulting sequences
#ifdef DEBUG_MODE
    int dist[4];
    int sum;
    cout << "nucl distribution of internal nodes" << endl;
    for (i=0; i<numInterNodes; i++) {
        cout << "[" << i << "] ";
        dist[0]=0;dist[1]=0;dist[2]=0;dist[3]=0;sum=0;
        for (j=0; j<seqLen; j++) {
            if (j>0 && siteCat[j]!=siteCat[j-1]) {
                cout << "(";
                for (n=0;n<4;n++) {
                    if (n>0)
                        cout << ",";
                    printf("%7.5f", (double) dist[n]/sum);
                }
                cout << ")";
                cout << " ";
                dist[0]=0;dist[1]=0;dist[2]=0;dist[3]=0;sum=0;
            }
            // cout << int2char[interSeqs[i*seqLen+j]];
            dist[interSeqs[i*seqLen+j]]++;sum++;
        }
        cout << "(";
        for (n=0;n<4;n++) {
            if (n>0)
                cout << ",";
            printf("%7.5f", (double) dist[n]/sum);
        }
        cout << ")";
        cout << endl;
    }
    
    cout << "nucl distribution of the leaf nodes" << endl;
    for (i=0; i<numSpecies; i++) {
        cout << "[" << i << "] ";
        dist[0]=0;dist[1]=0;dist[2]=0;dist[3]=0;sum=0;
        for (j=0; j<seqLen; j++) {
            if (j>0 && siteCat[j]!=siteCat[j-1]) {
                cout << "(";
                for (n=0;n<4;n++) {
                    if (n>0)
                        cout << ",";
                    printf("%7.5f", (double) dist[n]/sum);
                }
                cout << ")";
                cout << " ";
                dist[0]=0;dist[1]=0;dist[2]=0;dist[3]=0;sum=0;
            }
            // cout << int2char[leafSeqs[i*seqLen+j]];
            dist[leafSeqs[i*seqLen+j]]++;sum++;
        }
        cout << "(";
        for (n=0;n<4;n++) {
            if (n>0)
                cout << ",";
            printf("%7.5f", (double) dist[n]/sum);
        }
        cout << ")";
        cout << endl;
    }
#endif
     */
     
     
     // update the parameters such that
     // S6 in the rate matrix is set to 1
     allPS.updateContent(2);

     // print out the content of parameters
     string statusFileName = outPrefix + ".stat";
     ofstream statOut;
     statOut.open((char*) statusFileName.c_str());
     if (!statOut.is_open()) {
        cerr << "Error! The output status file : " << statusFileName << " cannot be created" << endl;
        exit(1);
     }
     vs.showContent(statOut);
     string str = "";
     allPS.showContent(str, topMatrix, &leafList);
     statOut << str;

     // update the parameters such that
     // Edge length is set to the rate of substitution
     allPS.updateContent(1);
     // update the edgeLen vector
     for (i=0; i<numCategories; i++) {
         for (j=0; j<numEdges; j++) {
            edgeLens[i*numEdges + j] = (allPS.ps[i])->t[j];
         }
     }

     string topStr = topMatrixToTreeFormat(topMatrix, &leafList, edgeLens, numCategories);
     statOut << topStr << endl;
     statOut.close();
        
    
     // generate an integer array from 0 to seqLen-1 in random order
    int* randIntArray = genRandIntArray(seqLen);

    /*
     // print out the integer array
     for (i=0; i<seqLen; i++)
         cerr << randIntArray[i] << endl;
     */
     // output the alignment of the leaves
     string outputFileName = outPrefix + ".out";
     ofstream fout;
     fout.open((char*) outputFileName.c_str());
     if (!fout.is_open()) {
         cerr << "Error! The output file : " << outputFileName << " cannot be created" << endl;
         exit(1);
     }
     
     int c;
     int space_len;
     if (outFormat==2) {
         fout << " " << numSpecies << " " << seqLen << endl;
         for (i=0; i<numSpecies; i++) {
             space_len = 10 - leafList[i].length();
             fout << leafList[i];
             for (j=0; j<space_len; j++)
                fout << " ";
             for (j=0; j<seqLen; j++) {
                c = randIntArray[j];
                fout << int2char[leafSeqs[i*seqLen+c]];
             }
             fout << endl;
         }
     } else {
         for (i=0; i<numSpecies; i++) {
             fout << ">" << leafList[i] << endl;
             for (j=0; j<seqLen; j++) {
                c = randIntArray[j];
                fout << int2char[leafSeqs[i*seqLen+c]];
             }
             fout << endl;
         }
     }
     fout.close();
     
     cout << endl;
     cout << "Output file: " << outputFileName << endl;
     cout << "Hetero 2 - Version " << VERSION << " finished" << endl;
    
    delete[] topMatrix;
    delete[] interSeqs;
    delete[] leafSeqs;
    delete[] randNuclGen;
    delete[] siteCat;
    delete[] siteNum;
    delete[] randIntArray;
}



// simulate PERFECT sequences according to the parameters resulting from RAL-RAS model
void simPerfectSeqFrRASParam(char* treeFile, char* siteInfoFile, char* paramFilelist, int seqLen, int outFormat, string outPrefix, int edgeRepresent) {
    
    char int2char[] = {'A','C','G','T'};
    
    // load the topology
    int* topMatrix;
    int numLineTopMat, numSpecies, numEdges, numInterNodes, numCategories;
    vector<string> leafList;
    vector<double> edgeLens;
    vector<string> nodeList;
    vector<string> siteCatList;
    topMatrix = genTopMatrix(treeFile, numLineTopMat, leafList, edgeLens, nodeList, siteCatList);
    
    numSpecies=(int)leafList.size();
    numEdges=numSpecies*2-2;
    numInterNodes=numSpecies-1; // # of internal nodes (including root)
    numCategories=siteCatList.size();
    
    // load the variables
    VariableSet vs(numCategories);
    vs.readSiteInfoFile(siteInfoFile, &siteCatList);
    
    // load the parameter list
    AllParameterSet allPS(numCategories);
    allPS.readParamFileList(paramFilelist, siteCatList, nodeList, edgeLens, numEdges, edgeRepresent);

    // compute all the eigen matrices
    allPS.computeAllEigenMatrix();
    // compute all the conditional probabilities
    int isReversible = 1;
    allPS.computeAllCondProb(isReversible);

    
#ifdef DEBUG_MODE
    // print out the topology matrix, parameters and variables
    cout << "topology matrix:" << endl;
    for (int i=0; i<numLineTopMat; i++) {
        for (int j=0; j<2; j++) {
            if (j>0)
                cout << ",";
            cout << topMatrix[i*2+j];
        }
        cout << endl;
    }
    allPS.showContent();
    allPS.printAllCondProb2();
    vs.showContent();
    int n;
#endif
    
    
    // initialize the variables
    int numChars = 4;
    char* interSeqs = new char[numInterNodes * seqLen];
    char* leafSeqs = new char[numSpecies * seqLen];
    char* currSeq;
    char* parentSeq;
    int i,j,k,l;
    
    // site categories
    //  0 : constant site
    // >0 : variable site
    int* siteNum = new int[vs.num_alpha+1];
    int* siteCat = new int[seqLen];
    int start = 0;
    int num = vs.beta*seqLen;
    siteNum[0] = num;
#ifdef DEBUG_MODE
    // show the site categories
    cout << "Sequence length: " << seqLen << endl;
#endif
    for (i=0; i<vs.num_alpha; i++) {
        start += num;
        if (i==vs.num_alpha-1)
            num = seqLen - start; // last variable group
        else
            num = vs.alpha[i]*seqLen;
        siteNum[i+1] = num;
        for (j=start; j<start+num && j<seqLen; j++) {
            siteCat[j] = i+1;
        }
    }

    // generate the corresponding sites of the root sequence
    currSeq = &(interSeqs[(numInterNodes-1)*seqLen]);
    start = 0;
    for (i=0; i<vs.num_alpha+1; i++) {
        num = siteNum[i];
        if (i==0) {
            // constant sites
            genSeq(currSeq, start, num, vs.probXGivenInv, numChars);
        } else {
            genSeq(currSeq, start, num, &(vs.rootNodeFreq[(i-1)*numChars]), numChars);
        }
        start += num;
    }
    
    // start simulations of the sequences for each internal and terminal node
    for (i=numLineTopMat-1; i>=0; i--) {
        parentSeq = &(interSeqs[i*seqLen]);
        /*
        // for testing only, all are A's
        parentSeq = &(interSeqs[(numInterNodes-1)*seqLen]);
         */
        
        for (j=0; j<2; j++) {
            // ========================================
            // consider the left node or the right node
            // ========================================
            int currEdgeID = i*2+j;
            int currNode = topMatrix[i*2+j];
            if (currNode < 0) {
                // leaf node
                int leafID = -currNode-1;
#ifdef DEBUG_MODE
                cout << "leafID : " << leafID << endl;
#endif
                currSeq = &(leafSeqs[leafID*seqLen]);
            } else {
                // internode
                int interNodeID = currNode-1;
#ifdef DEBUG_MODE
                cout << "interNodeID : " << interNodeID << endl;
#endif
                currSeq = &(interSeqs[interNodeID*seqLen]);
            }
            
            start = 0;
            for (k=0; k<=vs.num_alpha; k++) {
                if (k==0) {
                    // constant sites
                    for (l=start; l<start+siteNum[k]; l++) {
                        currSeq[l] = parentSeq[l];
                    }
                } else {
                    // variable sites
                     double* curr_prob = &(((allPS.ps)[k-1])->allCondProb[currEdgeID*numChars*numChars]);
                     genSeqFrParent(parentSeq, currSeq, start, siteNum[k], curr_prob, numChars);
                }
                start += siteNum[k];
            }
            
#ifdef DEBUG_MODE
            // show the distribution of nucl of parent sequence
            cout << "distribution of nucl of parent" << endl;
            int dist[4];
            for (k=0; k<4; k++)
                dist[k] = 0;
            for (k=0; k<seqLen; k++) {
                if (k>0 && siteCat[k-1]!=siteCat[k]) {
                    // print out the distribution
                    int tot = 0;
                    for (int p=0; p<4; p++) {
                        tot+=dist[p];
                    }
                    for (int p=0; p<4; p++) {
                        if (p>0)
                            cout << ",";
                        printf("%7.5f", (double)dist[p]/tot);
                    }
                    cout << " ";
                    // reset the numbers to zeros
                    for (int p=0; p<4; p++) {
                        dist[p]=0;
                    }
                }
                dist[parentSeq[k]]++;
            }
            // print out the distribution
            int tot = 0;
            for (int p=0; p<4; p++) {
                tot+=dist[p];
            }
            for (int p=0; p<4; p++) {
                if (p>0)
                    cout << ",";
                printf("%7.5f", (double)dist[p]/tot);
            }
            cout << endl;

            // show the distribution of nucl of current sequence
            cout << "distribution of nucl of current sequence" << endl;
            for (k=0; k<4; k++)
                dist[k] = 0;
            for (k=0; k<seqLen; k++) {
                if (k>0 && siteCat[k-1]!=siteCat[k]) {
                    // print out the distribution
                    int tot = 0;
                    for (int p=0; p<4; p++) {
                        tot+=dist[p];
                    }
                    for (int p=0; p<4; p++) {
                        if (p>0)
                            cout << ",";
                        printf("%7.5f", (double)dist[p]/tot);
                    }
                    cout << " ";
                    // reset the numbers to zeros
                    for (int p=0; p<4; p++) {
                        dist[p]=0;
                    }
                }
                dist[currSeq[k]]++;
            }
            // print out the distribution
            tot = 0;
            for (int p=0; p<4; p++) {
                tot+=dist[p];
            }
            for (int p=0; p<4; p++) {
                if (p>0)
                    cout << ",";
                printf("%7.5f", (double)dist[p]/tot);
            }
            cout << endl;
            
            
            // show the distribution of the conversion
            int conv[16];
            for (k=0; k<16; k++)
                conv[k]=0;
            for (k=0; k<seqLen; k++) {
                if (k>0 && siteCat[k-1]!=siteCat[k]) {
                    // print out the distribution
                    for (int p=0; p<4; p++) {
                        int tot = 0;
                        for (int q=0; q<4; q++) {
                            tot+=conv[p*4+q];
                        }
                        for (int q=0; q<4; q++) {
                            if (q>0)
                                cout << ",";
                            printf("%7.5f", (double)conv[p*4+q]/tot);
                        }
                        cout << endl;
                    }
                    cout << endl;
                    // reset the numbers to zeros
                    for (int p=0; p<4; p++) {
                        for (int q=0; q<4; q++) {
                            conv[p*4+q] = 0;
                        }
                    }
                }
                conv[parentSeq[k]*4+currSeq[k]]++;
            }
            // print out the distribution
            for (int p=0; p<4; p++) {
                int tot = 0;
                for (int q=0; q<4; q++) {
                    tot+=conv[p*4+q];
                }
                for (int q=0; q<4; q++) {
                    if (q>0)
                        cout << ",";
                    printf("%7.5f", (double)conv[p*4+q]/tot);
                }
                cout << endl;
            }
            cout << endl;
#endif
        }
    }
    
    /*
    // show the resulting sequences
#ifdef DEBUG_MODE
    int dist[4];
    int sum;
    cout << "nucl distribution of internal nodes" << endl;
    for (i=0; i<numInterNodes; i++) {
        cout << "[" << i << "] ";
        dist[0]=0;dist[1]=0;dist[2]=0;dist[3]=0;sum=0;
        for (j=0; j<seqLen; j++) {
            if (j>0 && siteCat[j]!=siteCat[j-1]) {
                cout << "(";
                for (n=0;n<4;n++) {
                    if (n>0)
                        cout << ",";
                    printf("%7.5f", (double) dist[n]/sum);
                }
                cout << ")";
                cout << " ";
                dist[0]=0;dist[1]=0;dist[2]=0;dist[3]=0;sum=0;
            }
            // cout << int2char[interSeqs[i*seqLen+j]];
            dist[interSeqs[i*seqLen+j]]++;sum++;
        }
        cout << "(";
        for (n=0;n<4;n++) {
            if (n>0)
                cout << ",";
            printf("%7.5f", (double) dist[n]/sum);
        }
        cout << ")";
        cout << endl;
    }
    
    cout << "nucl distribution of the leaf nodes" << endl;
    for (i=0; i<numSpecies; i++) {
        cout << "[" << i << "] ";
        dist[0]=0;dist[1]=0;dist[2]=0;dist[3]=0;sum=0;
        for (j=0; j<seqLen; j++) {
            if (j>0 && siteCat[j]!=siteCat[j-1]) {
                cout << "(";
                for (n=0;n<4;n++) {
                    if (n>0)
                        cout << ",";
                    printf("%7.5f", (double) dist[n]/sum);
                }
                cout << ")";
                cout << " ";
                dist[0]=0;dist[1]=0;dist[2]=0;dist[3]=0;sum=0;
            }
            // cout << int2char[leafSeqs[i*seqLen+j]];
            dist[leafSeqs[i*seqLen+j]]++;sum++;
        }
        cout << "(";
        for (n=0;n<4;n++) {
            if (n>0)
                cout << ",";
            printf("%7.5f", (double) dist[n]/sum);
        }
        cout << ")";
        cout << endl;
    }
#endif
     */
     
     
     // update the parameters such that
     // S6 in the rate matrix is set to 1
     allPS.updateContent(2);

     // print out the content of parameters
     string statusFileName = outPrefix + ".stat";
     ofstream statOut;
     statOut.open((char*) statusFileName.c_str());
     if (!statOut.is_open()) {
        cerr << "Error! The output status file : " << statusFileName << " cannot be created" << endl;
        exit(1);
     }
     vs.showContent(statOut);
     string str = "";
     allPS.showContent(str, topMatrix, &leafList);
     statOut << str;

     // update the parameters such that
     // Edge length is set to the rate of substitution
     allPS.updateContent(1);
     // update the edgeLen vector
     for (i=0; i<numCategories; i++) {
         for (j=0; j<numEdges; j++) {
            edgeLens[i*numEdges + j] = (allPS.ps[i])->t[j];
         }
     }

     string topStr = topMatrixToTreeFormat(topMatrix, &leafList, edgeLens, numCategories);
     statOut << topStr << endl;
     statOut.close();
        
    
     // generate an integer array from 0 to seqLen-1 in random order
    int* randIntArray = genRandIntArray(seqLen);

    /*
     // print out the integer array
     for (i=0; i<seqLen; i++)
         cerr << randIntArray[i] << endl;
     */
     // output the alignment of the leaves
     string outputFileName = outPrefix + ".out";
     ofstream fout;
     fout.open((char*) outputFileName.c_str());
     if (!fout.is_open()) {
         cerr << "Error! The output file : " << outputFileName << " cannot be created" << endl;
         exit(1);
     }
     
     int c;
     int space_len;
     if (outFormat==2) {
         fout << " " << numSpecies << " " << seqLen << endl;
         for (i=0; i<numSpecies; i++) {
             space_len = 10 - leafList[i].length();
             fout << leafList[i];
             for (j=0; j<space_len; j++)
                fout << " ";
             for (j=0; j<seqLen; j++) {
                c = randIntArray[j];
                fout << int2char[leafSeqs[i*seqLen+c]];
             }
             fout << endl;
         }
     } else {
         for (i=0; i<numSpecies; i++) {
             fout << ">" << leafList[i] << endl;
             for (j=0; j<seqLen; j++) {
                c = randIntArray[j];
                fout << int2char[leafSeqs[i*seqLen+c]];
             }
             fout << endl;
         }
     }
     fout.close();
     
     cout << endl;
     cout << "Output file: " << outputFileName << endl;
     cout << "Hetero 2 - Version " << VERSION << " finished" << endl;
    
    delete[] topMatrix;
    delete[] interSeqs;
    delete[] leafSeqs;
    delete[] siteNum;
    delete[] siteCat;
    delete[] randIntArray;
}
