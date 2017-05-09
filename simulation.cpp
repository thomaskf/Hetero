/*
 *
 * simulation.cpp
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

#include "simulation.h"

#define NEARLY_ZERO 1.0e-4
#define TINY_NUMBER 1.0e-5
#define TINY_TIME 1.0e-3


bool isEqual(double x, double y) {
    return (fabs(x-y) < TINY_NUMBER);
}

Sim::Sim() {
    size=0;
    q=NULL;
}

Sim::Sim(double q[], int size, int seed) {
    this->size = size;
    this->q = new double[size*size];
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            this->q[i*size+j] = q[i*size+j];
        }
    }
    if (!validate()) {
        cerr << "The values of Q are not valid" << endl;
        exit(1);
    }
#ifdef DEBUG_MODE
    cout << "Q is valid" << endl;
#endif
    srand(seed);
    
    this->acc = new int[size*size];
}

Sim::Sim(double q[], int size) {
    this->size = size;
    this->q = new double[size*size];
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            this->q[i*size+j] = q[i*size+j];
        }
    }
    if (!validate()) {
        cerr << "The values of Q are not valid" << endl;
        exit(1);
    }
#ifdef DEBUG_MODE
    cout << "Q is valid" << endl;
#endif
    this->acc = new int[size*size];
}

// initialize (for the default constructor)
void Sim::init(int size, int seed) {
    this->size = size;
    this->q = new double[size*size];
    srand(seed);
    this->acc = new int[size*size];
}

// initialize (for the default constructor)
void Sim::init(int size) {
    this->size = size;
    this->q = new double[size*size];
    this->acc = new int[size*size];
}

// to initialize the q matrix
void Sim::set_Q(double* q) {
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            this->q[i*size+j] = q[i*size+j];
        }
    }
    if (!validate()) {
        cerr << "The values of Q are not valid" << endl;
        exit(1);
    }
#ifdef DEBUG_MODE
    cout << "Q is valid" << endl;
#endif
}


Sim::~Sim() {
    delete (q);
    delete (acc);
}

// to validate the q matrix
bool Sim::validate() {
    // qii = - sum_j (qij) where i!=j
    for (int i=0; i<size; i++) {
        double sum=0.0;
        for (int j=0; j<size; j++) {
            if (i==j)
                continue;
            sum += q[i*size+j];
        }
        if (!isEqual(q[i*size+i],-sum)) {
            cerr << "q[" << i+1 << "," << i+1 << "] <> -sum_j {q_" << i+1 << "j where j<>" << i+1 << "}" << endl;
            cerr << "q[" << i+1 << "," << i+1 << "] = " << q[i*size+i] << endl;
            cerr << "-sum_j {q_" << i+1 << "j where j<>" << i+1 << "} = " << -sum << endl;
            return false;
        }
    }
    return true;
}

// set the time
void Sim::set_time(double time) {
    
    this->time = time;
    
    int a_constant = 1000;
    time_delta = TINY_TIME;
    if (time_delta > time)
        time_delta = time;
    base = (int) (1.0/time_delta * a_constant);
    while (base*10 > RAND_MAX && time_delta<0.1*time) {
        // the value of base is too large
        time_delta = time_delta*10;
        base = 1.0/time_delta * a_constant;
    }
    while (base*10 > RAND_MAX && a_constant>1) {
        a_constant = a_constant/10;
        base = 1.0/time_delta * a_constant;
    }
    // show the values of bases, time_delta and a_constant
#ifdef DEBUG_MODE
    cout << "base: " << base << " time_delta: " << time_delta << " a_constant: " << a_constant << " time: " << time << " ratio b/w time & time_delta : " << time/time_delta << endl;
#endif
    if (base*10 > RAND_MAX) {
        cerr << "base: " << base << " is too large" << endl;
        cerr << "time_delta: " << time_delta << " a_constant: " << a_constant << endl;
        exit(1);
    }
    // build the acc table
    for (int i=0; i<size; i++) {
        double acc_p = 0.0;
        for (int j=0; j<size; j++) {
            if (i==j) {
                acc[i*size+j]=-1;
            } else {
                acc_p += q[i*size+j]*time_delta*base;
                acc[i*size+j]=(int)acc_p;
            }
        }
    }
    // print out acc table
#ifdef DEBUG_MODE
    cout << "acc table:" << endl;
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if (j>0)
                cout << ",";
            cout << acc[i*size+j];
        }
        cout << endl;
    }
#endif
    // verify the acc table
    // to check: (1) none of them is zero; (2) diff bw the entry > 0
    for (int i=0; i<size; i++) {
        int pre_acc = 0;
        for (int j=0; j<size; j++) {
            if (i==j)
                continue;
            if (q[i*size+j] > NEARLY_ZERO && acc[i*size+j] <= pre_acc) {
                cerr << "acc table is not valid" << endl;
                // print out the Q table
                cout << "q table:" << endl;
                for (int k=0; k<size; k++) {
                    for (int l=0; l<size; l++) {
                        cout << q[k*size+l] << ",";
                    }
                    cout << endl;
                }
                exit(1);
            } else {
                pre_acc = acc[i*size+j];
            }
        }
    }
}

int Sim::simulate(int nucl) {
    // input: nucl (i.e. the nucl at time 0)
    // output: the simulated nucl after time
    
    double t=0.0;
    int rand_num;
    int i;
    while (t < time) {
        rand_num = rand() % base;
        for (i=0; i<size; i++) {
            if (i==nucl)
                continue;
            if (rand_num <= acc[nucl*size+i]) {
                nucl = i;
                break;
            }
        }
        t+= time_delta;
    }
    return nucl;
}

void Sim::showContent() {
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if (j>0)
                cout << ",";
            cout << q[i*size+j];
        }
        cout << endl;
    }
}


//==================================================
// RandNuclGenerator
//==================================================

RandNuclGenerator::RandNuclGenerator() {
    size = 0;
    acc = NULL;
}

RandNuclGenerator::~RandNuclGenerator() {
    if (acc != NULL)
        delete[] acc;
}

// input: probabilities (dimension: size x size),
//        representing the probility from nucl i to nucl j (where 0 <= i,j <= size-1)
void RandNuclGenerator::init(double* prob, int size) {

    int i,j;
    double subtotal;
    this->size = size;
    this->acc = new int[size*size];
    for (i=0; i<size; i++) {
        subtotal=0.0;
        for (j=0; j<size; j++) {
            subtotal += (double)RAND_MAX * prob[i*size+j];
            this->acc[i*size+j] = (int) round(subtotal);
        }
    }
    
}


// input a nucl i, where 0 <= i <= size-1
// return a random nucl j, where 0 <= j <= size-1
// according to the probabilies provided
// assumption: the random seed has been initialized
int RandNuclGenerator::randomNucl(int nucl) {

    int i;
    int randNum = rand();
    for (i=0; i<size-1; i++) {
        if (randNum < acc[nucl*size + i])
            return i;
    }
    return size-1;
}

void RandNuclGenerator::showAcc() {
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if (j>0)
                cout << ",";
            cout << acc[i*size+j];
        }
        cout << endl;
    }
}
