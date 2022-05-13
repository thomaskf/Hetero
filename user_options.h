/*
 *
 * user_options.h
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

#ifndef __RAL_RAS__user_options__
#define __RAL_RAS__user_options__

#include <vector>
#include <cstring>
#include <string>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include "tool_box.h"
#include "definitions.h"

using namespace std;

class GetOptions {
public:
    vector<char> flags;
    vector<string> values;
    
    int size();
    
    void read(int argc, char** argv);
};

class UserOptions {
public:
    
    // ======================================================
    // Variables
    // ======================================================
    
    // treeFile            : the tree file
    //                       (compulsory)
    char* treeFile;
    
    // siteIntoFile        : the site info file
    //                       (compulsory)
    char* siteIntoFile;
    
    // param file list     : the list of parameter files
    //                       (compulsory)
    char* paramFileList;
    
    // outputFormat        : the format of the output multiple sequence alignment file
    //                       1 - FASTA format
    //                       2 - sequential PHYLIP format
    //                       (default: 2)
    int outputFormat;
    
    // outputPrefix        : prefix for output files
    //                       (default: <tree file w/o .ext)
    string outputPrefix;

    // seqLen              : the length of sequences to be simulated
    //                       (default: 10000)
    int seqLen;
    
    // edgeRepresent       : edge length representation
    //                       1 - Average number of substitutions per site
    //                       2 - Time
    //                       (default: 1)
    int edgeRepresent;
    
    // to simulate a perfect sequence
    int perfectSeq;
    
    // input seed number
    unsigned long seed_number;

    // =================================
    // Functions
    // =================================
    
	// constructor
	UserOptions();
	
    // To show the usage
    void outputUsage(char* progName);

    // read the arguments and collect the user options
    int readArguments(int argc, char** argv, string* errMsg);
    
    // To show the summary of the parameters
    void showSummary();
    
    // by default, prefixOut = <alignment file> w/o .ext
    void setDefaultOutputPrefix();
};


#endif /* defined(__RAL_RAS__user_options__) */
