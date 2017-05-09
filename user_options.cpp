/*
 *
 * user_options.cpp
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


#include "user_options.h"

// read the arguments and collect the options
void GetOptions::read(int argc, char** argv) {
    
    int i;
    char* currArg;
    int currLen;
    char flag;
    string value;
    for (i=1; i<argc; i++) {
        // skip the first argument
        currArg = argv[i];
        currLen = (int) strlen(currArg);
        if (currArg[0]=='-' && currLen>1) {
            // get the flag
            flag=currArg[1];
            // get the value
            if (currLen>2) {
                value = (&currArg[2]);
            } else if (argc > i+1 && argv[i+1][0]!='-') {
                i++;
                value = argv[i];
            } else {
                value = "";
            }
            flags.push_back(flag);
            values.push_back(value);
        }
    }
    
}

int GetOptions::size() {
    
    return (int) flags.size();
}

// constructor
UserOptions::UserOptions() {
    //initialize the variables
    
    // outputFormat        : the format of the output multiple sequence alignment file
    //                       1 - FASTA format
    //                       2 - sequential PHYLIP format
    //                       (default: 2)
    outputFormat = 2;
    
    // seqLen              : the length of sequences to be simulated
    seqLen = DEFAULT_SEQ_LEN;

    // edgeRepresent       : edge length representation
    //                       1 - Average number of substitutions per site
    //                       2 - Time
    //                       (default: 1)
    edgeRepresent = 1;
    
    perfectSeq = 0;

}

// To show the usage of the program
void UserOptions::outputUsage(char* progName) {
    cout << "================================================================================" << endl;
    cout << "                  Welcome to Hetero 2 version " << VERSION << endl << endl;

    cout << "Syntax:" << endl;
    cout << "  " << progName << " <tree file> <site info file> <param file list> <other options>" << endl;
    cout << "  " << progName << " -h" << endl << endl;

    cout << "  <tree file>           : \"Tree file\" lists the tree of each site category, the" << endl;
    cout << "                          edge lengths and the labels of terminal/internal nodes" << endl << endl;

    cout << "  <site info file>      : \"Site info file\" lists the detailed information of" << endl;
    cout << "                          each site category, including the site proportion and" << endl;
    cout << "                          the nucleotide distribution at the root" << endl << endl;

    cout << "  <param file list>     : \"Param file list\" shows the name of the parameter file" << endl;
    cout << "                          of each variant site category" << endl << endl;

    cout << "Note: one may refer to the examples to understand the format of the files" << endl << endl;

    cout << "other options:" << endl << endl;

    cout << "  -l <sequence length>  : The length of sequences to be simulated" << endl;
    cout << "                          (default: " << DEFAULT_SEQ_LEN << ")" << endl << endl;

    cout << "  -f <output format>    : The format of simulated multiple sequence alignment" << endl;
    cout << "                          1 - FASTA format" << endl;
    cout << "                          2 - Sequential PHYLIP format (default)" << endl << endl;

    cout << "  -o <output prefix>    : Prefix for output files" << endl;
    cout << "                          (default: <tree file> w/o .ext)" << endl << endl;
    
    cout << "  -e < 1 or 2 >         : Representation of tree edges" << endl;
    cout << "                          1 - Average number of substitutions per site" << endl;
    cout << "                          2 - Time" << endl;
    cout << "                          (default: 1)" << endl << endl;

    cout << "  -p                    : Perfect sequence simulation" << endl << endl;

    cout << "  -h                    : This help page" << endl << endl;
    
    cout << "Output files:" << endl;
    cout << "  <output prefix>.out   : The simulated multiple sequence alignment file" << endl << endl;

    cout << "Example:" << endl;
    cout << "  Example tree file      : trees.txt" << endl;
    cout << "  Example site info file : site_info_file.txt" << endl;
    cout << "  Example param file list: param_file_list.txt" << endl;
    cout << "  Example parameter files: parameter_1.txt, parameter_2.txt" << endl;
    cout << "  To execute the Hetero2 program using the default parameters:" << endl;
    cout << "      $ " << progName << " trees.txt site_info_file.txt param_file_list.txt" << endl << endl;

    cout << "Contact: " << CONTACTPERSON << endl;
    cout << "================================================================================" << endl;
}


// read the arguments and collect the user options
int UserOptions::readArguments(int argc, char** argv, string* errMsg) {
    
    // return:
    //    0 if the user options are valid
    //    1 if there exists error message
    //    2 if a help menu is called
    
    int i;
    bool duplicateOption = false;
    bool emptyOption = false;
    vector<string> token;
    int minArgNum = 4; // minimum number of arguments
    
    bool l_option_assigned = false;
    bool f_option_assigned = false;
    bool o_option_assigned = false;
    bool e_option_assigned = false;
    
    *errMsg = "";
    
    if (argc < minArgNum) {
        // invoke the help menu
        return 2;
    }
    
    treeFile = argv[1];
    siteIntoFile = argv[2];
    paramFileList = argv[3];
    
    // get the other options if there is any
    if (argc > minArgNum) {
        GetOptions options;
        options.read(argc, argv);
        
        for (i=0; i<options.size(); i++) {
            
            char flag = options.flags[i];
            string value = options.values[i];
            
            switch (flag) {
                    
                case 'l':
                    if (l_option_assigned)
                        duplicateOption = true;
                    else if (value.length() == 0)
                        emptyOption = true;
                    else {
                        l_option_assigned = true;
                        int valueInt = atoi(value.c_str());
                        if (valueInt <= 0)
                            *errMsg = "Error! The value for '-l' has to be greater than 0";
                        else
                            seqLen = valueInt;
                    }
                    break;

                case 'f':
                    if (f_option_assigned)
                        duplicateOption = true;
                    else if (value.length() == 0)
                        emptyOption = true;
                    else {
                        f_option_assigned = true;
                        int valueInt = atoi(value.c_str());
                        if (valueInt < 1 || valueInt > 2)
                            *errMsg = "Error! The value for '-f' has to be 1 or 2";
                        else
                            outputFormat = valueInt;
                    }
                    break;

                case 'e':
                    if (e_option_assigned)
                        duplicateOption = true;
                    else if (value.length() == 0)
                        emptyOption = true;
                    else {
                        e_option_assigned = true;
                        int valueInt = atoi(value.c_str());
                        if (valueInt < 1 || valueInt > 2)
                            *errMsg = "Error! The value for '-e' has to be 1 or 2";
                        else
                            edgeRepresent = valueInt;
                    }
                    break;

                case 'o':
                    if (o_option_assigned)
                        duplicateOption = true;
                    else {
                        o_option_assigned = true;
                        if (value.length()==0)
                            emptyOption = true;
                        else
                            outputPrefix = value;
                    }
                    break;
					
                case 'p':
                    perfectSeq = 1;
                    break;

                    case 'h':
                    return 2;
                    break;
                    
                default:
                    *errMsg = "Unknown option '-" + string(1,flag);
                    break;
                    
            } // case
            
            if (duplicateOption) {
                *errMsg = "Error! Duplicate option -" + string(1,flag);
            } else if (emptyOption) {
                *errMsg = "Error! Empty value for the option -" + string(1,flag);
            }
            
            if (*errMsg != "") {
                return 1;
            }
        }
    }
    if (outputPrefix == "") {
        setDefaultOutputPrefix();
    }
    
    return 0;
}

// by default, prefixOut = <alignment file> w/o .ext
void UserOptions::setDefaultOutputPrefix() {
    
    outputPrefix = removeExtension(treeFile);
    
}

// To show the summary of the parameters
void UserOptions::showSummary() {
    cout << "================================================================================" << endl;

    cout << "                   Welcome to Hetero 2 - Version " << VERSION << endl << endl;
    
    cout << "Tree file ............................................. " << treeFile << endl;

    cout << "Site info file ........................................ " << siteIntoFile << endl;
    
    cout << "Param file list ....................................... " << paramFileList << endl;

    cout << "Length of sequences to be simulated ................... " << seqLen << endl;

    cout << "Format of output alignment ............................ ";
    switch (outputFormat) {
        case 1:
            cout << "FASTA" << endl;
            break;
        case 2:
            cout << "sequential PHYLIP" << endl;
            break;
    }
    
    cout << "Input parameters in ................................... " << outputPrefix << ".stat" << endl;

    cout << "Output file ........................................... " << outputPrefix << ".out" << endl;
    
    cout << "================================================================================" << endl;
}

