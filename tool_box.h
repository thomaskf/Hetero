/*
 *
 * tool_box.h
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
 
#ifndef __sim_tree__tool_box__
#define __sim_tree__tool_box__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>

using namespace std;


static void tokenizer(string seq, string separators, vector<string>* result) {
    // split the seq into many parts by "separators"
    // the vector<string> *result cannot be NULL
    result->clear();
    int startpos = (int) seq.find_first_not_of(separators);
    while (startpos != string::npos) {
        int endpos = (int) seq.find_first_of(separators, startpos);
        if (endpos != string::npos) {
            result->push_back(seq.substr(startpos, endpos-startpos));
            startpos = (int) seq.find_first_not_of(separators, endpos);
        } else {
            result->push_back(seq.substr(startpos));
            break;
        }
    }
}

static string removeSpacesInFront(string str) {
    // remove all the spaces in front
    size_t pos = 0;
    while (pos < str.length() && str[pos]==' ')
        pos++;
    if (pos < str.length())
        return str.substr(pos);
    else
        return "";
}

// return true if c is space or non-printable character
static bool isSpaceOrNonprintableChar(char c) {
    return (c <= 32 || c >= 127);
}


static void trim(string& str) {
    // remove the non-printable characters at the end of the string
    int i;
    for (i=(int)str.length()-1; i>=0; i--) {
        if (!isSpaceOrNonprintableChar(str[i]))
            break;
    }
    if (i==str.length()-1) {
        // do nothing
    } else if (i>=0) {
        // resize the string
        str.resize(i+1);
    } else {
        // empty the string
        str = "";
    }
}

static double sumProduct(double* vectorA, double* vectorB, int len) {
    // sum of the product of two vectors
    double sum = 0;
    for (int i=0; i<len; i++) {
        sum += vectorA[i]*vectorB[i];
    }
    return sum;
}

static map<string,int>* genReverseMap(vector<string>* strArray) {
    map<string,int>* newMap = new map<string,int>;
    for (int i=0; i<strArray->size(); i++) {
        newMap->insert(pair<string,int>(strArray->at(i),i));
    }
    return newMap;
}

// check the number of bits of the OS
// if the number of bits is not 64,
// then output the warning message
static void checkOSNumOfBits() {
    int numOfBits = sizeof(double)*8;
    
    if (numOfBits != 64) {
        cout << "================================================================================" << endl;
        cout << "                                    WARNING" << endl << endl;
        cout << "This is a " << numOfBits << "-bit system" << endl;
        cout << "To have a more accurate result, using a 64-bit system is highly recommended" << endl;
        cout << "================================================================================" << endl;
    }
}

static string removeExtension(string str) {
    size_t pos = str.find_last_of(".");
    if (pos != string::npos && pos > 0) {
        // found
        return str.substr(0,(int)pos);
    } else {
        return str;
    }
}


// integer to string
static string intToStr(int i) {
	static char* d2c = (char*) "0123456789"; // the array for digit to char
	if (i<0) {
		return "-" + intToStr(-i);
	} else if (i<10) {
		return string(1,d2c[i]);
	} else {
		return intToStr(i/10)+string(1,d2c[i%10]);
	}
}

// integer to string
// minLen : minimum length of the integer.
// If the number of digits < minimum length, then the integer will be displayed with leading zeros, unless the integer is zero
static string intToStr(int i, int minLen) {
	static char* d2c = (char*) "0123456789"; // the array for digit to char
	if (i<0) {
		return "-" + intToStr(-i, minLen);
	} else if (i<10) {
		if (minLen > 0)
			return string(minLen-1, '0') + string(1,d2c[i]);
		else
			return string(1,d2c[i]);
	} else {
		return intToStr(i/10, minLen-1)+string(1,d2c[i%10]);
	}
}

// double to string
static string doubleToStr(double d, double decimalPlace) {
	if (isnan(d)) {
		return "N/A";
	} else if (d < 0) {
		return "-" + doubleToStr(-d, decimalPlace);
	} else if (decimalPlace > 0) {
       d = round( d * pow(10.0, decimalPlace) ) / pow(10.0, decimalPlace);
		return intToStr((int)d) + "." + intToStr((int)((d - (int)d)*pow(10.0,decimalPlace)),decimalPlace);
	} else {
		return intToStr((int) round(d));
	}
}

static unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

#endif /* defined(__sim_tree__tool_box__) */
