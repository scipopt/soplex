/// Test of the Rational class

#include <math.h>
#include <iostream>
#include <string.h>

#include "rational.h"
//#include "rational.cpp"
#include "spxdefines.h"
//#include "spxdefines.cpp"
#include "spxalloc.h"
//#include "spxout.cpp"
#include "spxout.h"

using namespace soplex;

int main() {

int testint = 3;
double testdouble = 52.6;
long double testlongdouble= 6.35;

Rational testRational;
Rational fromint(testint);
Rational fromdouble(testdouble);
Rational fromlongdouble(testlongdouble);
Rational copyRational(fromint);
testRational = fromdouble;

bool equalityTest = (testRational == testRational);
bool doubleequalto = (testdouble == fromdouble);
bool equaltodouble = (fromdouble == testdouble);

bool greaterthanTest = (fromint > fromdouble);
bool doublegreaterthan = (testint > fromdouble);
bool greaterthandouble = (fromdouble > testint);


std::cout << "test construction from int " << testint << " : " << fromint << std::endl;
std::cout << "test construction from double " << testdouble << " : " << fromdouble << std::endl;
std::cout << "test construction from long double " << testlongdouble << " : " << fromlongdouble << std::endl;
std::cout << "test construction per copy constructor (should be " << fromint <<  ") : " << copyRational << std::endl;
std::cout << "test construction per = (should be " << fromdouble << " ) : " << testRational << std::endl;
std::cout << "test typecast to Real (should be " << testdouble << " ) : " << Real(fromdouble) << std::endl;

std::cout << "test equality operator (should be true) : " << equalityTest << std::endl;
std::cout << "test equality operator with first argument Real (should be true) : " << doubleequalto << std::endl;
std::cout << "test equality operator with second argument Real (should be true) : " << equaltodouble << std::endl;
std::cout << "test greater than operator (should be " << (testint > testdouble) << " ) :" << greaterthanTest << std::endl;
std::cout << "test greater than operator with first argument Real (should be " << (testint > testdouble) << " ) :" << doublegreaterthan << std::endl;
std::cout << "test greater than operator with second argument Real(should be " << (testdouble > testint) << " ) :" << greaterthandouble << std::endl;

testRational += testRational;

std::cout << "test += (should be " << testdouble + testdouble << " ) > " << testRational << " = " << Real(testRational) << std::endl;

testRational = testRational + testRational;

std::cout << "test + (should be " << 4*testdouble << " ) : " << testRational << " = " << Real(testRational) << std::endl;

fromint += 2.5;

std::cout << "test += with Real argument (should be " << testint + 2.5 << " ) : " << fromint << std::endl;

fromint = 4.5 + fromint;

std::cout << "test + with Real first argument (should be " << testint + 7 << " ) : " << fromint << std::endl;

fromint = fromint + 3.5;

std::cout << "test + with Real second argument (should be " << testint + 10.5 << " ) : " << fromint << std::endl;

Rational multTest(2.3);

multTest *= 3.4;

std::cout << "test *= with Real argument (should be " << 2.3 * 3.4  << " ) : " << multTest << " = " << Real(multTest) << std::endl;

multTest = 0.34 * multTest;

std::cout << "test * with Real first argument (should be " << 0.34 * 2.3 * 3.4 << " ) : " << multTest << " = " << Real(multTest) << std::endl;

multTest = multTest * 1.5;

std::cout << "test * with Real second argument (should be " << 0.34 * 2.3 * 3.4 * 1.5 << " ) : " << multTest << " = " << Real(multTest) << std::endl;



Rational fromstringone;
Rational fromstringtwo;
Rational fromstringthree;
Rational fromstringfour;

char* teststringone = "3310.";
char* teststringtwo = "+1.40015";
char* teststringthree = "3667.127e2";
char* teststringfour = "-2.2333e-12";

fromstringone.readString(teststringone);
fromstringtwo.readString(teststringtwo);
fromstringthree.readString(teststringthree);
fromstringfour.readString(teststringfour);

std::cout << "test reading " << teststringone << " from string (as member function): " << fromstringone << " = " << Real(fromstringone) << std::endl;
std::cout << "test reading " << teststringtwo << " from string (as member function): " << fromstringtwo << " = " << Real(fromstringtwo) << std::endl;
std::cout << "test reading " << teststringthree << " from string (as member function): " << fromstringthree << " = " << Real(fromstringthree) << std::endl;
std::cout << "test reading " << teststringfour << " from string (as member function): " << fromstringfour << " = " << Real(fromstringfour) << std::endl;

readStringRational(teststringone, fromstringone);
readStringRational(teststringtwo, fromstringtwo);
readStringRational(teststringthree, fromstringthree);
readStringRational(teststringfour, fromstringfour);

std::cout << "test reading " << teststringone << " from string : " << fromstringone << std::endl;
std::cout << "test reading " << teststringtwo << " from string : " << fromstringtwo << std::endl;
std::cout << "test reading " << teststringthree << " from string : " << fromstringthree << std::endl;
std::cout << "test reading " << teststringfour << " from string : " << fromstringfour << std::endl;

std::string rationalfromstring;
rationalfromstring = rationalToString(fromstringone, false);
std::cout << "test rationalToString (should be " << teststringone << " ) : " << rationalfromstring << std::endl;
rationalfromstring = rationalToString(fromstringtwo, false);
std::cout << "test rationalToString (should be " << teststringtwo << " ) : " << rationalfromstring << std::endl;
rationalfromstring = rationalToString(fromstringthree, false);
std::cout << "test rationalToString (should be " << teststringthree << " ) : " << rationalfromstring << std::endl;
rationalfromstring = rationalToString(fromstringfour, false);
std::cout << "test rationalToString (should be " << teststringfour << " ) : " << rationalfromstring << std::endl;

Rational fromdoubletwo;
double testdoubletwo = -2.5e3;
fromdoubletwo = testdoubletwo;

std::cout << "test assignment per $Rational = $double (should be " << testdoubletwo <<" ) : " << fromdoubletwo << " = " << Real(fromdoubletwo) << std::endl;

fromdoubletwo = abs(fromdoubletwo);
std::cout << "test absolute (should be " << abs(testdoubletwo) <<" ) : " << fromdoubletwo << " = " << Real(fromdoubletwo) << std::endl;

double testdoublethree = 42;
fromdoubletwo = testdoublethree;
fromdoubletwo = -fromdoubletwo;
std::cout << "test negation (should be " << -testdoublethree <<" ) : " << fromdoubletwo << " = " << Real(fromdoubletwo) << std::endl;

double testinfinity = infinity;
Rational frominfinity;
frominfinity = testinfinity;

std::cout << frominfinity << " = " << Real(frominfinity) <<std::endl;

frominfinity = -testinfinity;

std::cout << frominfinity << " = " << Real(frominfinity) <<std::endl;

Rational value = 1;
std::cout << "value = " << value << std::endl;

if( !value.readString("13.2") )
{
   std::cout << "something is wrong" << std::endl;
}

std::cout << "value = " << value << std::endl;


}
