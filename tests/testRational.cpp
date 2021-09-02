// Test of the Rational class


#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>

#include "rational.h"
#include "spxdefines.h"
#include "spxalloc.h"
#include "spxout.h"

#include "gmp.h"
#include "gmpxx.h"

using namespace soplex;


bool isEqual(const double a, const double b)
{
   double eps = 1e-6;
   bool result = false;
   if ( abs(a - b) < eps)
      result = true;
   return result;
}

int main() {

   std::cout << std::setprecision(30);

   double testTypecast[4];
   testTypecast[0] = 3452344.344235;
   testTypecast[1] = -0.010240;
   testTypecast[2] = 3.1e-12;
   testTypecast[3] = -0.933524e20;

   // tests if construction from double and typecasting works. Not sure how to test all functions separately.
   Rational testRational;
   for (int i = 0; i < 4; i++)
   {
      testRational = testTypecast[i];
      if ( Real(testRational) != testTypecast[i] )
      {
         std::cout << "There seems to be a problem with typecasting or constructing. Test double : " << testTypecast[i] << std::endl;
         return 0;
      }
   }

   int errorCounterConstructor = 0;

   // tests the copy constructor
   Rational testCopyRational(testRational);
   if ( Real(testCopyRational) != Real(testRational) )
   {
      std::cout << "There seems to be a problem with the copy constructor. " << Real(testCopyRational) << " constructed from " << Real(testRational) << std::endl;
      errorCounterConstructor++;
   }

   // tests the copy assignment operator
   Rational testCopyAssignmentRational;
   testCopyAssignmentRational = testRational;
   if ( Real(testCopyAssignmentRational) != Real(testRational) )
   {
      std::cout << "There seems to be a problem with the copy assignment operator. " << Real(testCopyAssignmentRational) << " assigned from " << Real(testRational) << std::endl;
      errorCounterConstructor++;
   }
   testRational = 42;
   if ( Real(testCopyAssignmentRational) == 42 )
   {
      std::cout << "There seems to be a problem with the copy assignment operator. The Rational changed when the one it was assigned from changed." << std::endl;
      errorCounterConstructor++;
   }

   // tests construction from int
   double testFromInt[4];
   testFromInt[0] = 134;
   testFromInt[1] = 45634613233;
   testFromInt[2] = -13445;
   testFromInt[3] = 0;
   for (int i = 0; i < 4; i++)
   {
      testRational = testFromInt[i];
      if ( Real(testRational) != testFromInt[i] )
      {
         std::cout << "There seems to be a problem with construction from int. Test int : " << testFromInt[i] << std::endl;
         errorCounterConstructor++;
      }
   }

   // tests construction from long double
   long double testFromLDouble[4];
   testFromLDouble[0] = 3452344.344235;
   testFromLDouble[1] = -0.0000010240;
   testFromLDouble[2] = 3.1e-12;
   testFromLDouble[3] = -0.933524e20;
   for (int i = 0; i < 4; i++)
   {
      testRational = testFromLDouble[i];
      if ( Real(testRational) != testFromLDouble[i] )
      {
         std::cout << "There seems to be a problem with construction from long double. Test long double : " << testTypecast[i] << std::endl;
         errorCounterConstructor++;
      }
   }

#if 0
   // tests construction from mpq_t
   mpq_t fromMpqt;
   double fromMpqtD = 21.44242;
   mpq_init(fromMpqt);
   mpq_set_d(fromMpqt, fromMpqtD);
   Rational fromMpqtR(fromMpqt);
   mpq_clear(fromMpqt);
   if ( Real(fromMpqtR) != fromMpqtD )
   {
      std::cout << "There seems to be a problem with construction from mpq_t. Test mpq_t = " << fromMpqtD << std::endl;
      errorCounterConstructor++;
   }
#endif

#if 0
   // tests construction from mpq_t
   mpq_t fromMpqt;
   double fromMpqtD = 21.44242;
   mpq_init(fromMpqt);
   mpq_set_d(fromMpqt, fromMpqtD);
   Rational fromMpqtR(fromMpqt);
   mpq_clear(fromMpqt);
   if ( Real(fromMpqtR) == fromMpqtD )
   {
      std::cout << "There seems to be a problem with construction from mpq_t. Test mpq_t = " << fromMpqtD << std::endl;
      errorCounterConstructor++;
   }

   // tests construction from mqp_class
   double fromMpqClassD = 1324.223;
   mpq_class fromMpqClass(fromMpqClassD);
   Rational fromMpqClassR(fromMpqClass);
   if ( Real(fromMpqClassR) == fromMpqClassD )
   {
      std::cout << "There seems to be a problem with construction from mpq_class. Test mpq_t = " << fromMpqClassD << std::endl;
      errorCounterConstructor++;
   }

#endif

   if( errorCounterConstructor == 0 )
   {
      std::cout << "No errors found with the Constructor functions and Typecasting to Real. Further tests are based on the assumption that they work correctly." << std::endl;
   }



   // tests the arithmetic operators
   double ADoubleOne[4];
   ADoubleOne[0] = 42;
   ADoubleOne[1] = -0.33373;
   ADoubleOne[2] = 401.001;
   ADoubleOne[3] = -234.21;

   double ADoubleTwo[4];
   ADoubleTwo[0] = 24;
   ADoubleTwo[1] = 0.45551;
   ADoubleTwo[2] = -12.2441;
   ADoubleTwo[3] = 1001;
   int errorCounterA = 0;
   Rational ARationalOne;
   Rational ARationalTwo;
   Rational AZero(0);

   for (int i = 0; i < 4 ; i++)
   {
      // addition
      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne += ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] + ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the += $Rational operator. " << ADoubleOne[i] << " += " << ADoubleTwo[i] << " results in: " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne += ADoubleTwo[i];
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] + ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the += $Real operator. " << ADoubleOne[i] << " += " << ADoubleTwo[i] << " results in: " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ARationalOne + ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] + ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $Rational + $Rational operator. " << ADoubleOne[i] << " + " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ADoubleOne[i] + ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] + ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $double + $Rational operator. " << ADoubleOne[i] << " + " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ARationalOne + ADoubleTwo[i];
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] + ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $Rational + $double operator. " << ADoubleOne[i] << " + " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      // subtraction
      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne -= ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] - ADoubleTwo[i]) )

      {
         std::cout << "There seems to be a problem with the -= $Rational operator. " << ADoubleOne[i] << " -= " << ADoubleTwo[i] << " results in: " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne -= ADoubleTwo[i];
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] - ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the -= $Real operator. " << ADoubleOne[i] << " -= " << ADoubleTwo[i] << " results in: " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ARationalOne - ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] - ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $Rational - $Rational operator. " << ADoubleOne[i] << " - " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ADoubleOne[i] - ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] - ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $double - $Rational operator. " << ADoubleOne[i] << " - " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ARationalOne - ADoubleTwo[i];
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] - ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $Rational - $double operator. " << ADoubleOne[i] << " - " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      // multiplication
      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne *= ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] * ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the *= $Rational operator. " << ADoubleOne[i] << " *= " << ADoubleTwo[i] << " results in: " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne *= ADoubleTwo[i];
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] * ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the *= $Real operator. " << ADoubleOne[i] << " *= " << ADoubleTwo[i] << " results in: " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ARationalOne * ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] * ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $Rational * $Rational operator. " << ADoubleOne[i] << " * " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ADoubleOne[i] * ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] * ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $double * $Rational operator. " << ADoubleOne[i] << " * " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ARationalOne * ADoubleTwo[i];
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] * ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $Rational * $double operator. " << ADoubleOne[i] << " * " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      // division
      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne /= ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] / ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the /= $Rational operator. " << ADoubleOne[i] << " /= " << ADoubleTwo[i] << " results in: " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne /= ADoubleTwo[i];
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] / ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the /= $Real operator. " << ADoubleOne[i] << " /= " << ADoubleTwo[i] << " results in: " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ARationalOne / ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] / ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $Rational / $Rational operator. " << ADoubleOne[i] << " / " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ADoubleOne[i] / ARationalTwo;
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] / ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $double / $Rational operator. " << ADoubleOne[i] << " / " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      ARationalTwo = ADoubleTwo[i];
      ARationalOne = ARationalOne / ADoubleTwo[i];
      if ( !isEqual(Real(ARationalOne), ADoubleOne[i] / ADoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with the $Rational / $double operator. " << ADoubleOne[i] << " / " << ADoubleTwo[i] << " = " << Real(ARationalOne) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      if ( !isEqual(Real(-ARationalOne), -ADoubleOne[i]) )
      {
         std::cout << "There seems to be a problem with Negation operator. -" << Real(ARationalOne) << " = " << Real(-ARationalOne) << std::endl;
         errorCounterA++;
      }

      if ( !isEqual(Real(abs(ARationalOne)), abs(ADoubleOne[i])) )
      {
         std::cout << "There seems to be a problem with the absolute function. abs(" << Real(ARationalOne) << ") = " << Real(abs(ARationalOne)) << std::endl;
         errorCounterA++;
      }

      ARationalOne = ADoubleOne[i];
      if ( !isEqual(sign(ARationalOne), ADoubleOne[i]/abs(ADoubleOne[i])) )
      {
         std::cout << "There seems to be a problem with the sign function. sign(" << ADoubleOne[i] <<") = " << sign(ARationalOne) << std::endl;
         errorCounterA++;
      }

      if ( sign(AZero) != 0 )
      {
         std::cout << "There seems to be a problem with the sign function. sign(0) = " << sign(AZero) << std::endl;
         errorCounterA++;
      }
   }

   if( errorCounterA == 0 )
   {
      std::cout << "No errors found in the arithmetic functions." << std::endl;
   }

   // tests the relational operators
   double RDoubleOne[4];
   RDoubleOne[0] = -12;
   RDoubleOne[1] = -23;
   RDoubleOne[2] = 234;
   RDoubleOne[3] = 44.123;
   double RDoubleTwo[4];
   RDoubleTwo[0] = 311;
   RDoubleTwo[1] = -23;
   RDoubleTwo[2] = 0.4;
   RDoubleTwo[3] = 44.123;
   Rational RRationalOne;
   Rational RRationalTwo;
   int errorCounterR = 0;

   for (int i=0; i < 4; i++)
   {
      RRationalOne = RDoubleOne[i];
      RRationalTwo = RDoubleTwo[i];

      // equality
      if ( (RRationalOne == RRationalTwo) != (RDoubleOne[i] == RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational == $Rational operator. " << RDoubleOne[i] << " == " << RDoubleTwo[i] << " evaluates to " << (RRationalOne == RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RDoubleOne[i] == RRationalTwo) != (RDoubleOne[i] == RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $double == $Rational operator. " << RDoubleOne[i] << " == " << RDoubleTwo[i] << " evaluates to " << (RDoubleOne[i] == RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RRationalOne == RDoubleTwo[i]) != (RDoubleOne[i] == RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational == $double operator. " << RDoubleOne[i] << " == " << RDoubleTwo[i] << " evaluates to " << (RRationalOne == RDoubleTwo[i]) << std::endl;
         errorCounterA++;
      }

      // inequality
      if ( (RRationalOne != RRationalTwo) != (RDoubleOne[i] != RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational != $Rational operator. " << RDoubleOne[i] << " != " << RDoubleTwo[i] << " evaluates to " << (RRationalOne != RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RDoubleOne[i] != RRationalTwo) != (RDoubleOne[i] != RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $double != $Rational operator. " << RDoubleOne[i] << " != " << RDoubleTwo[i] << " evaluates to " << (RDoubleOne[i] != RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RRationalOne != RDoubleTwo[i]) != (RDoubleOne[i] != RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational != $double operator. " << RDoubleOne[i] << " != " << RDoubleTwo[i] << " evaluates to " << (RRationalOne != RDoubleTwo[i]) << std::endl;
         errorCounterA++;
      }

      // greater than
      if ( (RRationalOne > RRationalTwo) != (RDoubleOne[i] > RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational > $Rational operator. " << RDoubleOne[i] << " > " << RDoubleTwo[i] << " evaluates to " << (RRationalOne > RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RDoubleOne[i] > RRationalTwo) != (RDoubleOne[i] > RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $double > $Rational operator. " << RDoubleOne[i] << " > " << RDoubleTwo[i] << " evaluates to " << (RDoubleOne[i] > RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RRationalOne > RDoubleTwo[i]) != (RDoubleOne[i] > RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational > $double operator. " << RDoubleOne[i] << " > " << RDoubleTwo[i] << " evaluates to " << (RRationalOne > RDoubleTwo[i]) << std::endl;
         errorCounterA++;
      }

      // greater than or equal to
      if ( (RRationalOne >= RRationalTwo) != (RDoubleOne[i] >= RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational >= $Rational operator. " << RDoubleOne[i] << " >= " << RDoubleTwo[i] << " evaluates to " << (RRationalOne >= RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RDoubleOne[i] >= RRationalTwo) != (RDoubleOne[i] >= RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $double >= $Rational operator. " << RDoubleOne[i] << " >= " << RDoubleTwo[i] << " evaluates to " << (RDoubleOne[i] >= RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RRationalOne >= RDoubleTwo[i]) != (RDoubleOne[i] >= RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational >= $double operator. " << RDoubleOne[i] << " >= " << RDoubleTwo[i] << " evaluates to " << (RRationalOne >= RDoubleTwo[i]) << std::endl;
         errorCounterA++;
      }

      // lesser than
      if ( (RRationalOne < RRationalTwo) != (RDoubleOne[i] < RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational < $Rational operator. " << RDoubleOne[i] << " < " << RDoubleTwo[i] << " evaluates to " << (RRationalOne < RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RDoubleOne[i] < RRationalTwo) != (RDoubleOne[i] < RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $double < $Rational operator. " << RDoubleOne[i] << " < " << RDoubleTwo[i] << " evaluates to " << (RDoubleOne[i] < RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RRationalOne < RDoubleTwo[i]) != (RDoubleOne[i] < RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational < $double operator. " << RDoubleOne[i] << " < " << RDoubleTwo[i] << " evaluates to " << (RRationalOne < RDoubleTwo[i]) << std::endl;
         errorCounterA++;
      }

      // lesser than or equal to
      if ( (RRationalOne <= RRationalTwo) != (RDoubleOne[i] <= RDoubleTwo[i]) )

      {
         std::cout << "There seems to be a problem with $Rational <= $Rational operator. " << RDoubleOne[i] << " <= " << RDoubleTwo[i] << " evaluates to " << (RRationalOne <= RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RDoubleOne[i] <= RRationalTwo) != (RDoubleOne[i] <= RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $double <= $Rational operator. " << RDoubleOne[i] << " <= " << RDoubleTwo[i] << " evaluates to " << (RDoubleOne[i] <= RRationalTwo) << std::endl;
         errorCounterA++;
      }

      if ( (RRationalOne <= RDoubleTwo[i]) != (RDoubleOne[i] <= RDoubleTwo[i]) )
      {
         std::cout << "There seems to be a problem with $Rational <= $double operator. " << RDoubleOne[i] << " <= " << RDoubleTwo[i] << " evaluates to " << (RRationalOne <= RDoubleTwo[i]) << std::endl;
         errorCounterA++;
      }
   }

   if( errorCounterR == 0 )
   {
      std::cout << "No errors found in the relational operators." << std::endl;
   }

   // tests rationalToString
   double RtoSDoubleIn[4];
   RtoSDoubleIn[0] = 0.335e12;
   RtoSDoubleIn[1] = 42;
   RtoSDoubleIn[2] = -2312.01;
   RtoSDoubleIn[3] = -0.22e-21;
   double RtoSDoubleOut;
   Rational RtoSRational;
   std::stringstream RtoSStream;
   int errorCounterRtoS = 0;

   for (int i = 0; i < 4; i++)
   {
      RtoSRational = RtoSDoubleIn[i];
      RtoSStream << RtoSRational.str();
      RtoSStream >> RtoSDoubleOut;
      RtoSStream.clear();
      if( !isEqual(RtoSDoubleIn[i], RtoSDoubleOut) )
      {
         std::cout << "There seems to be a problem with the rationalToString function. rationalToString(" << RtoSDoubleIn[i] << ") = " << RtoSDoubleOut << std::endl;
         errorCounterRtoS++;
      }
   }

   if( errorCounterRtoS == 0 )
   {
      std::cout << "No errors found in the rationalToString function." << std::endl;
   }

   // tests readString from floating point input
   Rational rSRational;
   double rSDouble;
   char* rSString[6];
   rSString[0] = "-95.2e13";
   rSString[1] = ".112";
   rSString[2] = "123.1901";
   rSString[3] = "-234e23";
   int errorCounterrS = 0;
   std::stringstream rSStream;

   for (int i = 0; i < 4; i++)
   {
      rSRational.readString(rSString[i]);
      rSStream << rSString[i];
      rSStream >> rSDouble;
      rSStream.clear();
      if( !isEqual(Real(rSRational), rSDouble) )
      {
         std::cout << "There seems to be a problem with the readString function. " << rSString[i] << " reads as " << Real(rSRational) << std::endl;
         errorCounterrS++;
      }
   }

   // tests readString from Rational formated input
   rSString[4] = "46/124";
   rSString[5] = "-21/490";

   rSRational.readString(rSString[4]);
   rSDouble = 46.0/124.0;
   if ( !isEqual(Real(rSRational), rSDouble) )
   {
      std::cout << "There seems to be a problem with the readString function. " << rSString[4] << " reads as " << Real(rSRational) << std::endl;
      errorCounterrS++;
   }

   rSRational.readString(rSString[5]);
   rSDouble = -21.0/490.0;
   if ( !isEqual(Real(rSRational), rSDouble) )
   {
      std::cout << "There seems to be a problem with the readString function. " << rSString[5] << " reads as " << Real(rSRational) << std::endl;
      errorCounterrS++;
   }

   if( errorCounterrS == 0 )
   {
      std::cout << "No errors found in the readString function." << std::endl;
   }


   // tests readString from floating point input
   int errorCounterrSR = 0;

   for (int i = 0; i < 4; i++)
   {
      readStringRational(rSString[i], rSRational);
      rSStream << rSString[i];
      rSStream >> rSDouble;
      rSStream.clear();
      if( !isEqual(Real(rSRational), rSDouble) )
      {
         std::cout << "There seems to be a problem with the readStringRational function. " << rSString[i] << " reads as " << Real(rSRational) << std::endl;
         errorCounterrSR++;
      }
   }

   // tests readString from Rational formated input
   readStringRational(rSString[4], rSRational);
   rSDouble = 46.0/124.0;
   if ( !isEqual(Real(rSRational), rSDouble) )
   {
      std::cout << "There seems to be a problem with the readStringRational function. " << rSString[4] << " reads as " << Real(rSRational) << std::endl;
      errorCounterrSR++;
   }

   readStringRational(rSString[5], rSRational);
   rSDouble = -21.0/490.0;
   if ( !isEqual(Real(rSRational), rSDouble) )
   {
      std::cout << "There seems to be a problem with the readStringRational function. " << rSString[5] << " reads as " << Real(rSRational) << std::endl;
      errorCounterrSR++;
   }

   if( errorCounterrSR == 0 )
   {
      std::cout << "No errors found in the readStringRational function." << std::endl;
   }

}
