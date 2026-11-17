// Regression tests for exact decimal and scientific-notation parsing.
#include <iostream>
#include <limits>
#include <string>

#include "soplex/rational.h"

int main()
{
   struct TestCase
   {
      const char* input;
      const char* expected;
   };
   const TestCase cases[] =
   {
      {"1e-1", "1/10"},
      {"-1E-1", "-1/10"},
      {"+1e-1", "1/10"},
      {"1.25e-2", "1/80"},
      {".125E+2", "25/2"},
      {"-12.50e+2", "-1250"},
      {"1e23", "100000000000000000000000"},
      {"1E+0", "1"},
      {"0e-10", "0"},
      {"0.1", "1/10"},
      {"1/10", "1/10"},
      {"42", "42"}
   };
   int failures = 0;

   for(const auto& test : cases)
   {
      const soplex::Rational actual = soplex::ratFromString(test.input);
      const soplex::Rational expected(test.expected);

      if(actual != expected)
      {
         std::cerr << test.input << ": expected " << expected << ", got " << actual << '\n';
         ++failures;
      }
   }

   // Compare against integer strings, independently of floating-point range.
   const std::string power = "1" + std::string(400, '0');

   if(soplex::ratFromString("1e-400") != soplex::Rational("1/" + power))
   {
      std::cerr << "1e-400 was not parsed exactly\n";
      ++failures;
   }

   // Do not continue a failing baseline into a floating-point overflow.
   if(failures != 0)
      return 1;

   if(soplex::ratFromString("1e400") != soplex::Rational(power))
   {
      std::cerr << "1e400 was not parsed exactly\n";
      return 1;
   }

   // Zero needs no large power, even for the most negative int exponent.
   const std::string zero = "0e" + std::to_string(std::numeric_limits<int>::min());

   if(soplex::ratFromString(zero.c_str()) != 0)
      return 1;

   return 0;
}
