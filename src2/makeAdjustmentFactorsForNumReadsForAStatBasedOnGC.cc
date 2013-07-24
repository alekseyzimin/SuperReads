// Takes a unitig_layout.txt file and reports coverage adjustment factors and GC
// content for each unitig (in [0,1]) to stdout; the output line is
// unitig# adjustmentFactor gcProportion
//
// You must call this with the name of the unitig_layout.txt file as an argument
//
// An example of the input file is
// /genome3/raid/tri/mouse/assembly_new_quorum-FROZEN/unitig_layout.txt
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <charb.hpp>
#include <misc.hpp>

const int windowSize = 200;

double calculateValueForPolynomial (double arg, std::vector<double> &vect);

int main (int argc, char **argv)
{
     double midVal = .55, minRange = .15, maxRange = .8;
     std::vector<double> polynomialCoefficients1 = { -3182.98379217935, -2095.35576418852, -645.53152074943, -92.1770222255863, 37.8653227363197 };
     std::vector<double> polynomialCoefficients2 = { -142.627135012101, 37.8653227363197 };
     // Used the derivative to find the x coord of the max value
     double maxVal = calculateValueForPolynomial (0.411073729214977-midVal, polynomialCoefficients1);
     if (argc < 2) {
	  fprintf (stderr, "You must have a unitig_layout.txt file as the input argument. Bye!\n");
	  exit (1); }
     char *unitigLayoutFile = argv[1];


     FILE *infile;
     charb line(100);

     infile = fopen (unitigLayoutFile, "r");
     int unitigNumber = -1; // Set to make the compiler happy; should always be in the file
     std::vector<int>GCcount;
     while (fgets (line, 100, infile)) {
	       if (strncmp (line, "unitig ", strlen("unitig ")) == 0) {
	       char *cptr = (char *)line + strlen ("unitig ");
	       unitigNumber = atoi (cptr);
	       continue; }
	  if (strncmp (line, "cns ", strlen("cns ")) != 0)
	       continue;
	  char *cptr = (char *)line + strlen ("cns ");
	  int len = strlen (cptr) - 1;
	  cptr[len] = 0;
	  GCcount.clear();
	  int numGCs = 0;
	  GCcount.push_back (numGCs);
	  for (int i=0; i<len; ++i) {
	       switch (cptr[i]) {
	       case 'c': case 'C': case 'g': case 'G': ++numGCs; break;
	       default:  break; // Do nothing
	       }
	       GCcount.push_back (numGCs);
	  }
	  double cg_content = (float) numGCs / len;
	  
	  double restrictedVal = cg_content;
	  if (restrictedVal < minRange)
	       restrictedVal = minRange;
	  if (restrictedVal > maxRange)
	       restrictedVal = maxRange;
	  double retVal;
	  if (restrictedVal <= midVal)
	       retVal = calculateValueForPolynomial (restrictedVal - midVal, polynomialCoefficients1);
	  else
	       retVal = calculateValueForPolynomial (restrictedVal - midVal, polynomialCoefficients2);
	  
	  double ratio = maxVal / retVal;
	  
//	  printf ("%lf %lf\n", retVal, ratio);
	  
	  
	  printf ("%d %f %lf\n", unitigNumber, ratio, cg_content);
     }
     fclose (infile);

     return (0);
}

double calculateValueForPolynomial (double arg, std::vector<double> &vect)
{
     double value = 0;
     for (unsigned int i=0; i<vect.size(); ++i) {
	  value *= arg;
//	  printf ("%d %lf\n", i, vect[i]);
	  value += vect[i]; }

     return (value);
}

