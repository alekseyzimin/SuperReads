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
     if (argc < 2) {
	  fprintf (stderr, "You must have a unitig_layout.txt file as the input argument. Bye!\n");
	  exit (1); }
     char *unitigLayoutFile = argv[1];

     FILE *infile;
     charb line(100);

     infile = fopen ("adjustmentFactorsForGCPct.txt", "r");
     std::vector<char *> flds;
     std::vector<double> dataPoints;
     while (fgets (line, 100, infile)) {
	  getFldsFromLine (line, flds);
	  dataPoints.push_back (atof (flds[1])); }
     fclose (infile);
	  
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

	  double ratio = calculateValueForPolynomial (cg_content, dataPoints);
	  double retVal = 1 / ratio;
	  
	  printf ("%d %f %lf\n", unitigNumber, ratio, cg_content);
     }
     fclose (infile);

     return (0);
}

double calculateValueForPolynomial (double arg, std::vector<double> &dataPoints)
{
     double value = 0;
     double arg100 = arg * 100;
     int iarg100 = arg100;
     double low = arg100 - iarg100;
     double high = 1 - low;
     if (iarg100 == 100)
	  value = dataPoints[100];
     else
	  value = dataPoints[iarg100] * high + dataPoints[iarg100+1] * low;

     return (value);
}

