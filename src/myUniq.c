#include <stdio.h>
#include <string.h>

char lines[2][10000000];

int main (int argc, char **argv)
{
     int curIndex = 0;

     fgets (lines[curIndex], 10000000, stdin);
     fputs (lines[curIndex], stdout);
     curIndex = 1 - curIndex;
     while (fgets(lines[curIndex], 10000000, stdin)) {
	  if (strcmp (lines[0], lines[1]) != 0)
	       fputs (lines[curIndex], stdout);
	  curIndex = 1 - curIndex; }
     
     return (0);
}
