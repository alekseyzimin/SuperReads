#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <charb.hpp>

int main(int argc, char *argv[])
{
  ssize_t bytes_read,i,j;
  charb my_string(100000),rev_sequence(100000);
  
  while(1)
    {
      if(fgets (my_string, 1, stdin)==NULL)
	break;
      else
	{
	  if(my_string[0]=='>')
	    {
	      printf("%s",(char*)my_string);
	    }
	  else
	    {
             bytes_read=strlen(my_string);
	      j=0;
	      for(i=bytes_read-2;i>=0;i--)
		{
		  switch (my_string[i])
		    {
		    case 'A': rev_sequence[j]='T';break;
		    case 'C': rev_sequence[j]='G';break;
		    case 'G': rev_sequence[j]='C';break;
		    case 'T': rev_sequence[j]='A';break;
		    case 'a': rev_sequence[j]='t';break;
		    case 'c': rev_sequence[j]='g';break;
		    case 'g': rev_sequence[j]='c';break;
		    case 't': rev_sequence[j]='a';break;
		    default:  rev_sequence[j]='N';
		    }
		  j++;
		}
              rev_sequence[j]='\0';
	      printf("%s\n",(char*)rev_sequence);
	    }
	}
    }
return(0);
}

