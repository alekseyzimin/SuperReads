#include<stdio.h>
#include<stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
  ssize_t bytes_read,i,j;
  size_t nbytes=10000000;
  char *my_string;
  char *rev_sequence;
  
  my_string = (char *)malloc (nbytes);
  rev_sequence = (char *)malloc (nbytes);
  
  while(1)
    {
      if(fgets (my_string, nbytes, stdin)==NULL)
	break;
      else
	{
	  if(my_string[0]=='>')
	    {
	      printf("%s",my_string);
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
	      printf("%s\n",rev_sequence);
	    }
	}
    }
return(0);
}

