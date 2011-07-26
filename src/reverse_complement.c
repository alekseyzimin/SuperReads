#include<stdio.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
  int bytes_read=0,i,j;
  int nbytes=10000000;
  char *my_string;
  char *rev_sequence;
  
  my_string = malloc (nbytes);
  
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
	     for(i=0;i<nbytes;i++)
	     {
		if(my_string[i]=='\0')
		{
		bytes_read=i-2;
		break;
		}
	     }
	      rev_sequence = malloc( bytes_read+2 );
	      j=0;
	      for(i=bytes_read;i>=0;i--)
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
	      printf("%s\n",rev_sequence);
	    }
	}
    }
return(0);
}

