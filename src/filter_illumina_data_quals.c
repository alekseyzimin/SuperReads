#include<stdio.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
  size_t nbytes = 10000000;
  ssize_t bytes_read,sequence_len,qual_len;
  int i;
  char *my_string;
  char *sequence;
  char *quals;
  
  my_string = malloc (nbytes);
  sequence = malloc(nbytes);
  quals=malloc(nbytes);
  
  if(argc==1)
    {
      printf("Need one agrument -- zero quality character!\n");
      return(1);
    }

  while(1)
    {
      bytes_read = getline (&my_string, &nbytes, stdin);
      if(bytes_read==-1)
	break;
      else
	{
	  if(my_string[0]=='@')
	    {
	      sequence_len = getline (&sequence, &nbytes, stdin);
	      qual_len = getline (&quals, &nbytes, stdin);
	      qual_len = getline (&quals, &nbytes, stdin);
	      if(sequence_len!=qual_len)
		{
		  printf("Number of bases and number of quals are different in read %s!!!",my_string+1);
		  return(1);
		}
	      else
		{
		  for(i=0;i<sequence_len-1;i++)
		    {
		      
		      if(quals[i]<=argv[1][0])
			{
			  sequence[i]='N';
			}
		    }
		  printf(">%s%s",my_string+1,sequence);
		}
	    }
	}
    }
return(0);
}

