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
      if(fgets(my_string, nbytes, stdin)==NULL)
	break;
      else
	{
	  if(my_string[0]=='@')
	    {
	      if(fgets(sequence, nbytes, stdin)==NULL)
		{
		return(1);
		}
              for(i=0;i<nbytes;i++)
		{
		if(sequence[i]=='\0')
			break;
		}
		sequence_len=i;

              if(fgets(quals, nbytes, stdin)==NULL)
                {
                return(1);
                }
              if(fgets(quals, nbytes, stdin)==NULL)
                {
                return(1);
                }
                for(i=0;i<nbytes;i++)
                {
                if(quals[i]=='\0')
                        break;
                }
                qual_len=i;


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

