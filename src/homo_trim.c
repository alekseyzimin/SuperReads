#include<stdio.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
  int i;
  size_t nbytes = 1000000;
  ssize_t sequence_len;
  char *my_string;
  char *sequence;
  int max_val=0,max_ind=0;
  double num_cg=0;
  char last_letter;
  int *index_array;
  
  my_string = malloc (nbytes);
  sequence = malloc(nbytes);
  index_array = malloc(sizeof(int)*1000000);
  
  if(argc==1)
    {
      printf("Need one agrument -- minimum peak height to trim!\n");
      return(1);
    }

  while(1)
    {
      if(fgets(my_string, nbytes, stdin)==NULL)
	break;
      else
	{
	  if(my_string[0]=='>')
	    {
	      max_val=0;
	      max_ind=0;
	      num_cg=0;
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

	      last_letter=(char)sequence[sequence_len-2];
   	      for(i=sequence_len-3;i>=0;i--)
		    {
		    if(sequence[i] != last_letter)
			{
			num_cg--;
			last_letter=sequence[i];
			}
		    else
			{
			num_cg++;
			}
			index_array[sequence_len-3-i]=num_cg;
		    }
	      for(i=0;i<sequence_len-2;i++)
		{
		if(index_array[i]>max_val)
			{
			max_val=index_array[i];
			max_ind=i+1;
			}
		}
		if(max_val<atoi(argv[1]))
			max_ind=0;
		sequence[sequence_len-max_ind-1]='\0';
		printf("%s%s\n",my_string,sequence);
	}
	}
}
return(0);
}
