#include<stdio.h>
#include<stdlib.h>
#include <src/charb.hpp>
int main(int argc, char *argv[])
{
  ssize_t i,sequence_len;
  charb my_string(100),sequence(100);
  size_t index_array[10000];
  size_t max_val=0,max_ind=0;
  double num_cg=0;
  char last_letter;
  
  if(argc==1)
    {
      printf("Need one agrument -- minimum peak height to trim!\n");
      return(1);
    }

  while(1)
    {
      if(!fgets(my_string,stdin)){ 
	break;
	}
      else
	{
	  if(my_string[0]=='>')
	    {
	      max_val=0;
	      max_ind=0;
	      num_cg=0;
              if(!fgets(sequence, stdin))
                {
                return(1);
                }
              sequence_len=strlen(sequence);

	      last_letter=sequence[sequence_len-1];
   	      for(i=sequence_len-2;i>=0;i--)
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
			index_array[sequence_len-2-i]=num_cg;
		    }
	      for(i=0;i<sequence_len-1;i++)
		{
		if(index_array[i]>max_val)
			{
			max_val=index_array[i];
			max_ind=i+1;
			}
		}
		if(max_val<(size_t)atoi(argv[1]))
			max_ind=0;
		sequence[sequence_len-max_ind]='\0';
		printf("%s%s\n",(char*)my_string,(char *)sequence);
	}
	}
}
return(0);
}
