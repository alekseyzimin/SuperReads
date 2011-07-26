#include<stdio.h>
#include<stdlib.h>
int main(int argc, char *argv[])
{
  FILE *FP;
  int i,max_read_num,n;
  int nbytes=10000000;
  char *my_string;
  char *quals;
  char *read_status;
  char *seq;
  
  my_string = malloc (nbytes);
  seq =  malloc (nbytes);
  quals = malloc (nbytes);
 
  if(argc<6)
	{
	printf("need five arguments: lib_id mean stdev number_of_reads fasta_file\n");
	return(1);
	}

max_read_num=atoi(argv[4]);
read_status=malloc(max_read_num+2);
for(i=0;i<=max_read_num;i++)
{
read_status[i]='N';
}

FP=fopen(argv[5],"r");
printf("{VER\nver:2\n}\n{LIB\nact:A\nacc:%s\nori:I\nmea:%s\nstd:%s\nsrc:\n.\nnft:1\nfea:\ndoNotOverlapTrim=1\n.\n}\n",argv[1],argv[2],argv[3]);
while(fgets (my_string, nbytes, FP)!=NULL)
    {
          if(my_string[0]=='>')
            {                       
             for(i=0;i<nbytes-1;i++)
             {
                if( my_string[i]==' ' || my_string[i]=='\n')
                {
                my_string[i]='\0';
                break;
                }
                else if( my_string[i]=='\0')
                {
                break;
                }
             }
             n=atoi(my_string+3);
             
	     fgets (seq, nbytes, FP);
             for(i=0;i<nbytes;i++)
             {
                if(seq[i]=='\n')
                {
                break;
                }
             }
             if(i>=64)
                {
                read_status[n]='Y';
                }
	/*printf("DEBUG %s %s %d %d\n", my_string,seq,n,i); */
        }
}
fclose(FP);
FP=fopen(argv[5],"r");

while(fgets (my_string, nbytes, FP)!=NULL)
    {
	  if(my_string[0]=='>')
	    {
	     for(i=0;i<nbytes-1;i++)
             {
                if( my_string[i]==' ' || my_string[i]=='\n')
		{
		my_string[i]='\0';
		break;
		}
		else if( my_string[i]=='\0')
		{
		break;
		}
	     }
             n=atoi(my_string+3);
	     if( read_status[n]=='N')
		continue;
 
	     fgets (seq, nbytes, FP);
             for(i=0;i<nbytes;i++)
             {
                if(seq[i]!='\n')
                {
                quals[i]='E';
                }
		else
		{
		quals[i]='\n';
		quals[i+1]='\0';
		break;	
		}
             }
	     if(n%2==0)
		{
		if(read_status[n]=='Y' && read_status[n+1]=='Y')
	              printf("{FRG\nact:A\nacc:%s\nrnd:1\nsta:G\nlib:%s\npla:0\nloc:0\nsrc:\n.\nseq:\n%s.\nqlt:\n%s.\nhps:\n.\nclr:0,%d\n}\n",my_string+1,argv[1],seq,quals,i);
		}
		else
                {
                if(read_status[n-1]=='Y' && read_status[n]=='Y')
                      printf("{FRG\nact:A\nacc:%s\nrnd:1\nsta:G\nlib:%s\npla:0\nloc:0\nsrc:\n.\nseq:\n%s.\nqlt:\n%s.\nhps:\n.\nclr:0,%d\n}\n",my_string+1,argv[1],seq,quals,i);
		}
	}
}
for(i=0;i<=max_read_num;i+=2)
{
if(read_status[i]=='Y' && read_status[i+1]=='Y')
{
printf("{LKG\nact:A\nfrg:%s%d\nfrg:%s%d\n}\n",argv[1],i,argv[1],i+1);
}
}
fclose(FP);
return(0);
}

