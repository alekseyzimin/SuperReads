#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>

#include<charb.hpp>
#define DEBUG 0

typedef ExpandingBuffer<int> int_buf;
typedef ExpandingBuffer<int_buf, remaper<int_buf> > twoD_int_buf;
typedef ExpandingBuffer<charb, remaper<charb> > twoD_charb_buf;

static int int_compare(const void *p1, const void *p2){
  if(*(const int *) p1 < *(const int *)p2)
    return(-1);
  else if(*(const int *) p1 > *(const int *)p2)
    return(1);
  else
    return(0);
}


void reverse_sr(const char * original, char * reversed){
  //will write the reversed sr into *reversed. Assume reversed is at
  //least as large as original.
  charb original_local=original;//because strtok_r changes input
  char * token, *savespace;
  int offset=strlen(original);
  
  reversed[offset]='\0';

  token=strtok_r(original_local,"_",&savespace);
  
  if(token[strlen(token)-1]=='F')
    token[strlen(token)-1]='R';
  else
    token[strlen(token)-1]='F';

  offset-=strlen(token);
  for(size_t i=0;i<strlen(token);i++)
    reversed[offset+i]=token[i];

  while(1){
    token=strtok_r(NULL,"_",&savespace);
    if(token==NULL)
      break;
    if(token[strlen(token)-1]=='F')
      token[strlen(token)-1]='R';
    else if(token[strlen(token)-1]=='R')
      token[strlen(token)-1]='F';
    reversed[--offset]='_';
    offset-=strlen(token);
    for(size_t i=0;i<strlen(token);i++)
      reversed[offset+i]=token[i];
  }
}

int main(int argc,char *argv[]){

  int i,j,k,l=0,lastKUnitigIndex,irreducibleSuperReadIndex=0;
  int_buf kUnitigsInSuperRead(1000);
  time_t time_start=time(NULL);
  int_buf candidates(100);
  twoD_charb_buf irreducibleSuperReadNames(10000000);
  charb line(1000000);
  char *token, *saveptr;
  //parse arguments here
  twoD_int_buf superReadIndicesForKUnitig(atoi(argv[1]));

  while(fgets(line,stdin)){
    l++;	
    if(l%500000==0){
      fprintf(stderr,"Processed %d super reads, irreducible %d, processing %d super reads per second\n",l,irreducibleSuperReadIndex,(int)floor(500000/difftime(time(NULL),time_start)));
      time_start=time(NULL);
    }
    //first we parse the super read line, space separated, to get the name
    token=strtok_r(line," ",&saveptr);
    if(!token) {
      fprintf(stderr, "Invalid empty line number %d, expected a SuperRead name\n", l);
      continue;
    }
#if DEBUG
    printf("Found super read %s\n",token);
#endif
    //we need to save the super read name
    charb superReadName=token;
    charb superReadName_save=token;
    //here we need to examine first and last k-unitig in the super read name
    j=0;
    token=strtok_r(superReadName,"FR_",&saveptr);
    for(i=0; token != NULL;i++){
      if(i%2==0)
	kUnitigsInSuperRead[j++]=atoi(token);
      token=strtok_r(NULL,"FR_",&saveptr); 
    }
    lastKUnitigIndex=j-1;
#if DEBUG
    printf("first k_unitig: %d last k_unitig %d, super read %s candidates for first %d candidates for last %d\n",kUnitigsInSuperRead[0],kUnitigsInSuperRead[lastKUnitigIndex],(char*)superReadName_save,(int)superReadIndicesForKUnitig[kUnitigsInSuperRead[0]].size(),(int)superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]].size());
#endif



    //now we try to reduce
    k=0;
    if(superReadIndicesForKUnitig[kUnitigsInSuperRead[0]].size()>0&&superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]].size()>0){
    int max_first_index=0;

    if(superReadIndicesForKUnitig[kUnitigsInSuperRead[0]][0]>superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]][0])
	max_first_index=(int)superReadIndicesForKUnitig[kUnitigsInSuperRead[0]][0];
    else
        max_first_index=(int)superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]][0];


    for(i=0;i<(int)superReadIndicesForKUnitig[kUnitigsInSuperRead[0]].size();i++)
      if(superReadIndicesForKUnitig[kUnitigsInSuperRead[0]][i]>=max_first_index)
	      candidates[k++]=superReadIndicesForKUnitig[kUnitigsInSuperRead[0]][i];

    for(i=0;i<(int)superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]].size();i++)
     if(superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]][i]>=max_first_index)
	      candidates[k++]=superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]][i];

#if DEBUG
    printf("Found %d candidates\n",k);
#endif

#if DEBUG
    for(i=0;i<k;i++)
      printf("Before qsort:Candidate %d, super read %s\n",candidates[i], (char*)irreducibleSuperReadNames[candidates[i]]);
#endif

    qsort(candidates,k,sizeof(int),int_compare);

    charb superReadName_reverse(superReadName_save.size());

    reverse_sr(superReadName_save,superReadName_reverse);
    //now we go through the sorted candidates and figure out which one matches
#if DEBUG
    for(i=0;i<k;i++)
      printf("Candidate %d, super read %s\n",candidates[i], (char*)irreducibleSuperReadNames[candidates[i]]);
#endif

    //now we go through the sorted candidates and figure out which one matches, we only look at the candidate if it is encountered twice
    int last_candidate=-1;
    int candidate_count=1;
    for(i=0;i<k;i++){
      if(candidates[i]==last_candidate && candidate_count==1) {
#if DEBUG
	printf("Checking candidate %s\n",(char*)irreducibleSuperReadNames[candidates[i]]);
#endif
	if(strstr(irreducibleSuperReadNames[candidates[i]],superReadName_save)!=NULL)
	  break;
	if(strstr(irreducibleSuperReadNames[candidates[i]],superReadName_reverse)!=NULL)
	  break;
	candidate_count++;
      } else {
	candidate_count=1;
      }
      last_candidate=candidates[i];
    }

#if DEBUG
    if(i<k){
      printf("Reduced %s to %s\n",(char*)superReadName_save,(char*)irreducibleSuperReadNames[candidates[i]]); 	
      continue;
    }
#else
    if(i<k){
      printf("%s %s\n",(char*)superReadName_save,(char*)irreducibleSuperReadNames[candidates[i]]);
      continue;
    }
#endif
}
#if DEBUG
    printf("Irreducible %s, index %d\n",(char*)superReadName_save,irreducibleSuperReadIndex);
#endif

    //if we got here, then the super read is irreducible :(
    //here is what we do with an irreducible super read
    irreducibleSuperReadNames[irreducibleSuperReadIndex]=superReadName_save;
    for(i=0;i<=lastKUnitigIndex;i++)
      superReadIndicesForKUnitig[kUnitigsInSuperRead[i]].push_back(irreducibleSuperReadIndex);
    irreducibleSuperReadIndex++;
  }
  return(0);
}
