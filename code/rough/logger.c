#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>


static int const buffer_size_write = 50;//1000;

char* logger(const char* prev_string ,const char* log_string, int force_write_flag,const char* filename)
{
  FILE *fp;
  int len_now, len_prev;
  char *workstr, *time_now;
  time_t time_;
  struct tm * timeinfo;
      
  //--------check for NULL --------//
    if (prev_string == NULL || log_string == NULL || filename == NULL){
           perror("\n Arguments cannot be null");
           return "error";
    }
  //--------------------------------

  len_prev = strlen(prev_string);
  len_now  = strlen(log_string);

  if (len_prev == 0 || prev_string[0]=='\0'){
     time_ = time(NULL);
     time_now = ctime(&time_);
     workstr = (char *)calloc(strlen(time_now) + len_now + 2, sizeof(char));
     strcpy(workstr, time_now);
  }  
  else{
     workstr = (char *)calloc(len_prev+len_now + 2, sizeof(char));
     strcpy(workstr, prev_string);
  }

  strcat(workstr, log_string);
  strcat(workstr,"\n");

  if  ( (strlen(workstr) > buffer_size_write) || force_write_flag){
       fp = fopen(filename, "a");
       if (fp != NULL)
            fwrite(workstr, 1, strlen(workstr), fp);
       fclose(fp);
  }    
  else
      return workstr;
     
}

int main(int argc, char** argv)
{
    int i;
    char *abc = "";
//    abc = logger("abc","efg123", 0, "abc.txt");
    do{
       abc = logger(abc,"efg123", 0, "abc.txt");
       printf("%s",abc);
       printf("-------------------------------\n");
       i = i+1;
       }while(i<20);
return 0;
}

