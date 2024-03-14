#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#include <mex.h>

static int const buffer_size_write = 1000;

char* logger(const char* prev_string ,const char* log_string, int force_write_flag,const char* filename)
{
  FILE *fp;
  int len_now, len_prev;
  char *workstr, *time_now;
  time_t time_;
  struct tm * timeinfo;
      
  /*--------check for NULL --------*/
    if (prev_string == NULL ||  filename == NULL){
           mexErrMsgTxt("\n Arguments cannot be null");
           return NULL;
    }
  /*--------------------------------*/

  len_prev = strlen(prev_string);
  len_now  = strlen(log_string);

  if (len_prev == 0 || prev_string[0]=='\0'){
     time_ = time(NULL);
     time_now = ctime(&time_);
     workstr = (char *)mxCalloc(strlen(time_now) + len_now + 2, sizeof(char));
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
       return "";
  }    
  else
      return workstr;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    char *prev_log, *new_str, *filename, *output_buff;
    int force_write_flag, *written_flag;
    size_t prev_log_len, new_str_len, filename_len;

 /*------------------------ check the i/p --------------------------------*/
   if(nrhs!=4) 
      mexErrMsgTxt( "usage: [appended_string, written_flag] = mexDolog <prev_string> <new_string> <force_write_flag(0/1)> <filename>.");
    else if(nlhs > 2) 
      mexErrMsgTxt( "error: Too many output arguments.");

 /*----------------------- get the i/p strings ---------------------------*/
    prev_log = mxArrayToString(prhs[0]);
    new_str  = mxArrayToString(prhs[1]);
    force_write_flag = (int)mxGetScalar(prhs[2]);
    filename = mxArrayToString(prhs[3]);

 /*----------------------- associate output ------------------------------*/
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    written_flag = (int*)plhs[1];

 /*-----------------------------call logger ------------------------------*/
    output_buff = logger(prev_log, new_str, force_write_flag, filename);
    if (output_buff == NULL)
         *written_flag = 1;
    plhs[0] = mxCreateString(output_buff);
 /* ------------------------free the memory ------------------------------*/
   mxFree(prev_log);
   mxFree(new_str);
   mxFree(filename);

 return;
}


