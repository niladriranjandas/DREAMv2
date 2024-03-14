#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include "XmlDomDocument.h"

char* sepStringsdelim(const char *str, const char *delim){
  std::string ans, ans_ ;
  char *pch, *ans_1, *tmp;

    tmp = (char *)str;
    
    pch = strtok(tmp,delim);
    while (pch != NULL){
      ans.append("'");
      ans.append(pch);
      ans.append("',");
      pch = strtok (NULL, delim);
    }
  ans_ = ans.substr(0,ans.length()-1);
  
  ans_1 = (char*) malloc(sizeof(ans_));
  strcpy(ans_1,ans_.c_str());
  return ans_1;
}

int main(int argc, char** argv)
{
    string value;
    XmlDomDocument* doc;
    //XmlDomDocument* doc = new XmlDomDocument("div_conquer_params.xml");
    if ( argc != 2 ){
         cout << "usage: " << argv[0] <<"<xml filename name>"<<endl;
         exit(-1);
    }
    else
       doc = new XmlDomDocument(argv[1]);

    if (doc) {
 
         // ---------- input file names and paths ------------------------//
             value = doc->getChildValue("inputs",0,"protein_name",0);
             printf("\n INPUTS.protein_name = %s;",sepStringsdelim(value.c_str(),","));
             value = doc->getChildValue("inputs",0,"seq_file",0);
             printf("\n INPUTS.seq_file = %s;",sepStringsdelim(value.c_str(),","));
             value = doc->getChildValue("inputs",0,"upl_file",0);
             printf("\n INPUTS.upl_file = {%s};",sepStringsdelim(value.c_str(),","));
             value = doc->getChildValue("inputs",0,"hbond_file",0);
            if(strlen(value.c_str()) == 0)
               printf("\n INPUTS.hbond_file = '';");         
            else
                printf("\n INPUTS.hbond_file = %s;",sepStringsdelim(value.c_str(),","));
             value = doc->getChildValue("inputs",0,"ang_file",0);
            if(strlen(value.c_str()) == 0)
               printf("\n INPUTS.ang_file = '';");         
            else
               printf("\n INPUTS.ang_file = %s;",sepStringsdelim(value.c_str(),","));
             value = doc->getChildValue("inputs",0,"protein_path",0);
             printf("\n INPUTS.protein_path = %s;",sepStringsdelim(value.c_str(),","));

         // ---------- Divide and conquer params -------------------------//
             value = doc->getChildValue("params",0,"hydrogen_omission",0);
             printf("\n PARAMS.hydrogen_omission = %s;", value.c_str());
             value = doc->getChildValue("params",0,"aug_bounds",0);
             if(strlen(value.c_str()) == 0)
                 printf("\n PARAMS.aug_bounds = 0;");         
             else
                 printf("\n PARAMS.aug_bounds = %s;", value.c_str());      
             value = doc->getChildValue("break_graph",0,"eta_lo",0);
             printf("\n PARAMS.eta_lo = %s;", value.c_str());      
             value = doc->getChildValue("break_graph",0,"eta_hi",0);
             printf("\n PARAMS.eta_hi = %s;", value.c_str());      
             value = doc->getChildValue("params",0,"include_neighbour",0);
             printf("\n PARAMS.include_neighbour = %s;", value.c_str());
             value = doc->getChildValue("multi_expand",0,"k_2",0);
             printf("\n PARAMS.multi_expand_k2 = %s;", value.c_str());      
             value = doc->getChildValue("multi_expand",0,"size_cutoff",0);
             printf("\n PARAMS.multi_expand_size_cutoff = %s;", value.c_str());      
             value = doc->getChildValue("multi_expand",0,"grp_min",0);
             printf("\n PARAMS.multi_expand_grp_min = %s;", value.c_str());      
             value = doc->getChildValue("multi_expand",0,"incr_min",0);
             printf("\n PARAMS.multi_incr_min = %s;", value.c_str());      
             value = doc->getChildValue("params",0,"grp_expand",0);
             printf("\n PARAMS.grp_expand = %s;", value.c_str());      

      delete doc;
    }

    exit(0);
 }
