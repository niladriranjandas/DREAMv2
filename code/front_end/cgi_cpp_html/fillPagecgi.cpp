#include <iostream>
#include <vector>
#include <string>

#include "cgicc/CgiDefs.h"
#include "cgicc/Cgicc.h"
#include "cgicc/HTTPHTMLHeader.h"
#include "cgicc/HTMLClasses.h"

#include <stdio.h>
#include <stdlib.h>

#include <fstream>

using namespace std;
using namespace cgicc;     // Or reference as cgicc::Cgicc formData; below in object instantiation.

int main(int argc, char **argv)
{
    try {
       Cgicc formData;
       string str, str2, protein_name;
       // Send HTTP header: Content-type: text/html
       cout << HTTPHTMLHeader() << endl;

       // Print: <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN" "http://www.w3.org/TR/REC-html40/strict.dtd">
       cout << HTMLDoctype(HTMLDoctype::eStrict) << endl;

       // Print: <html lang="en" dir="LTR">
       cout << html().set("lang", "EN").set("dir", "LTR") << endl;

       // Set up the HTML document
       cout << html() << head() << title("Cgicc example") << head() << endl;
       cout << body().set("bgcolor","#cccccc").set("text","#000000").set("link","#0000ff").set("vlink","#000080") << endl;

       cout << h1("This is a demonstration of the GNU CgiCC library") << endl;

 
//------------------------------inputs-------------------------------------------------------------------------------- 
       str = "<set_up>\n\t<inputs>\n";

       form_iterator fvalue1 = formData.getElement("protein_name");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()){ 
         str = str + "\t\t<protein_name>" + (**fvalue1).c_str() + "</protein_name>\n";
         protein_name = (**fvalue1).c_str();
       }
       else
          cout << "No text entered for value1" << endl;

       fvalue1 = formData.getElement("seq_file");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
          str = str + "\t\t<seq_file>" + (**fvalue1).c_str() + "</seq_file>\n";       
       else
          cout << "No text entered for value1" << endl;

       fvalue1 = formData.getElement("upl_file_1");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
          str = str + "\t\t<upl_file>" + (**fvalue1).c_str();
       else
          cout << "No text entered for value1" << endl;
       fvalue1 = formData.getElement("upl_file_2");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
          str = str + "," + (**fvalue1).c_str();
       fvalue1 = formData.getElement("upl_file_3");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
          str = str + "," + (**fvalue1).c_str();
       fvalue1 = formData.getElement("upl_file_4");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
          str = str + "," + (**fvalue1).c_str();
       str = str + "</upl_file>\n";

       fvalue1 = formData.getElement("H_bond_file");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
          str = str + "\t\t<hbond_file>" + (**fvalue1).c_str() + "</hbond_file>\n";
       else
          str = str + "\t\t<hbond_file></hbond_file>\n";       

       fvalue1 = formData.getElement("aco_file");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
          str = str + "\t\t<ang_file>" + (**fvalue1).c_str() + "</ang_file>\n";
       else          
          str = str + "\t\t<ang_file></ang_file>\n";       

       fvalue1 = formData.getElement("protein_path");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str = str + "\t\t<protein_path>" + (**fvalue1).c_str() + "</protein_path>\n";
       else
          cout << "No text entered for value1" << endl;

       str = str + "\t</inputs>\n";

//------------------------------params-------------------------------------------------------------------------------- 
      str2 = "\t<params>\n";
      fvalue1 = formData.getElement("h_ommission");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t<hydrogen_omission>" + (**fvalue1).c_str() + "</hydrogen_omission>\n";
       else
          cout << "No text entered for value1" << endl;


       fvalue1 = formData.getElement("aug_bounds");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t<aug_bounds>" + (**fvalue1).c_str() + "</aug_bounds>\n";
       else
          cout << "No text entered for value1" << endl;
       
       str2 = str2 + "\t\t<break_graph>\n";
       fvalue1 = formData.getElement("eta_lo");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t\t<eta_lo>" + (**fvalue1).c_str() + "</eta_lo>\n";
       else
          cout << "No text entered for value1" << endl;
       fvalue1 = formData.getElement("eta_hi");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t\t<eta_hi>" + (**fvalue1).c_str() + "</eta_hi>\n";
       else
          cout << "No text entered for value1" << endl;
       str2 = str2 + "\t\t</break_graph>\n";

       fvalue1 = formData.getElement("incl_neigh");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t<include_neighbour>" + (**fvalue1).c_str() + "</include_neighbour>\n";
       else
          cout << "No text entered for value1" << endl;

       str2 = str2 + "\t\t<multi_expand>\n";
       fvalue1 = formData.getElement("k_2");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t\t<k_2>" + (**fvalue1).c_str() + "</k_2>\n";
       else
          cout << "No text entered for value1" << endl;

       fvalue1 = formData.getElement("size_cutoff");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t\t<size_cutoff>" + (**fvalue1).c_str() + "</size_cutoff>\n";
       else
          cout << "No text entered for value1" << endl;

       fvalue1 = formData.getElement("grp_min");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t\t<grp_min>" + (**fvalue1).c_str() + "</grp_min>\n";
       else
          cout << "No text entered for value1" << endl;
       
       fvalue1 = formData.getElement("incr_min");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t\t<incr_min>" + (**fvalue1).c_str() + "</incr_min>\n";
       else
          cout << "No text entered for value1" << endl;
       str2 = str2 + "\t\t</multi_expand>\n";

       fvalue1 = formData.getElement("grp_expand");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()) 
           str2 = str2 + "\t\t<grp_expand>" + (**fvalue1).c_str() + "</grp_expand>\n";
       else
          cout << "No text entered for value1" << endl;

       str2 = str2 + "\t</params>\n";
       str2 = str2 + "</set_up>";
//------------------------------write to file------------------------------------------------------------------------- 
                          ofstream myfile;
                         //  myfile.open ("/home/niladri/Documents/xml_parser_try/cgi_try/example.txt");
                          string filename = "/tmp/" + protein_name + ".xml";
                          myfile.open (filename);
                         // myfile.open ("/tmp/example.xml");
                         
                          if (myfile.is_open()){
                                myfile << str;
                                myfile << str2;
                          }
                          else{
                                 cout << "error in opening the file"<<endl;
                                 cout << myfile.is_open()<<endl;
                          }
                           myfile.close();

       // Close the HTML document
       cout << body() << html();
    }
    catch(exception& e) {
       // handle any errors here.
       cout << "ERROR!!" << endl;
    }
    return 0;   // To avoid Apache errors.
}
