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

int main()
{
    
  Cgicc formData;       
       
   cout << "Content-type:text/html\r\n\r\n";
   cout << "<html>\n";
   cout << "<head>\n";
   cout << "<title>File Upload in CGI</title>\n";
   cout << "</head>\n";
   cout << "<body>\n";
 
//------------------------------inputs-------------------------------------------------------------------------------- 
       string str, str2, filename_tmp, protein_name;
       ofstream myfile_tmp;

       str = "<set_up>\n\t<inputs>\n";
       const_file_iterator file;

       form_iterator fvalue1 = formData.getElement("protein_name");
       if( !fvalue1->isEmpty() && fvalue1 != (*formData).end()){ 
         str = str + "\t\t<protein_name>" + (**fvalue1).c_str() + "</protein_name>\n";
         protein_name = (**fvalue1).c_str();
       }
       else
          cout << "No text entered for value1" << endl;
// -----------------------------------------------------------------------------------------------------------------
file = formData.getFile("seq_file");
       if(file != formData.getFiles().end()) {
          str = str + "\t\t<seq_file>" + file->getFilename() + "</seq_file>\n";
          filename_tmp = "/tmp/" + file->getFilename();
          myfile_tmp.open(filename_tmp);
          if (myfile_tmp.is_open()){
                 file->writeToStream(myfile_tmp);
                 cout << "Successfully written " << filename_tmp << endl;
                 myfile_tmp.close();
          }
          else
             cout << "Could not open " << filename_tmp << "for writing.";
       }
       else
          cout << " No sequence file uploaded" << endl;
       
//-----------------------------------------------------------------------       
       file = formData.getFile("upl_file1");
       if(file != formData.getFiles().end()) {
          str = str + "\t\t<upl_file>" + file->getFilename();
          filename_tmp = "/tmp/" + file->getFilename();
          myfile_tmp.open(filename_tmp);
          if (myfile_tmp.is_open()){
                 file->writeToStream(myfile_tmp);
                 cout << "Successfully written " << filename_tmp << endl;
                 myfile_tmp.close();
          }
          else
             cout << "Could not open " << filename_tmp << "for writing.";
       }
       else
          cout << " No sequence file uploaded" << endl;
       
       file = formData.getFile("upl_file2");
       if(file != formData.getFiles().end()) {
          str = str + "," + file->getFilename();
          filename_tmp = "/tmp/" + file->getFilename();
          myfile_tmp.open(filename_tmp);
          if (myfile_tmp.is_open()){
                 file->writeToStream(myfile_tmp);
                 cout << "Successfully written " << filename_tmp << endl;
                 myfile_tmp.close();
          }
          else
             cout << "Could not open " << filename_tmp << "for writing.";
       }

       file = formData.getFile("upl_file3");
       if(file != formData.getFiles().end()) {
          str = str + "," + file->getFilename();
          filename_tmp = "/tmp/" + file->getFilename();
          myfile_tmp.open(filename_tmp);
          if (myfile_tmp.is_open()){
                 file->writeToStream(myfile_tmp);
                // cout << "Successfully written " << filename_tmp << endl;
                 myfile_tmp.close();
          }
          else
             cout << "Could not open " << filename_tmp << "for writing.";
       }

       file = formData.getFile("upl_file4");
       if(file != formData.getFiles().end()) {
          str = str + "," + file->getFilename();
          filename_tmp = "/tmp/" + file->getFilename();
          myfile_tmp.open(filename_tmp);
          if (myfile_tmp.is_open()){
                 file->writeToStream(myfile_tmp);
                // cout << "Successfully written " << filename_tmp << endl;
                 myfile_tmp.close();
          }
          else
             cout << "Could not open " << filename_tmp << "for writing.";
       }

       str = str + "</upl_file>\n";

//-----------------------------------------------------------------------

       file = formData.getFile("h_bond_file");
       if(file != formData.getFiles().end()) {
          str = str + "\t\t<hbond_file>" + file->getFilename() + "</hbond_file>\n";
          filename_tmp = "/tmp/" + file->getFilename();
          myfile_tmp.open(filename_tmp);
          if (myfile_tmp.is_open()){
                 file->writeToStream(myfile_tmp);
               //  cout << "Successfully written " << filename_tmp << endl;
                 myfile_tmp.close();
          }
          else
             cout << "Could not open " << filename_tmp << "for writing.";
       }
       else
          str = str + "\t\t<hbond_file></hbond_file>\n";       

       file = formData.getFile("aco_file");
       if(file != formData.getFiles().end()) {
          str = str + "\t\t<ang_file>" + file->getFilename() + "</ang_file>\n";
          filename_tmp = "/tmp/" + file->getFilename();
          myfile_tmp.open(filename_tmp);
          if (myfile_tmp.is_open()){
                 file->writeToStream(myfile_tmp);
                 cout << "Successfully written " << filename_tmp << endl;
                 myfile_tmp.close();
          }
          else
             cout << "Could not open " << filename_tmp << "for writing.";
       }
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
                         //myfile.open ("example.txt");
                          if (myfile.is_open()){
                                myfile << str;
                                myfile << str2;
                                cout << "Successfully written ~/tmp/example.xml";
                          }
                          else{
                                 cout << "error in opening the file"<<endl;
                                 cout << myfile.is_open()<<endl;
                          }
                           myfile.close();

       // Close the HTML document
   cout << "</body>\n";
   cout << "</html>\n";
    
    return 0;   // To avoid Apache errors.
}
