#include <map>
#include <vector>
#include <algorithm> 
#include <iostream>

#include "extractconstraints.h"

std::multimap<int, int> makeMultimap(std::vector<int> row_i){
 using namespace std;
   int i=0;    
   multimap<int, int> map_rowi;
   vector<int>::iterator it = row_i.begin();   

   for(; it != row_i.end() ; i++,++it)
       map_rowi.insert(pair<int, int>(*it,i));

  return(map_rowi);
}

std::map<int,int> makeMap(std::vector<int> row_i){
  using namespace std;
     int i=0;
     map<int,int> map_rowi;
     vector<int>::iterator it = row_i.begin();

     for(;it != row_i.end(); i++,++it)
         map_rowi.insert(pair<int, int>(*it,i));
   
   return(map_rowi);
}

std::vector<int> makeVector(std::multimap<int, int>::const_iterator values, int count){
  using namespace std;
    vector<int> retvect;
    for(int i=0;i<count;i++,++values)
      retvect.push_back(values->second);

  return(retvect);
}

std::vector<std::pair<int,int>> reIndxByAtomGrp(std::vector<int> row1, std::vector<int> row2, std::vector<int> atom_grp){
   using namespace std;
        map<int,int> map_atomgrp = makeMap(atom_grp);
        pair<int,int> tmp;
        vector<pair<int,int>> ret;
        map<int,int>::const_iterator value1, value2;
        vector<int>::iterator it1,it2;

        for(it1 = row1.begin(), it2 = row2.begin(); (it1 != row1.end()) && (it2 != row2.end()); ++it1,++it2){
            value1 = map_atomgrp.find(*it1);
            value2 = map_atomgrp.find(*it2);
            ret.push_back(pair<int,int>(value1->second, value2->second));      
        }
   return ret;
}

std::vector<std::pair<std::pair<int,int>,double>> getExtractedCons(std::vector<int> dist_i, std::vector<int> dist_j, std::vector<double> dist_ij, std::vector<int> grp_atoms){
  using namespace std;
        vector<int> member_indx_1, member_indx_2, pos, tmp;
       
 /* ------ initialize empty vector ----------- */
    member_indx_1 = vector<int>();
    member_indx_2 = vector<int>();

  // use multimap to hasp map with chain dist_i and dist_j as (dist_i, indexes (in case of repetation in dist_i) same for j
    multimap<int, int> row_1, row_2;
    multimap<int, int>::const_iterator values;
    row_1 = makeMultimap(dist_i);
    row_2 = makeMultimap(dist_j);

   // querry the multimap each for dist_i and dist_j to get indices for each element of grp_atoms
    for(vector<int>::iterator it = grp_atoms.begin(); it != grp_atoms.end(); ++it){
        values = row_1.find(*it);
        tmp = makeVector(values,row_1.count(*it));
        member_indx_1.insert(member_indx_1.end(), tmp.begin(), tmp.end());                
      
        values = row_2.find(*it);
        tmp = makeVector(values,row_2.count(*it));
        member_indx_2.insert(member_indx_2.end(), tmp.begin(), tmp.end());   
    }

   // sort the querry results from prov step and do set intersection after sorting
    sort(member_indx_1.begin(), member_indx_1.end());
    sort(member_indx_2.begin(), member_indx_2.end());
    set_intersection(member_indx_1.begin(), member_indx_1.end(), member_indx_2.begin(), member_indx_2.end(), back_inserter(pos));   

    vector<int> dist_i_pos, dist_j_pos;
    for(vector<int>::iterator it=pos.begin(); it != pos.end(); ++it){
           dist_i_pos.push_back(dist_i.at(*it));
           dist_j_pos.push_back(dist_j.at(*it));
    }
    
   // reindex the entry as per grp_atoms    
    vector<pair<int,int>> extract_cons = reIndxByAtomGrp(dist_i_pos,dist_j_pos,grp_atoms);

   // store the final result in an vector
    vector<pair<int,int>>::iterator it1=extract_cons.begin();
    vector<int>::iterator it2 = pos.begin();
    vector<pair<pair<int,int>,double>> finalresult;
   
    for(; it1 != extract_cons.end(); ++it1,++it2){
         pair<int,int> dist_i_dist_j = make_pair( (it1->first)+1, (it1->second)+1);                    //done +1 because lower index in matlab is 1 instead of 0
         pair<pair<int,int>,double> extracted_dist_cons = make_pair(dist_i_dist_j,dist_ij.at(*it2));
         finalresult.push_back(extracted_dist_cons);
    }

  return finalresult;
}

std::vector<std::pair<std::pair<int,int>,double>> getExtractedConsAnchorAnchor(std::vector<int> dist_i, std::vector<int> dist_j, std::vector<double> dist_ij, std::vector<int> anchors){
 using namespace std;
        vector<int> member_indx_1, member_indx_2, pos_a_v, pos_v_a, tmp;
       
 /* ------ initialize empty vector ----------- */
    member_indx_1 = vector<int>();
    member_indx_2 = vector<int>();

  // use multimap to hasp map with chain dist_i and dist_j as (dist_i, indexes (in case of repetation in dist_i) same for j
    multimap<int, int> row_1, row_2;
    multimap<int, int>::const_iterator values;
    row_1 = makeMultimap(dist_i);
    row_2 = makeMultimap(dist_j);

   // querry the multimap each for dist_i and dist_j to get indices for each element of grp_atoms
    for(vector<int>::iterator it = anchors.begin(); it != anchors.end(); ++it){
        values = row_1.find(*it);
        tmp = makeVector(values,row_1.count(*it));
        member_indx_1.insert(member_indx_1.end(), tmp.begin(), tmp.end());                
      
        values = row_2.find(*it);
        tmp = makeVector(values,row_2.count(*it));
        member_indx_2.insert(member_indx_2.end(), tmp.begin(), tmp.end());   
    }

   // sort the querry results from prov step 
    sort(member_indx_1.begin(), member_indx_1.end());
    sort(member_indx_2.begin(), member_indx_2.end());
   // get the members with (anchors_i, vars_j , d_ij)
    set_difference(member_indx_1.begin(), member_indx_1.end(), member_indx_2.begin(), member_indx_2.end(), back_inserter(pos_a_v));
    set_difference(member_indx_2.begin(), member_indx_2.end(), member_indx_1.begin(), member_indx_1.end(), back_inserter(pos_v_a));

  //
    vector<pair<pair<int,int>,double>> finalresult;
    vector<int>::iterator it_a_v = pos_a_v.begin();
    vector<int>::iterator it_v_a = pos_v_a.begin();

    for(;it_a_v != pos_a_v.end(); ++it_a_v){
       pair<int,int> a_v = make_pair(dist_j.at(*it_a_v), dist_i.at(*it_a_v));
       pair<pair<int,int>,double> a_v_dij = make_pair(a_v, dist_ij.at(*it_a_v));
       finalresult.push_back(a_v_dij);
    }

    for(;it_v_a != pos_v_a.end(); ++it_v_a){
        pair<int,int> v_a = make_pair(dist_i.at(*it_v_a), dist_j.at(*it_v_a));
        pair<pair<int,int>,double> v_a_dij = make_pair(v_a, dist_ij.at(*it_v_a));
        finalresult.push_back(v_a_dij);
    }
  
  return finalresult;  
}

std::vector<std::pair<std::pair<int,int>,double>> getExtractedConsVarsVars(std::vector<int> dist_i, std::vector<int> dist_j, std::vector<double> dist_ij, std::vector<int> vars){
 using namespace std;
        vector<int> member_indx_1, member_indx_2, pos, tmp;
       
 /* ------ initialize empty vector ----------- */
    member_indx_1 = vector<int>();
    member_indx_2 = vector<int>();

  // use multimap to hasp map with chain dist_i and dist_j as (dist_i, indexes (in case of repetation in dist_i) same for j
    multimap<int, int> row_1, row_2;
    multimap<int, int>::const_iterator values;
    row_1 = makeMultimap(dist_i);
    row_2 = makeMultimap(dist_j);

   // querry the multimap each for dist_i and dist_j to get indices for each element of grp_atoms
    for(vector<int>::iterator it = vars.begin(); it != vars.end(); ++it){
        values = row_1.find(*it);
        tmp = makeVector(values,row_1.count(*it));
        member_indx_1.insert(member_indx_1.end(), tmp.begin(), tmp.end());                
      
        values = row_2.find(*it);
        tmp = makeVector(values,row_2.count(*it));
        member_indx_2.insert(member_indx_2.end(), tmp.begin(), tmp.end());   
    }

   // sort the querry results from prov step 
    sort(member_indx_1.begin(), member_indx_1.end());
    sort(member_indx_2.begin(), member_indx_2.end());
    set_intersection(member_indx_1.begin(), member_indx_1.end(), member_indx_2.begin(), member_indx_2.end(), back_inserter(pos));   

   //
    vector<pair<pair<int,int>,double>> finalresult;
    vector<int>::iterator it_v_v = pos.begin();

    for(;it_v_v != pos.end(); ++it_v_v){
       pair<int,int> v_v = make_pair(dist_i.at(*it_v_v), dist_j.at(*it_v_v));
       pair<pair<int,int>,double> v_v_dij = make_pair(v_v, dist_ij.at(*it_v_v));
       finalresult.push_back(v_v_dij);
    }

  return finalresult;
}
