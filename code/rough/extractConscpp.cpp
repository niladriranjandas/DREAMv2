#include <map>
#include <vector>
#include <algorithm> 
#include <iostream>

using namespace std;

multimap<int, int> makeMultimap(vector<int> row_i){

   int i=0;    
   multimap<int, int> map_rowi;
   vector<int>::iterator it = row_i.begin();   

   for(; it != row_i.end() ; i++,++it)
       map_rowi.insert(pair<int, int>(*it,i));

  return(map_rowi);
}

map<int,int> makeMap(vector<int> row_i){
     int i=0;
     map<int,int> map_rowi;
     vector<int>::iterator it = row_i.begin();

     for(;it != row_i.end(); i++,++it)
         map_rowi.insert(pair<int, int>(*it,i));
   
   return(map_rowi);
}

vector<int> makeVector(multimap<int, int>::const_iterator values, int count){
    vector<int> retvect;
    for(int i=0;i<count;i++,++values)
      retvect.push_back(values->second);

  return(retvect);
}

vector<pair<int,int>> reIndxByAtomGrp(vector<int> row1, vector<int> row2, vector<int> atom_grp){
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


int main(){
   vector<int> dist_i, dist_j, dist_ij;
   vector<int> grp_atoms, member_indx_1, member_indx_2, pos, tmp;

    dist_i = {1,1,2,2,3,3,4,4,4,5,6,7,10,10};
    dist_j = {2,3,4,5,5,6,6,7,8,1,2,3,5,8};
    dist_ij = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    grp_atoms = {1,3,5,10};

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
           printf("(i,j): (%d,%d)\n",dist_i.at(*it), dist_j.at(*it));
    }
    
   // reindex the entry as per grp_atoms
    cout<<endl<<"reindexed"<<endl;
    vector<pair<int,int>> extract_cons = reIndxByAtomGrp(dist_i_pos,dist_j_pos,grp_atoms);

   // store the final result in an vector
    vector<pair<int,int>>::iterator it1=extract_cons.begin();
    vector<int>::iterator it2 = pos.begin();
    vector<pair<pair<int,int>,float>> finalresult;
   
    for(; it1 != extract_cons.end(); ++it1,++it2){
         pair<int,int> dist_i_dist_j = make_pair( (it1->first)+1, (it1->second)+1);
         pair<pair<int,int>,float> extracted_dist_cons = make_pair(dist_i_dist_j,dist_ij.at(*it2));
         finalresult.push_back(extracted_dist_cons);
    }

    for(vector<pair<pair<int,int>,float>>::iterator it=finalresult.begin(); it!=finalresult.end(); ++it){
           pair<int,int> a = it->first;
           float b = it->second;

           printf("\n(i,j,d_ij): (%d,%d,%f) ",a.first,a.second, b);
    }
    cout<<endl;
         
}
