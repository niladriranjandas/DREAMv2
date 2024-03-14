/* header file for extract cons. can also be used to get index of duplicate list */

#ifndef EXTRACTCONSTRAINTS_H
#define EXTRACTCONSTRAINTS_H

#include<map>
#include<vector>

std::multimap<int, int> makeMultimap(std::vector<int> row_i); //function to make multimap with <row_i,index in row_i>
std::map<int,int> makeMap(std::vector<int> row_i);            //function to make map with <row_i,index in row_i>

std::vector<int> makeVector(std::multimap<int, int>::const_iterator values, int count); //function to make vector out of the results returned by multimap search

std::vector<std::pair<int,int>> reIndxByAtomGrp(std::vector<int> row1, std::vector<int> row2, std::vector<int> atom_grp);  // reindex row1 and row2 according to the list atom_grp

std::vector<std::pair<std::pair<int,int>,double>> getExtractedCons(std::vector<int> dist_i, std::vector<int> dist_j, std::vector<double> dist_ij, std::vector<int> grp_atoms); // extracted cons reindexed as per atom_grp

#endif
