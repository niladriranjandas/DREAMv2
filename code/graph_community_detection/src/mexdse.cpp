#include <iostream>
#include "dse.h"
#include <math.h>
#include <parallel/algorithm>
#include <Eigen/SparseCore>
#include <Eigen/Core>

#include <mex.h>

/* Definitions to keep compatibility with earlier versions of ML */
// #ifndef MWSIZE_MAX
// typedef int mwSize;
// typedef int mwIndex;
// typedef int mwSignedIndex;
// 
// #if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
// /* Currently 2^48 based on hardware limitations */
// # define MWSIZE_MAX    281474976710655UL
// # define MWINDEX_MAX   281474976710655UL
// # define MWSINDEX_MAX  281474976710655L
// # define MWSINDEX_MIN -281474976710655L
// #else
// # define MWSIZE_MAX    2147483647UL
// # define MWINDEX_MAX   2147483647UL
// # define MWSINDEX_MAX  2147483647L
// # define MWSINDEX_MIN -2147483647L
// #endif
// #define MWSIZE_MIN    0UL
// #define MWINDEX_MIN   0UL
// #endif

#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;
#endif

#define DEBUG 0

mxArray * getMexArray(const std::vector<int>& v){
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}

// the supplied index of row and cols should be upper triangular

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     
 //declare variables
    mxArray *a_in_m, *b_in_m, *nodes_in_m, *c_out_m, *d_out_m, *param_density, *ret_cell, *cell_ele;
    const mwSize *dimsr , *dimsc;
    double *a, *b, *num_nodes, *c, *d, *density_param, *cell_ele_ptr;
    int dimr, dimc, numdimsr, numdimsc;
  //  int i,j;
    
  //associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    b_in_m = mxDuplicateArray(prhs[1]);
    nodes_in_m = mxDuplicateArray(prhs[2]);
    param_density = mxDuplicateArray(prhs[3]);
    
  //figure out dimensions
    dimsr = mxGetDimensions(prhs[0]);
    numdimsr = mxGetNumberOfDimensions(prhs[0]);

    dimsc = mxGetDimensions(prhs[1]);
    numdimsc = mxGetNumberOfDimensions(prhs[1]);

    dimr = (int)dimsr[0]; dimc = (int)dimsc[0];

//associate outputs
  //   c_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

//associate pointers
    a = mxGetPr(a_in_m);
    b = mxGetPr(b_in_m);
    num_nodes = mxGetPr(nodes_in_m);
    density_param = mxGetPr(param_density);

  //  c = mxGetPr(c_out_m);
    
//---------------check the inputs---------------------
   if (nrhs !=4){
      mexErrMsgIdAndTxt("Graph_community_master:mexdse:nrhs",
                      "Four inputs required.");   
   }
    
   if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])){
        mexErrMsgIdAndTxt("Graph_community_master:mexdse:prhs[0]",
                      "1st input argument should be double.");  
   }
   
   if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])){
        mexErrMsgIdAndTxt("Graph_community_master:mexdse:prhs[1]",
                      "2nd input argument should be double.");  
   }
    
   if( mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1])){
          mexErrMsgIdAndTxt("Graph_community_master:mexdse:prhs[1]",
                      "Number of elements in 1st and 2nd input argument should be equal."); 
   }
    
   if( !mxIsDouble(prhs[2]) || 
          mxIsComplex(prhs[2]) ||
          mxGetNumberOfElements(prhs[2]) != 1 ) {
          mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                       "3rd input argument, number of nodes must be a scalar.");
    }
    
    
   if( !mxIsDouble(prhs[3]) || 
          mxIsComplex(prhs[3]) ||
          mxGetNumberOfElements(prhs[3]) != 1 ||
          *density_param < 0  || *density_param >1){
           mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar",
                       "4th input argument, must be a scalar within the range 0 to 1."); 
   }
     
//-----------------------------------------------------  

// The actual logic begins
    std::cerr << "reading graph...";
//Eigen::SparseMatrix<double> g = readGraph(std::string("sample.graph"));
    Eigen::SparseMatrix<double> g = readMyFunc(a,dimr,b,dimc,*num_nodes);
    std::cerr << "done." << std::endl;

#if DEBUG
    std::cout << "Adjacency Matrix" << std::endl;
    std::cout << g << std::endl;
#endif

    std::vector<C_tuple*> C;
    std::cerr << "getting C vector...";
    getC(g, C);

    std::cerr << "sorting...";
    //__gnu_parallel::sort(C.begin(), C.end(), C_tuple_compare);   // this causes linking error
                                                                   // while compiling mex
    std::sort(C.begin(), C.end(), C_tuple_compare);
//     std::vector<C_tuple> C_test = modifyC(C,10-1);
//         std::cout<<std::endl;
//         for (uint64_t i = 0; i < C.size(); i++)
//         {
//            std::cout << (C[i]->i)+1 << ", " << (C[i]->j)+1 << ", " << C[i]->value<< std::endl;
//         }
//         std::cout<<"-------modified C---"<<std::endl;
//         for (uint64_t i = 0; i < C_test.size(); i++)
//         {
//            C[i] = &C_test[i];
//         }
//         std::cout<<"----new C-----"<<std::endl;
//         for (uint64_t i = 0; i < C.size(); i++)
//         {
//            std::cout << (C[i]->i)+1 << ", " << (C[i]->j)+1 << ", " << C[i]->value<< std::endl;
//         }
    
    std::cerr << "done." << std::endl;

#if DEBUG
    for (uint64_t i = 0; i < C.size(); i++)
    {
        std::cout << C[i]->i << ", " << C[i]->j << ", " << C[i]->value<< std::endl;
    }
#endif

    std::cerr << "creating T...";
    node* root = createT(g, C);
    std::cerr << "done." << std::endl;

#if DEBUG
    //postorder(printNode, root);
    levelorder(printNode, root);

    node* left = root->leftChild;
    while (left->leftChild) left = left->leftChild;

    node* right = root->rightChild;
    while (right->rightChild) right = right->rightChild;

    node* lca = LCA(left, right);
    std::cout << "lca for " << left->vertex << ", " << right->vertex << ": " << lca->vertex << std::endl;

    left = right->parent->leftChild;
    lca = LCA(left, right);
    std::cout << "lca for " << left->vertex << ", " << right->vertex << ": " << lca->vertex << std::endl;
#endif

    std::cerr << "counting vertices and edges...";
    countVerticesAndEdges(g, root);
    std::cerr << "done." << std::endl;
    std::cerr << "computing density...";
    computeDensity(root);
    std::cerr << "done." << std::endl;

#if DEBUG
    std::cout << "\nPrinting after the countVerticesAndEdges\n" << std::endl;
    postorder(printNode, root);
#endif

    std::cerr << "extracting subgraphs...";
    //extractSubgraphs(root, 0.75);    
    std::vector< std::vector<int> > dense_grp;
    extractSubgraphs(root,*density_param, &dense_grp);
    std::cout<<"---size dense_grp-----"<<dense_grp.size()<<std::endl;
    //alocate the cell to return
    ret_cell = mxCreateCellMatrix(dense_grp.size(),1);
                 for (std::vector< std::vector<int> >::size_type i = 0; i<dense_grp.size(); ++i)
                 {
                   std::vector<int> tmp  ;
                   size_t size_arr = dense_grp[i].size();
                  // cell_ele = mxCreateNumericArray(1,size_arr,mxDOUBLE_CLASS, mxREAL);
                  // cell_ele_ptr = mxGetPr(cell_ele);
                   for(std::vector<int>::size_type j=0;j<dense_grp[i].size();++j)
                   {
                     std::cout<<dense_grp[i][j]<<" ";
                   //  cell_ele_ptr[j] = dense_grp[i][j];          
                   }
                   cell_ele = getMexArray(dense_grp[i]);
                   std::cout<<std::endl;
                   mxSetCell(ret_cell,i,mxDuplicateArray(cell_ele));
                   mxDestroyArray(cell_ele);
                 }
    std::cerr << "done." << std::endl;
    plhs[0] = ret_cell;
   
//-------------------------------------------------------------------------

// // for testing purpose extractSubgraphsWithNode is called for vertex 10 hardcorded
//     std::cerr << "extracting subgraphs...";
//     //extractSubgraphs(root, 0.75);    
//     std::vector< std::vector<int> > dense_grp_node;
//     extractSubgraphsWithNode(root,*density_param, 10-1,&dense_grp_node);
//     std::cout<<"---size dense_grp-----"<<dense_grp_node.size()<<std::endl;
//     //alocate the cell to return
//     ret_cell = mxCreateCellMatrix(dense_grp_node.size(),1);
//                  for (std::vector< std::vector<int> >::size_type i = 0; i<dense_grp_node.size(); ++i)
//                  {
//                    std::vector<int> tmp  ;
//                    size_t size_arr = dense_grp_node[i].size();
//                   // cell_ele = mxCreateNumericArray(1,size_arr,mxDOUBLE_CLASS, mxREAL);
//                   // cell_ele_ptr = mxGetPr(cell_ele);
//                    for(std::vector<int>::size_type j=0;j<dense_grp_node[i].size();++j)
//                    {
//                      std::cout<<dense_grp_node[i][j]<<" ";
//                    //  cell_ele_ptr[j] = dense_grp[i][j];          
//                    }
//                    cell_ele = getMexArray(dense_grp_node[i]);
//                    std::cout<<std::endl;
//                    mxSetCell(ret_cell,i,mxDuplicateArray(cell_ele));
//                    mxDestroyArray(cell_ele);
//                  }
//     std::cerr << "done." << std::endl;
//     plhs[1] = ret_cell;
//-------------------------------------------------------------------------    

}

