//g++ -c makeConfPlotOfEnsemble.cpp -I/usr/include/python2.7 -L/usr/lib -lm -lpython2.7 -std=c++11
//g++ makeConfPlotOfEnsemble.o -o makeConfPlotOfEnsemble -I/usr/include/python2.7 -L/usr/lib -lm -lpython2.7 -std=c++11
//./makeConfPlotOfEnsemble 1a7m test1.tsv test2.tsv test3.tsv test4.tsv test5.tsv

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <Python.h>

// trim from both ends (in place)
std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

std::vector<double> getPhiOrPsi(std::string file_name,int ind){
	// ind=1: return phi	ind=2: return psi
	std::string each_line;
	std::ifstream inFile;

	inFile.open(file_name);

	char *cr, cr2[50];
	int i;

	std::vector<double> torsion;
	if(inFile.is_open()){
		while(std::getline(inFile, each_line)){
			cr=strtok(strcpy(cr2,each_line.c_str())," \t");
			i=0;
			for(i=0;cr!=NULL;cr=std::strtok(NULL," \t"),i++){
				if (ind==1 && i==1)
					torsion.push_back(stof(trim(cr)));
				if (ind==2 && i==2)
				       	torsion.push_back(stof(trim(cr)));
			}
		}
	} 
	else {
		std::cout<<"\nFile open error " << file_name << "\n";	
		torsion.push_back(1000);
	}
 return torsion;
}

int runPythonForConfplot(std::string pythoncode, std::string pythonfunc, std::vector<std::vector<double>> phi, std::vector<std::vector<double>> psi, std::string protein_name, std::string file_2ary){

    std::vector<double> phi_for_i, psi_for_i;
    
	PyObject *pName, *pModule, *pDict, *pFunc;
	PyObject *pArgTuple, *pValue, *pXVec, *pYVec;
    int i;
    
	Py_Initialize();
	
	// Set the path to include the current directory in case the module is located there. Found from
	// http://stackoverflow.com/questions/7624529/python-c-api-doesnt-load-module
	// and http://stackoverflow.com/questions/7283964/embedding-python-into-c-importing-modules
	//
	PyObject *sys = PyImport_ImportModule("sys");
	PyObject *path = PyObject_GetAttrString(sys, "path");
	PyList_Append(path, PyString_FromString("."));
	
//	pName = PyString_FromString(argv[1]);   //Get the name of the module
	pName = PyString_FromString(pythoncode.c_str());
	pModule = PyImport_Import(pName);     //Get the module
	
	Py_DECREF(pName);
	
	if (pModule != NULL) {
		pFunc = PyObject_GetAttrString(pModule, pythonfunc.c_str());  //Get the function by its name
		/* pFunc is a new reference */
		
		if (pFunc && PyCallable_Check(pFunc)) {

			//Set up a tuple that will contain the function arguments. In this case, the
			//function requires two tuples, so we set up a tuple of size 2.
            if(!file_2ary.empty())
                pArgTuple = PyTuple_New(phi.size()+psi.size()+4);
            else
                pArgTuple = PyTuple_New(phi.size()+psi.size()+3);
            
            int count=0;
            PyTuple_SetItem(pArgTuple, count, PyString_FromString(protein_name.c_str()));
            count+=1;
            PyTuple_SetItem(pArgTuple, count, PyString_FromString("phi"));
            count+=1;
            for(std::vector<std::vector<double>>::iterator it=phi.begin(); it != phi.end(); ++it,++count){
                phi_for_i = *it;
                // transfer each cpp vector to python tuple
                pXVec = PyTuple_New(phi_for_i.size());
                for (i = 0; i < phi_for_i.size(); ++i) {
                    pValue = PyFloat_FromDouble(phi_for_i[i]);
                    if (!pValue) {
                        Py_DECREF(pXVec);
                        Py_DECREF(pModule);
                        fprintf(stderr, "Cannot convert array value\n");
                        return 1;
                    }
                    PyTuple_SetItem(pXVec, i, pValue);
                }
                
                //Set the argument tuple to contain the tuple
                PyTuple_SetItem(pArgTuple, count, pXVec);
            }
            
            PyTuple_SetItem(pArgTuple, count, PyString_FromString("psi"));
            count+=1;
            for(std::vector<std::vector<double>>::iterator it=psi.begin(); it != psi.end(); ++it,++count){
                psi_for_i = *it;
                // transfer each cpp vector to python tuple
                pXVec = PyTuple_New(psi_for_i.size());
                for (i = 0; i < psi_for_i.size(); ++i) {
                    pValue = PyFloat_FromDouble(psi_for_i[i]);
                    if (!pValue) {
                        Py_DECREF(pXVec);
                        Py_DECREF(pModule);
                        fprintf(stderr, "Cannot convert array value\n");
                        return 1;
                    }
                    PyTuple_SetItem(pXVec, i, pValue);
                }
                
                //Set the argument tuple to contain the tuple
                PyTuple_SetItem(pArgTuple, count, pXVec);
            }
	    
	        if(!file_2ary.empty()){
                PyTuple_SetItem(pArgTuple, count, PyString_FromString(file_2ary.c_str()));
                count++;
            }

            //Call the python function
			pValue = PyObject_CallObject(pFunc, pArgTuple);
			Py_DECREF(pArgTuple);
			Py_DECREF(pXVec);

			if (pValue != NULL) {
				printf("Result of call: %ld\n", PyInt_AsLong(pValue));
				Py_DECREF(pValue);
			}
			//Some error catching
			else {
				Py_DECREF(pFunc);
				Py_DECREF(pModule);
				PyErr_Print();
				fprintf(stderr,"Call failed\n");
				return 1;
			}
			
		}
		else {
			if (PyErr_Occurred())
				PyErr_Print();
			fprintf(stderr, "Cannot find function \"%s\"\n", pythonfunc);
		}
		Py_XDECREF(pFunc);
		Py_DECREF(pModule);
	}
	else {
		PyErr_Print();
		fprintf(stderr, "Failed to load \"%s\"\n", pythoncode);
		return 1;
	}
	Py_Finalize();
	return 0;	

}


int main(int argc, char *argv[]){

	std::vector<double> phi;
	std::vector<double> psi;
    
    std::vector<std::vector<double>> ensemble_phi, ensemble_psi;
    int flag=0, maxitter=argc;

	if ( argc < 3 ){
		std::cout<<"Usage: "<<argv[0]<<"protein_name file1 file2 ...\n ";
		exit (EXIT_FAILURE);
	}

	std::size_t found = std::string(argv[argc-1]).find("-2ary");    
    if(found!=std::string::npos){
        flag=1;
        maxitter=argc-1;
    }
        
	for(int i=2;i<maxitter;i++){
		phi = getPhiOrPsi(argv[i],1);
		if(phi.size() == 1)
			std::exit(-1);
		psi = getPhiOrPsi(argv[i],2);
		if(psi.size() == 1)
			std::exit(-1);
		ensemble_phi.push_back(phi);
		ensemble_psi.push_back(psi);
	}
       
    if(flag)        
        runPythonForConfplot("pythonConfPlot","plotStdVectors",ensemble_phi,ensemble_psi,argv[1],std::string(argv[argc-1]).substr(6));
    else
        runPythonForConfplot("pythonConfPlot","plotStdVectors",ensemble_phi,ensemble_psi,argv[1],"");        
    
    return 0;
}


