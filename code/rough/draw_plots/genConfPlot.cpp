/*
 *The program expects a txt file having the following fields
 *<4-letter-pdb-code>:<chain name>:<residue 3 letter><residue_no>  <phi>  <psi>
 *e.g.
 *1HMP:chain0:ASP1	9999.000000	164.009064
 * compile with: g++ genConfPlot.cpp -std=c++11 -I/usr/include/python2.7 -I/home/niladri/Downloads/softwares/matplotlib-cpp-master -lpython2.7
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>

#include "matplotlibcpp.h"

namespace plt=matplotlibcpp;

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

void plotTorsions(std::vector<double> vals){
	std::vector<double> x(vals.size());
	std::iota(std::begin(x),std::end(x),0);
	plt::ylim(-180,180);
	plt::plot(x,vals,"r-o");
	plt::show();
}

void plotTorsionsEnsemble(std::vector<std::vector<double>> ensemble_torsions, std::string plotname){
	std::vector<double> torsion_for_i;
	plt::ylim(-180,180);
	int i=1;
	std::string name;

	plt::title("Confplot for torsional angles in ensemble");
	plt::xlabel("residues");
	plt::ylabel(plotname);
	for(std::vector<std::vector<double>>::iterator it=ensemble_torsions.begin(); it != ensemble_torsions.end(); ++it,++i){
		torsion_for_i = *it;
		std::vector<double> x(torsion_for_i.size());
		std::iota(std::begin(x),std::end(x),0);
		name = "ensemble: " + std::to_string(i);
		plt::named_plot(name, x, torsion_for_i, "-o");
	}
	plt::legend();
	plt::show();
}

void plotPhiPsiSubplotEnsemble(std::vector<std::vector<double>> phi_ensemble, std::vector<std::vector<double>> psi_ensemble, std::string protein_name){
	std::vector<double> phi, psi;
	int i=1;
	std::string name;

//	plt::title("Confplot for torsional angle in ensemble");
	
	plt::subplot(2,1,1);
	plt::title("Confplot for torsional angle in ensemble");
	plt::ylim(-180,180);
	plt::xlabel("residues");
	plt::ylabel("phi");
	        for(std::vector<std::vector<double>>::iterator it=phi_ensemble.begin(); it != phi_ensemble.end(); ++it,++i){
                phi = *it;
                std::vector<double> x(phi.size());
                std::iota(std::begin(x),std::end(x),0);
                name = "ensemble: " + std::to_string(i);
                plt::named_plot(name, x, phi, "-o");
        }
        plt::legend();

	i=1;
	plt::subplot(2,1,2);
	plt::ylim(-180,180);
	plt::xlabel("residues");
        plt::ylabel("psi");
                for(std::vector<std::vector<double>>::iterator it=psi_ensemble.begin(); it != psi_ensemble.end(); ++it,++i){
                psi = *it;
                std::vector<double> x(psi.size());
                std::iota(std::begin(x),std::end(x),0);
                name = "ensemble: " + std::to_string(i);
                plt::named_plot(name, x, psi, "-o");
        }
        plt::legend();

	plt::show();
//	std::string img_name = protein_name + ".eps";
//	plt::save("img_name");
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
 return torsion;
}

int main(int argc, char *argv[]){

	std::vector<double> phi;
	std::vector<double> psi;
	std::vector<std::vector<double>> ensemble_phi, ensemble_psi;

	if ( argc < 2 ){
		std::cout<<"Usage: "<<argv[0]<<" file1 file2 ...\n ";
		exit (EXIT_FAILURE);
	}
		

	for(int i=1;i<argc;i++){
		phi = getPhiOrPsi(argv[i],1);
		psi = getPhiOrPsi(argv[i],2);

		ensemble_phi.push_back(phi);
		ensemble_psi.push_back(psi);
	}

//	plotTorsionsEnsemble(ensemble_phi,"phi");
//	plotTorsionsEnsemble(ensemble_psi,"psi");
	plotPhiPsiSubplotEnsemble(ensemble_phi,ensemble_psi,"abc");
}
