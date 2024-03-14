function [upl_densities, upl_densities_strong, upl_densities_medium, upl_densities_weak, outfile] = subGraphDensityFromUpl(protein_name, in_seq_file, in_upl_file, in_hbond_file, in_ang_file)
%%
%[upl_densities, outfile] = subGraphDensityFromUpl('5x1x_ambi5r1_noe','5x1x.seq', {'5x1x_ambi5r1_noe_concatupl.upl'}, '', '5x1x_concat_dihed.aco');
%%	
 
    STRONG = 2.7
    MEDIUM = 3.5
    WEAK   = 6
	outfile = protein_name + ".csv"

	INPUTS.protein_name  = protein_name ;
	INPUTS.seq_file = in_seq_file  ;
	INPUTS.upl_file = in_upl_file  ;
	INPUTS.hbond_file = in_hbond_file;
	INPUTS.ang_file = in_ang_file  ;

	hydrogen_omission = 0;
	omit.percent = 0;
	LOG.div_n_conquer = "abc"
    
    [TOL_DENSITY_STRONG, TOL_DENSITY_MEDIUM, TOL_DENSITY_LOW, DESTDIR] = setParams() % set the tolerances for density_ij check
%%	
	%preprocess  ---- preprocess already run in div_conquer_preloc_test
	div_conquer_preloc_test
%%	

	num_bounds = length(raw_up{1});
	upl_densities = zeros(num_bounds,10);	

	fp = fopen(outfile,'w');
	fprintf(fp,'uplindex,resii,resij,expdist,densityi,densityj,densityij')
	for l=1:num_bounds
		resi_i = raw_up{1}(l).sres;
		resi_j = raw_up{1}(l).tres;
		dist_ij = raw_up{1}(l).dist;
        %test.ret_densities
        %test.CI_test_resi
        resi_i
        resi_j
        test
        [density_i, density_j, density_ij] = getSubgraphDensity(resi_i, resi_j, test.ret_densities, test.cI_test_resi);
        
        fprintf(fp,'\n%d,%d,%d,%f,%f,%f,%f', l, resi_i, resi_j, dist_ij, density_i, density_j, density_ij);
        upl_densities(l,:) = [double(l), double(resi_i), double(resi_j), dist_ij, density_i, density_j, density_ij, abs(density_i-density_ij), abs(density_j-density_ij), max(abs(density_i-density_ij), abs(density_j-density_ij))];
	end
	fclose(fp)
	   
%%
	%upl_densities(:,8)  = abs( upl_densities(:,5) - upl_densities(:,7));
	%upl_densities(:,9)  = abs( upl_densities(:,6) - upl_densities(:,7));	
	%upl_densities(:,10) = max( upl_densities(:,8), upl_densities(:,9));

%% 
	upl_densities_strong = upl_densities(upl_densities(:,4)>0 & upl_densities(:,4)<=STRONG,:);
	upl_densities_medium = upl_densities(upl_densities(:,4)>STRONG & upl_densities(:,4)<=MEDIUM,:);
	upl_densities_weak   = upl_densities(upl_densities(:,4)>MEDIUM & upl_densities(:,4)<=WEAK,:);

    strong_diff_no0  = upl_densities_strong(:,end);
    strong_diff_no0(strong_diff_no0==0) = NaN;
    strong_diff_std  = nanstd(strong_diff_no0);
    strong_diff_mean = nanmean(strong_diff_no0);
    upl_densities_strong(:,end+1) = upl_densities_strong(:,end) >= (strong_diff_mean - strong_diff_std);
    strong_ij = upl_densities_strong(:,7);
    strong_ij(strong_ij==0) = NaN;
    strong_ij_mean = nanmean(strong_ij);
    strong_ij_std  = nanstd(strong_ij);
    strong_density_ij_thresshold = min(TOL_DENSITY_STRONG, strong_diff_mean - strong_diff_std);
    upl_densities_strong(:,end+1) = upl_densities_strong(:,7) <= strong_density_ij_thresshold;
    upl_densities_strong(:,end+1) = upl_densities_strong(:,end-1) | upl_densities_strong(:,end);
    writeUpls(upl_densities_strong, strcat(DESTDIR,'/',protein_name,'_strong.csv'));

    medium_diff_no0  = upl_densities_medium(:,end);
    medium_diff_no0(medium_diff_no0==0) = NaN;
    medium_diff_std  = nanstd(medium_diff_no0);
    medium_diff_mean = nanmean(medium_diff_no0);
    upl_densities_medium(:,end+1) = upl_densities_medium(:,end) >= (medium_diff_mean - medium_diff_std);
    medium_ij = upl_densities_medium(:,7);
    medium_ij(medium_ij==0) = NaN;
    medium_ij_mean = nanmean(medium_ij);
    medium_ij_std  = nanstd(medium_ij);
    medium_density_ij_thresshold = min(TOL_DENSITY_MEDIUM, medium_diff_mean - medium_diff_std);
    upl_densities_medium(:,end+1) = upl_densities_medium(:,7) <= medium_density_ij_thresshold;    
    upl_densities_medium(:,end+1) = upl_densities_medium(:,end-1) | upl_densities_medium(:,end);
    writeUpls(upl_densities_medium, strcat(DESTDIR,'/',protein_name,'_medium.csv'));

    weak_diff_no0  = upl_densities_weak(:,end);
    weak_diff_no0(weak_diff_no0==0) = NaN;
    weak_diff_std  = nanstd(weak_diff_no0);
    weak_diff_mean = nanmean(weak_diff_no0);
	upl_densities_weak(:,end+1) = upl_densities_weak(:,end) >= (weak_diff_mean - weak_diff_std);
    weak_ij = upl_densities_weak(:,7);
    weak_ij(weak_ij==0) = NaN;
    weak_ij_mean = nanmean(weak_ij);
    weak_ij_std  = nanstd(weak_ij);
    weak_density_ij_thresshold = min(TOL_DENSITY_LOW, weak_diff_mean - weak_diff_std);
    upl_densities_weak(:,end+1) = upl_densities_weak(:,7) <= weak_density_ij_thresshold;
    upl_densities_weak(:,end+1) = upl_densities_weak(:,end-1) | upl_densities_weak(:,end);
    writeUpls(upl_densities_weak, strcat(DESTDIR,'/',protein_name,'_weak.csv'));

end

function [TOL_DENSITY_STRONG, TOL_DENSITY_MEDIUM, TOL_DENSITY_LOW, DESTDIR] = setParams()
%%
%   set parameters for the tolerance values
%       checked: min(sample_mean - sample_std, tolerance) for various categories
%%
  
        TOL_DENSITY_STRONG = 0.2;
        TOL_DENSITY_MEDIUM = 0.25;
        TOL_DENSITY_LOW    = 0.3;

        DESTDIR = 'ambiguous_bounds';
end

function writeUpls(ipmatrix, opfile)
%% 

%%
	[r,c] = size(ipmatrix);
	op = fopen(opfile,'w');
	fprintf(opfile)

%
	fprintf(op,'uplindex,resii,resij,distij,densityi,densityj,densityij,diffdensity_i_ij,diffdensity_j_ij,diffdensity_max_i_j,diffdensitycheck,density_ij_check,density_together');
    for i=1:r
    	fprintf(op,'\n%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d',ipmatrix(i,1),ipmatrix(i,2),ipmatrix(i,3),ipmatrix(i,4),ipmatrix(i,5),ipmatrix(i,6),ipmatrix(i,7),ipmatrix(i,8),ipmatrix(i,9),ipmatrix(i,10),ipmatrix(i,11),ipmatrix(i,12),ipmatrix(i,13));
    end	
    fclose(op);

end