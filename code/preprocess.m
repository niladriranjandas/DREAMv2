%%--------------------------Our algo---------------------------------------
%    preprocess: called by resiwise_distributed.m
%             read the i/p files, create the upper, lower bounds
%             calls preprocesses to create Comp and wh_Comp structures
%--------------------------------------------------------------------------
%%

addpath('helperfunctions');
addpath('inputreaders');
addpath('sdpstuff');
addpath('refinement');
addpath('../packages/hanso');

%% ----------------- set the inputs and parameters ------------------------

load '../data/ainfo.mat'

% Drop side-chain hydrogen atoms or not
%hydrogen_omission = 1;%0;%1;


% Number of BFGS iterations
gd_tol  = 10^-9;
cg_iter = 500;
f = [10 10 10 10 10]; % [f_hb, f_tau, f_tal, f_vdw, f_tas]

fprintf('==========================================================================\n')
fprintf('------------Our algo:Protein Structure Determination----------------------\n')

tstart = tic;

protein_name = INPUTS.protein_name;
in_max_res = [];
in_min_res = [];

in_seq_file   = INPUTS.seq_file ;
in_upl_file   = INPUTS.upl_file ;
in_hbond_file = INPUTS.hbond_file;
in_ang_file   = INPUTS.ang_file ;

%% ----------- read the i/p files

protein_path = ['../protein/' protein_name '/'];
fprintf('*************************************************************************\n');
fprintf('Protein: %s\n', protein_name);
fprintf('-Reading input files...\n');
% Reading input data
%==========================================================================
seq_file = [protein_path in_seq_file];
[seq, num] = seq_reader(seq_file);
   org_num = num; %added by me for prepForFillGap
max_res = max(num);
min_res = min(num);
    
num_upl = size(in_upl_file, 2);
upl_file = cell(1, num_upl);
for i = 1:num_upl
    upl_file{i} = [protein_path in_upl_file{i}];
end
if ~isempty(in_hbond_file)
    hbond_file = [protein_path in_hbond_file];
    hbond_write_file = [protein_path protein_name '_hbo.upl'];
    %hbond_reader(hbond_file,hbond_write_file);
    hbond_reader_me(hbond_file,hbond_write_file);
    num_upl = num_upl + 1;
    upl_file{num_upl} = hbond_write_file;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
raw_up    = cell(1, num_upl);
raw_up_ho = cell(1, num_upl);
temp_min_res = +inf;
temp_max_res = -inf;
for i = 1:num_upl
    raw_up{i} = dist_reader(upl_file{i}, A);
    if hydrogen_omission
        raw_up_ho{i} = dist_reader(upl_file{i}, A, hydrogen_omission);
    end
    temp_min_res = min(temp_min_res, min([raw_up{i}.tres raw_up{i}.sres]));
    temp_max_res = max(temp_max_res, max([raw_up{i}.tres raw_up{i}.sres]));
end
if isempty(in_min_res)
    min_res = max(temp_min_res-1, min_res);
else
    min_res = in_min_res;
end
if isempty(in_max_res)
    max_res = min(temp_max_res+1, max_res);
else
    max_res = in_max_res;
end
    
% Remove informationless (w/o any constraints)
% parts from and N- and C-terminus
ind_del_N = num < min_res;
ind_del_C = num > max_res;
ind_del   = ind_del_N | ind_del_C;
seq(ind_del) = [];
num(ind_del) = [];
    
% dihedral angle constraints
ang_file = [protein_path in_ang_file];
[phi_cons, psi_cons] = ang_reader(ang_file, num);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pdb_write_file = [protein_path protein_name '.pdb'];
%==========================================================================
tread = toc(tstart);
fprintf('\tdone: %4.1f sec\n', tread)

fprintf('-Sampling a random molecule...\n')
% Generating a random structure
%==========================================================================
[phi, psi] = ang_sampler(seq,phi_cons,psi_cons);
switch hydrogen_omission
    case 0
        [rand_X, Comp] = ibuildprot(seq,num,phi,psi,A);
        wh_rand_X = rand_X;
        wh_Comp   = Comp;
    case 1
        [wh_rand_X, wh_Comp, rand_X, Comp] = ibuildprot(seq,num,phi,psi,A);
end

if exist(ang_file, 'file')
    [ang_lo_cons, ang_up_cons] = ang_dist_conmaker(phi_cons,psi_cons,wh_Comp);
end
%==========================================================================
trand = toc(tstart);
fprintf('\tdone: %4.1f sec\n',trand - tread)

fprintf('-Reducing the SDP problem...\n')

% Perform the reduction
%==========================================================================
switch hydrogen_omission
    case 0
        [U, Comp.cliq_dims] = reducer(rand_X,Comp);
        wh_Comp = Comp;
    case 1
        dont_compute_U = 1;
        [U, Comp.cliq_dims] = reducer(rand_X,Comp);
        [wh_U, wh_Comp.cliq_dims] = reducer(wh_rand_X,wh_Comp,dont_compute_U);
end
%==========================================================================
treduc = toc(tstart);
fprintf('\tdone: %4.1f sec\n',treduc - trand)
fprintf('-Forming constraints...\n')
% Generating upper and lower bounds constraints
%==========================================================================
start = 0;
wh_up_bounds = nan(50000,4);
for i = 1:num_upl
    temp_upl = upper_maker(raw_up{i}, wh_Comp);
    wh_up_bounds(start+1:start+size(temp_upl,1),:) = temp_upl;
    start = start + size(temp_upl,1);
end
wh_up_bounds(isnan(wh_up_bounds(:,1)), :) = [];

if hydrogen_omission
    start = 0;
    ho_up_bounds = nan(50000,4);
    for i = 1:num_upl
        temp_upl = upper_maker(raw_up_ho{i}, wh_Comp);
        ho_up_bounds(start+1:start+size(temp_upl,1),:) = temp_upl;
        start = start + size(temp_upl,1);
    end
    ho_up_bounds(isnan(ho_up_bounds(:,1)), :) = [];
end

% adding torsion-angle constraints
if exist(ang_file, 'file')
    wh_up_bounds = [wh_up_bounds; ang_up_cons];
    wh_sdp_lo_bounds = ang_lo_cons;
end
    
switch hydrogen_omission
    case 0
        % equality cons
        wh_eq_cons = equality_con_former(rand_X,Comp);
        eq_cons = wh_eq_cons;
        % upper bounds
        up_bounds = wh_up_bounds;
        % vdw bounds
        wh_vdw_bounds = vdw_bound_maker(wh_Comp);
        vdw_bounds    = wh_vdw_bounds;
        sdp_lo_bounds = wh_sdp_lo_bounds;
    case 1
        % equality cons
        wh_eq_cons = equality_con_former(wh_rand_X,wh_Comp);
        eq_cons    = equality_con_former(rand_X,Comp);
        % upper bounds
        up_bounds = map_bounds(ho_up_bounds,Comp.atoms_map);
        if exist(ang_file,'file')
            ang_up_bounds = map_bounds(ang_up_cons,Comp.atoms_map);
            up_bounds = [up_bounds; ang_up_bounds];
        end
        % vdw bounds
        wh_vdw_bounds = vdw_bound_maker(wh_Comp);
        vdw_bounds    = vdw_bound_maker(Comp);
        sdp_lo_bounds = map_bounds(ang_lo_cons,Comp.atoms_map);
        % a bug that need to be fixed sometimes
        % when TYR is HH it is mapped to its OH
        up_bounds(:,4) = wh_up_bounds(:,4);
end
    
lo_bounds     = vdw_bounds;
wh_lo_bounds  = wh_vdw_bounds;
if ~exist('sdp_lo_bounds', 'var')
    sdp_lo_bounds = [];
else
    lo_bounds    = [lo_bounds; sdp_lo_bounds];
    wh_lo_bounds = [wh_lo_bounds; wh_sdp_lo_bounds];
end
%==========================================================================
tbounds = toc(tstart);
fprintf('\tdone: %4.1f sec\n',tbounds - treduc)

%==========================================================================
if omit.percent ~= 0
    %omit.percent = 30;
    omit.num_upbounds = round((omit.percent/100)*length(up_bounds));
    omit.lo_indx = randi(length(up_bounds),1, omit.num_upbounds);

    up_bounds_bkp = up_bounds;
    fprintf('\n Deleted %d bounds from up bounds ',omit.num_upbounds);
    omit.include_indx = setdiff(1:length(up_bounds),omit.lo_indx);
    
    up_bounds = up_bounds(omit.include_indx,:);
end
%up_bounds = up_bounds(omit.include_indx,:);
