function [resi_region,part_tot_resi] = getRamachandranReigionDistri(dihed_angles)
%% Function to get the resi's falling in the various regions in the ramachandran map.
%  Calls the module isPinsidePoly() mex code.
%  Note that the region is as specified in the file ramachandranRegions.m
%  (part of Bioinfo toolbox) (Needs Bioinfo Toolbox to function.)
%  Input:  dihed_angles: k x 3 array containing the (phi,psi,chi) triplet
% Output:   resi_region: array of structure containing the fields
%                         region : 'Core Beta'
%                                   'Core R-Alpha'
%                                   'Core L-Alpha'
%                                   'Allowed 1'
%                                   'Allowed 2'
%                                   'Allowed 3'
%                                   'Allowed 4'
%                                   'Allowed 5'
%                                   'Allowed 6'
%                                   'Disallowed'
%                  (same as that of ramachandranRegions.m except the last)%                                  
%                           resi : resi falling in the region
%         part_tot_resi: part of the total resi available
%%          
   if nargin~=1
      error('Module: getRamachandranReigionDistri: Only 1 input needed.');
   end
   
   [r,c] = size(dihed_angles);
   if c~=3
      error('Module: getRamachandranReigionDistri: Incorrect input format.');
   end
%% get the informations about the various regions
   region_info = feval('ramachandranRegions');  % needs Bioinfo toolbox (code ramchandranRegions.m)  
                                               % fields Name,Color,Patch
   region_num  = size(region_info,2);
    
%%  begin checking for every residue with the dihedral angle and the 
%   bounded poly describing a particular region in ramachandran map
   resi_region = initializeResiRegion(extractfield(region_info,'Name'));   
   flag = zeros(1,r);  %0: no data 1: disallowed region 2: one of the region
   
   for i=1:r       
     if ~(isnan(dihed_angles(i,1)) || isnan(dihed_angles(i,2)))
        if(dihed_angles(i,1)<=180 && dihed_angles(i,1)>=-180) && ...
               (dihed_angles(i,2)<=180 && dihed_angles(i,2)>=-180) 
                 flag(i) =1;
                 for j=1:region_num
                   if isPinsidePoly(region_info(j).Patch,dihed_angles(i,1:2))
                      resi_region(j).resi =  [resi_region(j).resi,i];
                      flag(i) = 2;
                      break;   % so that resi 'i' isn't counted in more than 1 region
                              % to work strictly allowed regions must
                              % appear before partially allowed in
                              % region_info
                              
                   end
                 end          
        else
          fprintf('\n %d has incorrect dihedral \phi or \psi value.',i);
          error('Module: getRamachandranReigionDistri: Incorrect dihedral angle found.')
        end
     end
   end
   
   resi_region(end).resi = find(flag==1);
   part_tot_resi         = (r-numel(find(flag==0))) / r;
   
end

function resi_region = initializeResiRegion(region_names)

    resi_region = repmat(struct('Name',[],'resi',[]),length(region_names),1);
   for i=1:length(region_names)
      resi_region(i).Name = region_names{i};
      resi_region(i).resi = [];
   end
   
   resi_region(i+1).Name = 'Disallowed';
   resi_region(i+1).resi = [];
end