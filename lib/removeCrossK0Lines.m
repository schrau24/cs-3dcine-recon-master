function [mrecon] = removeCrossK0Lines(mrecon)

%% Define mrecon points to be used
noisemeas = sum(mrecon.Parameter.Labels.Index.typ == 5);
phasecorrmeas = sum(mrecon.Parameter.Labels.Index.typ == 3);
data_in_use_index = false(size(mrecon.Parameter.Labels.Index.rtop));
data_in_use_index(1+noisemeas+phasecorrmeas:mrecon.Parameter.Labels.OriginalLabelLength)=true;  %rm: you already only read in data of typ ==1?

% channels
chan=mrecon.Parameter.Labels.CoilNrs(:,1);
nchan=length(chan);
% flow encoding directions
extr1=mrecon.Parameter.Parameter2Read.extr1;
tmp =mrecon.Parameter.Labels.Index.ky(data_in_use_index)==0 & mrecon.Parameter.Labels.Index.kz(data_in_use_index) == 0 & mrecon.Parameter.Labels.Index.chan(data_in_use_index) == chan(1) & mrecon.Parameter.Labels.Index.extr1(data_in_use_index)== extr1(1);
sg_sampling_index = find(tmp);

% k0 indices for all coils
k0_m_index = []; k0_p_index = []; k0_s_index = [];
count = 0;
for ii=1:length(sg_sampling_index)
    ind = (sg_sampling_index(ii)-1)+1 : (sg_sampling_index(ii)-1)+nchan;
    count = count + 1;
    
    % also the k0_m, k0_p, and k0_s indices
    switch count
        case 1  % measurement direction
            k0_m_index = cat(1,k0_m_index,ind');
            
        case 2  % phase direction
            k0_p_index = cat(1,k0_p_index,ind');
            
        case 3  % slice direction
            k0_s_index = cat(1,k0_s_index,ind');
    end
    
    if count == 3
        count = 0;
    end
end


%     % set k0_p and k0_s to typ=2
startInd = noisemeas+phasecorrmeas;
    mrecon.Parameter.Labels.Index.typ(startInd+sort(cat(1,k0_p_index,k0_s_index))) = 2;

% % save startInd, typ==2, typ==4 indices
% ind = mrecon.Parameter.Labels.Index;
% typ2 = find(ind.typ==2) - startInd;
% typ4 = sort(cat(2,k0_p_index,k0_s_index));
% 
% % set k0_p and k0_s to typ=2
% mrecon.Parameter.Labels.Index.typ(ind.typ==2) = 1;
end