function [kspace,averages] = fill_kspace3D(ksp_4D, share)

% This function creates 2 arrays
% (1) the 3D kspace data sorted into the correct DCE frames and phase-encoding positions
% (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes

% Required input:
%
% ksp_4D                = sorted k-space data from MRecon, with dimensions:
    % [dimx, dimy, dimz, nCoils, nPhases]
    % dimz                  = 2nd phase encoding dimension
    % dimy                  = dimensions of the images: dimy (phase encoding)
    % dimx                  = dimensions of the images: dimx (readout)

sorted_kspace = ksp_4D;
sorted_kspace = permute(sorted_kspace, [4 5 3 2 1]);     % to match how Gustav does it
% [dimx, dimy, dimz, nCoils, nPhases, nResp] = size(sorted_kspace);
[nCoils,nPhases,dimz,dimy,dimx] = size(sorted_kspace);
sorted_averages = double(sorted_kspace ~= 0);

% temp new k-space copy
new_kspace = sorted_kspace;
new_averages = sorted_averages;

% find center of k-space
% to find lev, row, col, we average over cardiac and resp bins, and do a
% sum-of-squares over coils, then pick the point with high signal
tmp = mean(sorted_kspace,2);
tmp = squeeze(sqrt(sum(tmp.*tmp,1)));
% we only need dimx/2 point to get the point in our 2D (ky, kz) matrix
tmp = abs(squeeze(tmp));
[~,II] = max(tmp(:));
[lev, row, col] = ind2sub(size(tmp),II);

% Weighted view sharing
if (share > 0) && (nPhases > 1) 
    
    disp('View sharing ...');
    
    % respiratory of cardiac frames
    nrframes = nPhases;
    
    % determine share range. 
    % question: why work with both nPhases and nResp at the same time?
    maxshare = round(max([nPhases])/2); % maximum number of shares
    share(share > maxshare) = maxshare;
    weights = gauss(1:share+1,share,0);
    weights = weights/max(weights);
    
    % define ellipsoid regions
    % question:  why divide by share when this concerns the time
    % dimensions?
    Rz = round(dimz/share/2);
    Ry = round(dimy/share/2);
    Rx = round(dimx/share/2);
    [Z,Y,X] = ndgrid(1:dimz,1:dimy,1:dimx);
    for i = 1:share
        L(i,:,:,:) = sqrt( ((lev-Z)/(Rz*i)).^2 + ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
    end
    C(1,:,:,:) = L(1,:,:,:);
    if share > 1
        for i = 2:share
            C(i,:,:,:) = L(i,:,:,:) - L(i-1,:,:,:);
        end
    end
    
    % weights
    for i = 1:share
        for j = 1:share
            weights(i,j) = gauss(i+j-1,share,0);
        end
    end
    weights = 0.5*weights/max(weights(:));
    
    % apply sharing to k-space
    for frame = 1:nrframes
        disp(' ')
        disp(['Frame=' num2str(frame)])
        for i = -share:share
            disp(['share=' num2str(i)])
            sharedframe = frame + i;
            sharedframe(sharedframe < 1) = nrframes - sharedframe - 1;
            sharedframe(sharedframe > nrframes) = sharedframe - nrframes;
            
            if i~=0
                
                for j = 1:share
                    
                    ROI = repmat(reshape(squeeze(C(j,:,:,:)),[1 1 dimz dimy dimx 1])*weights(j,abs(i)),[nCoils 1 1 1 1]);
                    new_kspace(:,frame,:,:,:)   = new_kspace(:,frame,:,:,:)   + sorted_kspace(:,sharedframe,:,:,:)   .* ROI;
                    new_averages(:,frame,:,:,:) = new_averages(:,frame,:,:,:) + sorted_averages(:,sharedframe,:,:,:) .* ROI;
                    
                end
                
            end
            
        end
        
    end
end

% Normalize by number of averages
new_kspace = new_kspace./new_averages;            
new_kspace(isnan(new_kspace)) = complex(0);     % correct for NaN because of division by zero in case of missing k-lines
new_kspace(isinf(new_kspace)) = complex(0);

%% Apply a circular Tukey filter
filterwidth = 0.1;
flt = circtukey3D(dimz,dimy,dimx,lev,row,col,filterwidth);
tukeyfilter(1,1,:,:,:) = flt;

% Report back
kspace = ipermute(new_kspace.*tukeyfilter,[4 5 3 2 1]);
averages = ipermute(new_averages,[4 5 3 2 1]);

end