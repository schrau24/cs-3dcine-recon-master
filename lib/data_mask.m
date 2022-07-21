function data = data_mask(data,NoiseClipValue)
% data = data_mask(data)
%
% mask the data based on a noise level

% take mean over entire data set (x,y,z, : )
dims = size(data);
th_data = mean(abs(reshape(data,[dims(1:3),prod(dims(4:end))])),4);

% make mask by threshold / NoiseClipValue
th_mask = th_data > NoiseClipValue;

% apply mask to data
th_mask_all = repmat(th_mask,[1,1,1,dims(4:end)]);
data = data .* th_mask_all;
end