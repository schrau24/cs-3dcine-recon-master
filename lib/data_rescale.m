function data = data_rescale(data)
% out = data_rescale(in)
% 
% scale data to 2^12 = 4096

datamax = ceil(max(abs(data(:))));
data = data ./ datamax .* 2^12;
end
