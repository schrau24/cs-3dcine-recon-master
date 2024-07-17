function kdataNew = preWhiten(filename, kdata)
MRn = MRecon(filename);
MRn.Parameter.Parameter2Read.typ = 5;
MRn.ReadData;
data_noise = MRn.Data;

data_ksp = kdata;
sizeMRDATA = size(data_ksp);
Ncoils = size(data_ksp,4);
Nsamples = sizeMRDATA(1)*sizeMRDATA(2)*sizeMRDATA(3);

psi = (1/(Nsamples-1))*(data_noise' * data_noise);
L = chol(psi,'lower');
L_inv = inv(L);

% loop over non-singleton dimensions (6, 10, 11, 12) to save memory
tmpdata = reshape(data_ksp,Nsamples,Ncoils,[]);
data_corr2 = zeros(size(tmpdata));
for ii = 1:size(tmpdata,3)
    data = tmpdata(:,:,ii).';
    data = conj(L_inv) * data;
    data_corr2(:,:,ii) = data.';
end
kdataNew = reshape(data_corr2,sizeMRDATA);
clear tmpdata data_corr2