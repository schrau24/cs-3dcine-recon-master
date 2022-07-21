function k = trajectoryForMotionCorrection(mrecon)
% build the 3D coordinates for a PROUD scan, to be use for
% motion/respiratory correction
% output: k, nFE x nReadouts x 3, normalized to +/- 0.5
% Eric Schrauben, Amsterdam UMC, location AMC, 2021-03-15

% first create and normalize k-space locations for x,y,z directions
nFE = size(mrecon.Data,1);
kxRange = mrecon.Parameter.Encoding.KxRange;
kxVector = linspace(kxRange(1),kxRange(2), nFE);
kx = repmat(kxVector,[size(mrecon.Data,2) 1])';
kx = 0.5*kx/max(abs(kx(:)));    % normalize
maxKyRange = double(max(abs(mrecon.Parameter.Encoding.KyRange)));
ky = double(repmat(mrecon.Parameter.Labels.Index.ky(mrecon.Parameter.Labels.Index.typ==1),[1 nFE])');
ky = 0.5*ky/maxKyRange;         % normalize to true kspace extent (with zero-filling)
maxKzRange = double(max(abs(mrecon.Parameter.Encoding.KzRange)));
kz = double(repmat(mrecon.Parameter.Labels.Index.kz(mrecon.Parameter.Labels.Index.typ==1),[1 nFE])');
kz = 0.5*kz/maxKzRange;         % normalize to true kspace extent (with zero-filling)
% build 3D kspace trajectory to correct
k = cat(3,kx,ky,kz); clear kx ky kz