function imstack2d = imstackmontage(imstack3d)
% imstack2d = imstackmontage(imstack3d)
% imstack3d should be 3D only
%
% 2017-06-17, l.m.gottwald@amc.uva.nl

framesx = ceil(sqrt(size(imstack3d,3)));
framesy = ceil(sqrt(size(imstack3d,3)));
dimx = size(imstack3d,1);
dimy = size(imstack3d,2);
imstack2d = zeros(dimx*framesx,dimy*framesy,class(imstack3d));
slice = 1;
for ix = 1:framesx
    for iy = 1:framesy
        if slice <= size(imstack3d,3)
            imstack2d(1+(ix-1)*dimx:dimx*ix,1+(iy-1)*dimy:dimy*iy) = imstack3d(:,:,slice);
        end
        slice = slice+1;
    end
end
        