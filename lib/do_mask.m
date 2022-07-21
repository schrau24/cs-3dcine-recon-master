function do_mask(mrecon, filepath)
%DO_MASK Summary of this function goes here
%   Detailed explanation goes here

[path, filename] = fileparts(filepath);
mask = load_mat(path, filename, 'mask');

s = size(mrecon.Data);

% pad size to match case with averaging dimension
s = padarray(s, [0 12-length(s)], 1, 'post');

% previously 7:10, now trying to catch avg dim too

if length(size(mask)) == length(size(mrecon.Data))
	MASK = mask;
else
	if length(size(mask)) == 3
		% 2D+t
		% r : [y z t (1) (4:5) (7:12)]
		% p : [(1) y z (4:5) t (7:12)]
		rm = [1 1 1 s(1) s(4:5) s(7:12)];
		pm = [(1)+3 1 2 (4:5)+1 3 (7:12)];
		% pad mask to match (potential) oversampled size
		mask = padarray(mask, [s(2)-size(mask,1), s(3)-size(mask,2)]./2 );
	elseif length(size(mask)) == 4
		% 3D+t
		% r : [x y z t (4:5) (7:12)]
		% p : [x y z (4:5) t (7:12)]
		rm = [1 1 1 1 s(4:5) s(7:12)];
		pm = [1 2 3 (4:5)+1 4 (7:12)];
	end
	MASK = permute(repmat(mask,rm),pm);
end

%disp(size(mrecon.Data)); disp(size(mask)); disp(size(MASK)); % DEBUG
%mrecon.RemoveOversampling; MASK = padarray(MASK,[0 0 0 0 0 1],0,'post'); %
%SUPER HACKY HARDCODE FIX, hopefully obsolete

mrecon.Data = mrecon.Data .* MASK;

end
function value = load_mat(path, filename, variable )
%LOAD_MAT Summary of this function goes here
%   Detailed explanation goes here

str = load(fullfile(path, [filename '.mat']), variable);
value = str.(variable);

end