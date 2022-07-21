function save_mask( profile, mask )
%SAVE_MASK Summary of this function goes here
%   Detailed explanation goes here

save_mat(dir_out(profile.data_dir), [profile.name '.mask'], 'mask', mask);

end

