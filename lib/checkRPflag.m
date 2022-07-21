function out = checkRPflag( RP, field )
% check if profile parameter exists and is true
% 'field' needs to be a string

out = isfield(RP,field);
    if out
        eval(sprintf('out = (RP.%s == 1);',field));     
    end
end