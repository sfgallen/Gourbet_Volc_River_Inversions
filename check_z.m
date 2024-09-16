function z = check_z(S,z)
% this function makes sure the river doesn't flow backwards
%
% Author: Sean F. Gallen
% Contact: sean.gallen[at]colostate.edu
% Date modified: 09/28/2020
    Six = S.ix;                         % donors
    Sixc = S.ixc;                       % recievers
    
    % ma
    for lp = numel(Six):-1:1
        z_pre = z(Sixc(lp));
        z_cur = z(Six(lp));
        if z_cur <= z_pre
            z(Six(lp)) = z_pre+0.1;
        end
    end
end