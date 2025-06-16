function [E_tot,E_sca] = tot2sca(E_simu,U_inp_end,opt)
%UNTITLED4 Summary of this function goes here

switch opt
    case 'Vol'
E_tot = squeeze(E_simu(:,:,end));
    case 'out'
E_tot = E_simu;
end

E_inc = U_inp_end;
E_sca = E_tot-E_inc;

end

