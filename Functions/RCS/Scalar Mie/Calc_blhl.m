function [ blhl ] = Calc_blhl( l, wn, ref, rad )
%% SPH_EXP_COEF calculates expansion coefficients for sphere
% -------------------------------------------------------------------------
%% INPUT:
% -------------------------------------------------------------------------
% n_max  - cutoff for summation
% wn     - wavenumber, scalar, inverse microns
% ref    - refractive index of a sphere
% rad    - radius of a sphere, microns
% -------------------------------------------------------------------------
%% OUTPUT:
% -------------------------------------------------------------------------
% an, bn - expansion coefficients
% -------------------------------------------------------------------------
%% COPYRIGHT
% -------------------------------------------------------------------------
% Copyright 2018 Ilia Rasskazov, 
% University of Illinois at Urbana-Champaign &
% University of Rochester
% -------------------------------------------------------------------------
% Author:        Ilia Rasskazov
%                iliar@illinois.edu
%                irasskaz@ur.rochester.edu
% -------------------------------------------------------------------------
% Organizations: Beckman Institute for Advance Science and Technology,
%                University of Illinois at Urbana-Champaign 
%                http://beckman.illinois.edu/
%
%                The Institute of Optics, University of Rochester
%                http://www.hajim.rochester.edu/optics/
% -------------------------------------------------------------------------
%% PRE-ALLOCATING DATA
% -------------------------------------------------------------------------
                                                 % terms for summation
  x_vac = 2*pi*wn*rad;                                                      % size parameter 'in vacuum'
  x_par = x_vac*ref;                                                        % size parameter 'in sphere'
 jn_vac =  sbesselj( l, x_vac );
 hn_vac =  sbesselh( l, x_vac );
djn_vac = dsbesselj( l, x_vac );
dhn_vac = dsbesselh( l, x_vac );
 jn_par =  sbesselj( l, x_par );
djn_par = dsbesselj( l, x_par );
% -------------------------------------------------------------------------
%% CALCULATING EXPANSION COEFFICIENTS
alphal = Calc_alphal( l );
blhl = alphal.*( ref.*djn_par.*jn_vac -      jn_par.*djn_vac )./...
     (      dhn_vac.*jn_par - ref.*hn_vac.*djn_par ).*hn_vac;
% -------------------------------------------------------------------------
end