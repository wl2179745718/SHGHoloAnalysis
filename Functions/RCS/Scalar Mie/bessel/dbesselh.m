function dH = dbesselh( nu, k, z )
%% DBESSELH calculates first derivative of the Bessel function of the third kind
% -------------------------------------------------------------------------
%% INPUT:
% -------------------------------------------------------------------------
% nu - the order of Ricatti-Bessel function
% k  - 1 or 2, default = 1
% z  - argument
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% dH - first derivative of the Bessel function of the third kind
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
%% CHECKING INPUT
% -------------------------------------------------------------------------
if nargin == 2
    z = k;
    k = 1;
end
% -------------------------------------------------------------------------
%% CALCULATING dH
% -------------------------------------------------------------------------
dH = 0.5*( besselh( nu-1, k, z ) - besselh( nu+1, k, z) );
% -------------------------------------------------------------------------
end