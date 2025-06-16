function j = sbesselj( nu, z )
%% SBESSELJ calculates spherical Bessel function of the first kind
% -------------------------------------------------------------------------
%% INPUT:
% -------------------------------------------------------------------------
% nu - the order of spherical Bessel function
% z  - argument
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% j - spherical Bessel function of the first kind
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
%% CALCULATING j
% -------------------------------------------------------------------------
j = ricbesj( nu, z )./z;
% -------------------------------------------------------------------------
end