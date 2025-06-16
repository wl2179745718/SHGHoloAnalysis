function S = ricbesj( nu, z )
%% RICBESJ calculates the Riccati-Bessel function of the first kind
% -------------------------------------------------------------------------
%% INPUT:
% -------------------------------------------------------------------------
% nu - the order of Ricatti-Bessel function
% z  - argument
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% S - Ricatti-Bessel function of the first kind
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
%% CALCULATING S
% -------------------------------------------------------------------------
S = sqrt( pi*z/2 ).*besselj( nu+0.5, z );
% -------------------------------------------------------------------------
end