function [ am, bm ] = sph_exp_coef_mie( m_max, wn, ref, rad )
%% SPH_EXP_COEF_MIE calculates expansion coefficients for sphere
% -------------------------------------------------------------------------
%% INPUT:
% -------------------------------------------------------------------------
% m_max  - cutoff for summation
% wn     - wavenumber, scalar, inverse microns
% ref    - refractive index of a sphere
% rad    - radius of a sphere, microns
% -------------------------------------------------------------------------
%% OUTPUT:
% -------------------------------------------------------------------------
% am, bm - expansion coefficients
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
        m = ( 1 : m_max );    
    x_vac = 2*pi*wn*rad;                                                    % size parameter 'in vacuum'
    x_par = x_vac*ref;                                                      % size parameter 'in sphere'
   Sn_vac =  ricbesj( m, x_vac );
  dSn_vac = dricbesj( m, x_vac );
   Sn_par =  ricbesj( m, x_par );
  dSn_par = dricbesj( m, x_par );
 xixn_vac =  ricbesh( m, x_vac );
dxixn_vac = dricbesh( m, x_vac );
% -------------------------------------------------------------------------
am = ( ref.*Sn_par.*dSn_vac   -        Sn_vac.*dSn_par )./ ...
     ( ref.*Sn_par.*dxixn_vac -      xixn_vac.*dSn_par );
bm = (      Sn_par.*dSn_vac   - ref.*  Sn_vac.*dSn_par )./ ...
     (      Sn_par.*dxixn_vac - ref.*xixn_vac.*dSn_par );
% -------------------------------------------------------------------
end