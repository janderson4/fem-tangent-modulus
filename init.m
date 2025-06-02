function [SIG,Ce,ALPHA,R,mu,H]=init(sigma_y0, E, ET, nu, n_int)
% FUNCTION TO INITIALIZE VARIABLES
% Args:
%   sigma_y0: the initial yield stress
%   E: the elastic modulus
%   ET: the plastic modulus
%   nu: the Poisson ratio
%   n_int: the number of integration points
% Returns:
%   SIG: the initialized stress at all integration points (6 x n_int)
%   Ce: the initialized elastic tangent operator (6 x 6)
%   ALPHA: the initialized alpha at all integration points (6 x n_int)
%   SIG: the initialized R at all integration points (1 x n_int)
%   SIG: the initialized stress at all integration points (6 x n_int)
%   SIG: the initialized stress at all integration points (6 x n_int)

SIG=zeros(6,n_int);
ALPHA=zeros(6,n_int);

% INITIALIZE C

bulk=E/(3*(1-2*nu));
mu=E/(2*(1+nu));

% elastic tangent tensor assuming tensorial shear strain
Ce=2*mu*eye(6);
Ce(1:3,1:3)=Ce(1:3,1:3)+(bulk-2/3*mu).*ones(3);

% INITIALIZE R

R = sigma_y0*sqrt(2/3)*ones(1,n_int);

H_prime=ET/(1-ET/E);

H=2/3*H_prime;
   
return
