function [bound] =RobustnessBound(d,n)
% For a pair of prime dimension MUBS d, returns the upper bound on the 
%possible (white noise) visibility that one could create with rank n 
%pseduo-measurements. Based on the bound from:
%https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.126.200404 

%Inputs:
%d : a PRIME NUMBER dimension.
%n : the compression dimension one is allowed to use to simulate the MUB
%pair.

%Outputs:
% bound: the maximum visibility allowed by the bound of the maximum
% robustness for a fixed compression dimension.


bound = ((d + sqrt(d)-1)*sqrt(n) - 1)/((d - 1)*(sqrt(n) + 1));