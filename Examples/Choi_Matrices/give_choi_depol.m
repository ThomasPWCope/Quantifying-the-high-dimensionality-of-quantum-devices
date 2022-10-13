function [ChoiDepol] = give_choi_depol(p,dim)
%This function returns the Choi matrix of the depolarising channel: 
%a useful example for the channel dimension measure.

%Input:
%p- The NOISE of the channel (not the visibility)
%dim - the dimension of the channel.

%Output:
%ChoiDepol - The choi matrix of the depolarising channel.

ChoiDepol = IsotropicState(dim,1-p);