function [ChoiErase] = give_choi_erase(ev)
%This function returns the Choi matrix of the erasure channel: 
%a useful example for the channel dimension measure.

%Input:
%ev - the probability of erasure.
%ChoiErase - the choi matrix of the erasure channel (applied to the second
%system)
ME=IsotropicState(2,1);
Kraus1=sqrt(1-ev)*[[1,0];[0,1];[0,0]];
Kraus2=sqrt(ev)*[[0,0];[0,0];[1,0]];
Kraus3=sqrt(ev)*[[0,0];[0,0];[0,1]];
ChoiErase = kron(eye(2),Kraus1)*ME*kron(eye(2),transpose(Kraus1))+kron(eye(2),Kraus2)*ME*kron(eye(2),transpose(Kraus2))++kron(eye(2),Kraus3)*ME*kron(eye(2),transpose(Kraus3));