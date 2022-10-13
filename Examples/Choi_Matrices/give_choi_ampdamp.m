function [ChoiAD] = give_choi_ampdamp(gamma)
%This function returns the Choi matrix of the amplitude damping channel: 
%a useful example for the channel dimension measure.

%Input:
%gamma: the damping parameter of the channel
%ChoiAD: the choi matrix of the amplitude damping channel.
ME=IsotropicState(2,1);
Kraus1=[[1,0];[0,sqrt(1-gamma)]];
Kraus2=[[0,sqrt(gamma)];[0,0]];
ChoiAD = kron(eye(2),Kraus1)*ME*kron(eye(2),transpose(Kraus1))+kron(eye(2),Kraus2)*ME*kron(eye(2),transpose(Kraus2));