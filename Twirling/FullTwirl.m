function [p,Assemblage_Out]=FullTwirl(d,Assemblage)
%This code takes either a:
%steering assemblage
%pseudo-measurement assemblage
%With 2 inputs and d outputs, in the form (d,d,2*d).
%and returns a new assemblage (with the same trace), that consists of the
%canonical X,Z MUBs mixed with white noise. 

%p: The visibility parameter
%Assemblage_Out: the twirled assemblage.

AssemXZ=XZTwist(d,Assemblage);
AssemP=PermTwist(d,AssemXZ);
Assemblage_Out=FourierTwist(d,AssemP);
Marginal=sum(Assemblage_Out(:,:,1:d),3);

if abs(trace(Marginal)-1)<=10^-5
    p = d*(Assemblage_Out(1,1,d+1)-Assemblage_Out(d,d,d+1));
elseif abs(trace(Marginal)-d)<=10^-5
    p = (Assemblage_Out(1,1,d+1)-Assemblage_Out(d,d,d+1));
else
    Assemblage_Out=d*Assemblage_Out/(trace(Marginal));
    p = (Assemblage_Out(1,1,d+1)-Assemblage_Out(d,d,d+1));
end
end
    
    



