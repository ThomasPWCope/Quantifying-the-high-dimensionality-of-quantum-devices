function [p,AssemOut]=See_Saw_Approach_Gen(d,MUBSet,k,startno,runno)
%This is code based on work from H. Chau Nyugen and Jonathan Steinberg. It
%attempts to find a lower bound for the visibility of a set of MUBs (as an
%assemblage)
%possible to simulate via lower dimensional assemblages.

%Inputs: 
%d  - the dimension (must be prime)
% MUBSet: the set of MUBs to be considered (a subset of 1:d+1)
%k: the allowed dimension to simulate with.
%startno: the number of initial starting points chosen.
%runno: the number of sew-saws the algorithm makes.

%Outputs p: a lower bound of the maximum possible visibility one can simulate.
%AssemOut: the lower rank assemblage which, when twirled over in the
%correct way, will reproduce the noisy MUB assemblage with visibility p.
Assem3=MUBAssemprime(d,1);

Assem3=Assem3(:,MUBSet,:,:);
largestvalue=0;
for run=1:startno
currentvalue=0;
lentuse=length(MUBSet);
EAssem=zeros(k,k,d,lentuse);
for povm=1:length(MUBSet)
RPOVM=RandomPOVM(k,d);
%RPOVM2=RandomPOVM(k,d);
for out=1:d
EAssem(:,:,out,povm)=RPOVM{out};
%EAssem(:,:,out,2)=RPOVM2{out};
end
end
for i=1:runno
Operator=zeros(d*k,d*k);
for xvcount=1:length(MUBSet)
    xv=MUBSet(xvcount);
    for av=1:d
       Operator=Operator+kron(squeeze(EAssem(:,:,av,xvcount)),squeeze(Assem3(av,xvcount,:,:)));
    end
end


[T,V]=eigs(Operator);
MaxVector=T(:,1);
MaxState=MaxVector*ctranspose(MaxVector);
Operator=zeros(d*k,d*k);
numberin=length(MUBSet);
cvx_begin sdp quiet
variable Enew(k,k,d,numberin) Hermitian semidefinite
for xvcount=1:numberin
    xv=MUBSet(xvcount);
    for av=1:d 
       Operator=Operator+kron(Enew(:,:,av,xvcount),squeeze(Assem3(av,xvcount,:,:)));
    end
end
maximise real(trace(Operator*MaxState))
subject to
for xvcount=1:numberin
sum(Enew(:,:,:,xvcount),3)==eye(k);
end
cvx_end
EAssem=Enew;
currentvalue=cvx_optval;
end
if currentvalue>largestvalue
    largestvalue=currentvalue;
    BestState=MaxState;
    BestAssem=EAssem;
end
end
largestvalue;
p=((largestvalue)-(length(MUBSet)/d^2))/((length(MUBSet)/d)-(length(MUBSet)/d^2));
AssemOut=zeros(d,d,d,numberin);
for xvcount=1:length(MUBSet)
    xv=MUBSet(xvcount);
    for av=1:d
        AssemOut(:,:,av,xvcount)=PartialTrace(MaxState*kron(EAssem(:,:,av,xvcount),eye(d)),1,k,d);
    end
end





