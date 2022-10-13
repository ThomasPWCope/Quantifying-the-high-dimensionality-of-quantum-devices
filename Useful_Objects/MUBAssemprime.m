function[assemblage]=MUBAssemprime(p,weight)
%This function creates a steering assemblages whose components are
%the MUB measurements in prime dimension applied to a convex mixture of
%the maximally entangled and maximally mixed state.

%Inputs:
%p: The dimension of the measurements (must be prime)
%weight: the convex weight of the maximally entangled state in the mixture.

%Outputs:
%Assemblages: A steering assemblages of dimensions (p,p+1,p,p), ordered as
%outputs, inputs,dim,dim. It corresponds to a white noisy maximally
%entangled state measured with MUB measurements (in prime dim).

state=MaxEntangled(p)*transpose(MaxEntangled(p));
assemblages=zeros(p,p+1,p,p);
barr=MUBPrime(p);
for measurement=1:p+1
for outcome=1:p
assemblage(outcome,measurement,:,:)=weight*PartialTrace(kron(conj(squeeze(barr(measurement,outcome,:)))*transpose(squeeze(barr(measurement,outcome,:))),eye(p))*state,1,p)+(1-weight)*eye(p)/p^2;
end
end



