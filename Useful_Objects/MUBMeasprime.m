function[MeasSet]=MUBMeasprime(p,weight)

%This function creates a measurement assemblages whose components are
%the MUB measurements in prime dimension mixed with white noise.

%Inputs:
%p: The dimension of the measurements (must be prime)
%weight: the visibility of the operators (high weight, low noise).

%Outputs: 
%MeasSet: A measurement assemblages of dimensions (p,p+1,p,p), ordered as
%outputs, inputs,dim,dim. It corresponds to a white noisy  MUB measurements (in prime dim).

MeasSet=zeros(p,p+1,p,p);
barr=MUBPrime(p);
for measurement=1:p+1
for outcome=1:p
MeasSet(outcome,measurement,:,:)=(weight)*conj(squeeze(barr(measurement,outcome,:)))*transpose(squeeze(barr(measurement,outcome,:)))+((1-weight)/p)*eye(p);
end
end



