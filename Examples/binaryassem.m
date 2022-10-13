function [assemout]=binaryassem(d)
%binaryassem: for a given dimension (d), constructs a steering assemblage whose
%true Schmidt measure is log2(2)=1, but the upper bound given by the greedy
%algorithm (steering_schmidt_meas) will be log2(d). 
assemout=zeros(d,2,d,d);
for i=1:d
    assemout(i,1,i,i)=1/d;
    plusvec=zeros(d,1);
    plusvec(mod(i-1,d)+1)=1/sqrt(2);
    plusvec(mod(i,d)+1)=1/sqrt(2);
    minusvec=zeros(d);
    minusvec(mod(i,d)+1)=1/sqrt(2);
    minusvec(mod(i+1,d)+1)=-1/sqrt(2);
    assemout(i,2,:,:)=(plusvec*transpose(plusvec)+minusvec*transpose(minusvec))/(2*d);
end

    
    
