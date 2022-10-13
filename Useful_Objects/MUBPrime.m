function [basisarray]=MUBPrime(p);
%This function calculates the p+1 mutually biased bases based via the Weyl
%construction, without requiring an extra MATLAB package which has the
%discrete Fourier transform:

%Inputs: 
%p, a prime number. 

%Outputs: 
%basisarray - a (p+1,p,p) array, where basisarray(n,k,:) gives the k^th
%basis vector of basis n. 






%First we define an array with indices (basis,vector,element)
basisarray=zeros(p+1,p,p);
if p==2
basisarray(1,1,:)=[1,0];
basisarray(1,2,:)=[0,1];
basisarray(2,1,:)=[1/sqrt(2),1/sqrt(2)];
basisarray(2,2,:)=[1/sqrt(2),-1/sqrt(2)];
basisarray(3,1,:)=[1/sqrt(2),1j/sqrt(2)];
basisarray(3,2,:)=[1/sqrt(2),-1j/sqrt(2)];
else
for cbasis=1:p
basisarray(p+1,cbasis,cbasis)=1;
end
for basisno=1:p
BMat=zeros(p,p);
for BMatR=1:p %gives basisvecno
for BMatC=1:p %givescompno
BMat(BMatR,BMatC)=exp((((BMatR-1)*(p-(BMatC-1)))-((basisno-1)*(((p-1)*(p)/2)-((BMatC-1)*(BMatC)/2))))*2*pi*1j/p);
end
end
BMat=BMat/sqrt(p);
basisarray(basisno,:,:)=BMat;
end
end
end

    

    
    
    
    

