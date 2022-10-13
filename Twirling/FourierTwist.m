function [FullAssemblages] = FourierTwist(d,InitialAssem)
%This code takes either a:
%steering assemblage
%pseudo-measurement assemblage
%With 2 inputs and d outputs, in the form (d,d,2*d).
%and returns a new assemblage whose operators have been
%averaged over the fourier group of d (prime).

%If an assemblage has the following: (e.g. output of XZTwist)
%input 1 has operators diagonal in the X basis with cycled eigenvalues, 
% where (d-1) of them co-incide.
%input 1 has operators diagonal in the Z basis with cycled eigenvalues, 
% where (d-1) of them co-incide.

%Then the returned operator will satisfy:
%input 1 has operators diagonal in the X basis with cycled eigenvalues, 
% where (d-1) of them co-incide.
%input 1 has operators diagonal in the Z basis with cycled eigenvalues, 
% where (d-1) of them co-incide.

%The two sets of eigenvalues will co-incide.



%Full assemblage keeps track of the averaged assemblage.
FullAssemblages=zeros(d,d,2*d);

%One can check for each Unitary which operators are contributing, via permutations.
Contribute=zeros(2*d,d-1);

%The Fourier unitary has order 4: we now create these four unitaries:
%requires a package containing the discrete fourier transform.
FUnitaries=zeros(d,d,4);
FPermutations=zeros(2*d,4);
for FNum=0:3
FUnitaries(:,:,FNum+1)=mpower(dftmtx(d)/sqrt(d),FNum);
end
FPermutations(:,1)=[1:2*d];
FPermutations(:,2)=[d+1,flip(d+2:2*d),1:d];
FPermutations(:,3)=[1,flip(2:d),d+1,flip(d+2:2*d)];
FPermutations(:,4)=[d+1:2*d,1,flip(2:d)];

%used to keep track of the contributions.
count=1;

%We now run through the four possible Fourier matrix powers.
for FVal=1:4
    
%We update the assemblage:
for i=1:2*d
CurrentAssem(:,:,i)=FUnitaries(:,:,FVal)*InitialAssem(:,:,i)*ctranspose(FUnitaries(:,:,FVal));
end
%And apply the induced permutation.
CurrentAssem=CurrentAssem(:,:,FPermutations(:,FVal));

%We then add this contribution to the total assemblage to be outputted:
FullAssemblages=FullAssemblages+CurrentAssem;

%Here we keep track of the applied permutation.
Order=[1:2*d];
Order=Order(FPermutations(:,FVal));
Contribute(:,count)=Order;
count=count+1;
end

FullAssemblages=FullAssemblages/4;
end


