function [FullAssemblages] = PermTwist(d,InitialAssem)
%This code takes either a:
%steering assemblage
%pseudo-measurement assemblage
%With 2 inputs and d outputs, in the form (d,d,2*d).
%and returns a new assemblage whose operators have been
%averaged over the automorphism group of d (prime).

%If an assemblage has the following: (e.g. output of XZTwist)
%input 1 has operators diagonal in the X basis with cycled eigenvalues.
%input 2 has operators diagonal in the Z basis with cycled eigenvalues.

%Then the returned operator will satisfy:
%input 1 has operators diagonal in the X basis with cycled eigenvalues, 
% where (d-1) of them co-incide.
%input 1 has operators diagonal in the Z basis with cycled eigenvalues, 
% where (d-1) of them co-incide.

%Full assemblage keeps track of the averaged assemblage.
FullAssemblages=zeros(d,d,2*d);
%One can check for each Unitary which operators are contributing, via permutations.
Contribute=zeros(2*d,d-1);


%We now create the automorphism unitaries: The permute the eigenvectors of
%both X and Z, leaving |z_0> and |x_0> invariant.
ShiftUnitaries=zeros(d,d,d-1);
ShiftPermutations=zeros(2*d,d-1);
for shift=0:d-2
for r=1:d
ShiftUnitaries(mod((r-1)+(r-1)*shift,d)+1,r,shift+1)=1;
ShiftPermutations(r,shift+1)=mod((r-1)+(r-1)*shift,d)+1;
ShiftPermutations(mod((r-1)+(r-1)*shift,d)+1+d,shift+1)=r+d;
end
end

%used to keep track of the contributions.
count=1;

%The size of the automorphism group is (d-1) (= |multiplicative group|).
for PVal=1:d-1
    
%We update the assemblage:
for i=1:2*d
CurrentAssem(:,:,i)=ShiftUnitaries(:,:,PVal)*InitialAssem(:,:,i)*ctranspose(ShiftUnitaries(:,:,PVal));
end
%And apply the induced permutation.
CurrentAssem=CurrentAssem(:,:,ShiftPermutations(:,PVal));


%We then add this contribution to the total:
FullAssemblages=FullAssemblages+CurrentAssem;

%Here we keep track of the applied permutation.
Order=[1:2*d];
Order=Order(ShiftPermutations(:,PVal));
%Contribute(:,count)=Order;
count=count+1;
end

%Finally we normalise the assemblage to preserve the trace.
FullAssemblages=FullAssemblages/(d-1);
end