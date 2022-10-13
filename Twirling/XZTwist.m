function [FullAssemblages] = XZTwist(d,InitialAssem)

%This code takes either a:
%steering assemblage
%pseudo-measurement assemblage
%With 2 inputs and d outputs, in the form (d,d,2*d).
%and returns a new assemblage (with the same trace), that has been averaged
%over the generalised Pauli X and Pauli Z symmetry for d(prime).

%This means:
%input 1 has operators diagonal in the X basis with cycled eigenvalues.
%output 2 has operators diagonal in the Z basis with cycled eigenvalues.


%Full assemblage keeps track of the averaged assemblage.
FullAssemblages=zeros(d,d,2*d);

%One can check for each X^iZ^j which operators are contributing, either via
%norm or via permutations.
Contribute=zeros(2*d,d^2);

%Unitaries X^i, and the corresponding permutations on the measurements
XUnitaries=zeros(d,d,d);
XPermutations=zeros(2*d,d);
%Unitaries Z^j, and the corresponding permutations on the measurements
ZUnitaries=zeros(d,d,d);
ZPermutations=zeros(2*d,d);
for Index=0:d-1
XUnitaries(:,:,Index+1)=GenPauli(Index,0,d,1);
XPermutations(1:d,Index+1)=[1:d];
XPermutations(d+1:2*d,Index+1)=mod([0:d-1]-Index,d)+d+1;
ZUnitaries(:,:,Index+1)=GenPauli(0,Index,d,1);
ZPermutations(1:d,Index+1)=mod([0:d-1]-Index,d)+1;
ZPermutations(d+1:2*d,Index+1)=[d+1:2*d];
end

%used to keep track of the contributions.
count=1;

%We run two loops: first we keep track of the power of X.
for XVal=1:d

%We update the assemblage.
for i=1:2*d
CurrentAssem(:,:,i)=XUnitaries(:,:,XVal)*InitialAssem(:,:,i)*ctranspose(XUnitaries(:,:,XVal));
end
%we then rearrange the assemblage via the permutation induced by the X
%unitary. 
CurrentAssemX=CurrentAssem(:,:,XPermutations(:,XVal));
%InitialAssem=CurrentAssem;

for ZVal=1:d
    
%We now update the assemblage using the Z unitary - note the use of
%CurrentAssemX, the post X unitary assemblage.
for k=1:2*d
CurrentAssem(:,:,k)=ZUnitaries(:,:,ZVal)*CurrentAssemX(:,:,k)*ctranspose(ZUnitaries(:,:,ZVal));
end
%The permutation then induced by the Z unitary is applied.
CurrentAssem=CurrentAssem(:,:,ZPermutations(:,ZVal));

%We now calculate the largest eigenvalue of each operator - this indicates
%where our original operators have been permuted to:
for j=1:2*d
Contribute(j,count)=norm(CurrentAssem(:,:,j),2);
end
%We then add this contribution to the total assemblage to be outputted:
FullAssemblages=FullAssemblages+CurrentAssem;

%Alternatively, one can use contribute to keep track of the overall
%permutation applied, for each X, Z pair.
Order=[1:2*d];
Order=Order(XPermutations(:,XVal));
Order=Order(ZPermutations(:,ZVal));
%Contribute(:,count)=Order;
count=count+1;

end
end
%Finally, the assemblage is normalised to preserve the trace:
FullAssemblages=FullAssemblages/d^2;
end




