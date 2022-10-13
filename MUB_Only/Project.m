function [InitialAssem,Proj] = Project(d,n,varargin)
%This function takes:
%d: a prime dimension,
%n: a number <=d, a rank constraint.
%varargin:
% a length n vector of numbers <=d - the choice of non-zero outputs for
% measurement 1
%a length n vector of numbers <=d - the choice of non-zero outputs for
% measurement 2

%The function outputs a steering assemblage, which we conjecture is the
%best possible for constructing white noise MUBs from low dimensions.

%It can also output the projector, which is used to ensure the rank is
%correct.

%Technical Note: Although the output of this function is an assemblage, 
%this construction is actually aimed to reproduce the set of noisy MUB
%measurements. since the the noisy MUB steering assemblage and noisy MUB
%measurement assemblage are equivalent up to transposition; one can convert
%this construction to the steering scenario without too much trouble - the
%twirling functions should still work; and the unitary set is closed under
%transposition.


% If no choices are supplied, we choose the first k outputs;
%this is known NOT to be optimal for d=5;
if nargin==2
    firstchoice=[1:n];
    secondchoice=[1:n];
elseif nargin==3
    firstchoice=varargin{1};
    secondchoice=[1:n];
elseif nargin==4
    firstchoice=varargin{1};
    secondchoice=varargin{2};
else 
    error('Wrong number of variables!')
end
    
    
%We first choose our MUBs, in prime dimensions.
PrimeMUB=MUBPrime(d);
Marginal=zeros(d,d);

%Our construction begins by totalling the chosen output operators.
for no=1:n
no1=firstchoice(no);
no2=secondchoice(no);
    Marginal=Marginal+ctranspose(PrimeMUB(1,:,no1))*PrimeMUB(1,:,no1);
    Marginal=Marginal+ctranspose(PrimeMUB(d+1,:,no2))*PrimeMUB(d+1,:,no2);
end

%We then make the marginal numerically more stable:
Marginal=(Marginal+ctranspose(Marginal))/2;

%Subsequently, we take the largest k eigenvectors and project onto these.
[T,V]=eigs(Marginal);
Proj=zeros(d,d);
for no=1:n
    Proj=Proj+T(:,no)*ctranspose(T(:,no));
end

%Our non-zero choices are then projected onto this subspace, whilst the
%others are kept as 0.
InitialAssem=zeros(d,d,2*d);
for no=1:n
    no1=firstchoice(no);
    no2=secondchoice(no);
    InitialAssem(:,:,no1)=Proj*ctranspose(PrimeMUB(1,:,no1))*PrimeMUB(1,:,no1)*Proj;
    InitialAssem(:,:,d+no2)=Proj*ctranspose(PrimeMUB(d+1,:,no2))*PrimeMUB(d+1,:,no2)*Proj;
end

%Finally, the assemblage is normalised to have marginal 1.
InitialAssem=2*InitialAssem/(trace(sum(InitialAssem,3)));

