function [Schmidt_Meas,SepState,EntState] = state_schmidt_meas(State,varargin)
%This function gives:
%the Schmidt measure of a quantum state on  H_2 x H_2
%                                           H_2 x H_3
%                                           H_3 x H_2


%Inputs:
%State - a density operator (:,:) matrix
% varargin - allows two additional arguments, specifying the dimensions of
% each system. Default: tries to split the State into equal size spaces.

%Outputs:
%Schmidt_Meas - the exact Schmidt measure for states on the above spaces.
%SepState: - the seperable component
%EntState: - the entangled compontent




%REQUIRES:
%CVX 
%QETLAB.
[R,C]=size(State);
if R~=C
    error("Not square!")
end
if abs(trace(State)-1)>10^-5
    error('Normalise State')
end
if length(varargin)>2
    error('too many inputs')
elseif length(varargin)==2
    D1=varargin{1};
    D2=varargin{2};
    if D1*D2~=R
        error('Wrong Dimensions')
    end
elseif length(varargin)==1
    D1=varargin{1};
    D2=fix(R/D1);
    if D1*D2~=R
        error('Wrong Dimensions')
    end
else
if fix(sqrt(R))^2==R
    D1=fix(sqrt(R));
    D2=fix(sqrt(R));
else
    error('Specify input/output dimensions')
end
end

disp(join(['State belongs to H_',num2str(D1),' X H_',num2str(D2)]))



cvx_begin sdp quiet
variable alphaV
variable sig0(R,C) Hermitian semidefinite 
minimise 1-trace(sig0)
subject to
1*sig0 >=0
State-sig0>=0
PartialTranspose(sig0,2,D1)>=0
cvx_end
if and(min([D1,D2])==2,max([D1,D2])<=3)
Schmidt_Meas=cvx_optval;
else
    error('Seperable and PPT not equivalent')
    %Schmidt_Meas=log2(min([D1,D2]))*cvx_optval
end
SepState=sig0;
EntState=State-sig0;
