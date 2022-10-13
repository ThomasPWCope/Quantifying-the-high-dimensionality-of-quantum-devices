function [dim_meas,BreakState,EntState] = channel_dim_meas(ChannelChoi,varargin)
%This function gives:
%the dimension measure of a quantum channel H_2 -> H_2
%                                           H_2 -> H_3
%                                           H_3 -> H_2
%An upper bound on the dimension measure of a quantum channel otherwise.

%Inputs:
%ChannelChoi: The Choi matrix of the channel (normalised)
% varargin - allows two additional arguments, specifying the dimensions of
% each system (input then output). Default: tries to split the State into equal size spaces.

%Outputs:,,
%dim_meas - the exact dimension measure of the channel
%BreakState: - the entanglement breaking channel component (as a Choi
%matrix)
%EntState: - the entanglement preserving component (as a Choi
%matrix)



%REQUIRES:
%CVX 
%QETLAB.
[R,C]=size(ChannelChoi);
if R~=C
    error("Not square!")
end
if abs(trace(ChannelChoi)-1)>10^-5
    error('Normalise Choi State')
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

disp(join(['Channel is H_',num2str(D1),' to H_',num2str(D2)]))


EBreakingState=0*ChannelChoi;
EPresState=0*ChannelChoi;
cvx_begin sdp quiet
variable alphaV
variable sig0(R,C) Hermitian semidefinite 
minimise alphaV
subject to
1*sig0 >=0
ChannelChoi-sig0>=0
PartialTranspose(sig0,2,D1)>=0
(alphaV-1)*eye(D1)+D1*PartialTrace(sig0,2,D1)>=0
cvx_end
if and(min([D1,D2])==2,max([D1,D2])<=3)
dim_meas=alphaV;
else
    error('Seperable and PPT do not co-incide')
    %dim_meas=log2(min([D1,D2]))*alphaV
end
BreakState=sig0;
EntState=ChannelChoi-sig0;
