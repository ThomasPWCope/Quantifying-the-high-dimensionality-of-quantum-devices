function [dim_meas,JMeas,InComp,LHS] = measurement_dim_meas(MeasSet)
%This function takes a measurement assemblage, and gives an upper bound on 
%the dimension measure, by maximising the weight of the jointly measureable
%part (greedy algorithm) and assigning maximal dimension to the rest.

%Note: If the assemblage is a a set of qubit measurements, this code gives the TRUE
%value (not an upper bound).

%InputAssem should be a measurement ensemble of the form 
%MeasSet(output,input,:,:) is a dim x dim matrix. (so a 4d array). 
% It returns: 
%JMeas: the jointly measureable pseudo-measurement component.
%InComp: the incompatible pseudo-measurement component.
%LHS: Additionally, one can request the individual deterministic
%pseudo-measurements making up the Jointly measurement component.

[ANum,XNum,Dim1,Dim2]=size(MeasSet);
if Dim1~=Dim2
    error("Not square!")
else
    Dim=Dim1;
end
JMeas=zeros(ANum,XNum,Dim,Dim);
InComp=zeros(ANum,XNum,Dim,Dim);
cvx_begin sdp quiet
variable LHS(Dim,Dim,ANum^XNum) Hermitian semidefinite %Note the reverse of the other arrays.
variable alphav
minimise alphav
subject to
IdCon=0;
for xcon=1:XNum
    for acon=1:ANum
    Constraint=0.5*(squeeze(MeasSet(acon,xcon,:,:))+conj(transpose(squeeze(MeasSet(acon,xcon,:,:)))));
    for lcon= 1:ANum^XNum
        Achoices=dec2base(lcon-1,ANum,XNum);
        if Achoices(xcon)==num2str(acon-1)
            Constraint=Constraint-LHS(:,:,lcon);
        end
    end
    Constraint >= 0;
    end
end
for lcon= 1:ANum^XNum
        IdCon=IdCon+LHS(:,:,lcon);
end
(alphav-1)*eye(Dim)+IdCon>=0;
cvx_end

if Dim>2
   disp('Only an upper bound!')
end
dim_meas=log2(Dim)*alphav;
for xcon= 1:XNum
    for acon= 1:ANum
        InComp(acon,xcon,:,:)=squeeze(MeasSet(acon,xcon,:,:));
        for lcon= 1:ANum^XNum
        Achoices=dec2base(lcon-1,ANum,XNum);
        if Achoices(xcon)==num2str(acon-1)
            InComp(acon,xcon,:,:)=squeeze(InComp(acon,xcon,:,:))-squeeze(LHS(:,:,lcon));
            JMeas(acon,xcon,:,:)=squeeze(JMeas(acon,xcon,:,:))+squeeze(LHS(:,:,lcon));
        end
        end
    end
end