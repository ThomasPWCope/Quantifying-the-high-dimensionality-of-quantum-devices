function [dim_meas,JMeas,InComp] = measurement_incomp_weight(MeasSet,varargin)
%This function takes a set of measurements and decomposes them into a
%convex combination of jointly compatible and incompatible measurements.
%If chosen (or left as default), this is then weighted by the logarithm of
%the dimension, for comparision with the dimension measure. As stated, this
%quantity is NOT equal to the dimension measure.

%Input:
%MeasSet:
%Measurement Set should be a measurement ensemble of the form 
%MeasSet(output,input,:,:) is a dim x dim matrix. (so a 4d array)
%varargin
if len(varargin)==1
    if or(0,1)
    weightingyn=varargin;    
        
    else
       error('Only 0 (no log_2 d) or 1 (log_2 d) Possible')
    end
    elseif len(varargin)>1
        error('Too Many Arguments!')
else
        weightingyn=1
end        

%Outputs:
%dim_meas - if qubit measurements, the exact convex weight. Otherwise, an
%upper bound.
%JMeas - the jointly compatible measurement set component.
%InComp - the incompatible measurement set component.


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
(alphav-1)*eye(Dim)+IdCon==0;
cvx_end

if Dim>2
   disp('Only an upper bound!')
end

if weightingyn==0
    dim_meas=alphav;
else
    dim_meas=alphav*log2(Dim);
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