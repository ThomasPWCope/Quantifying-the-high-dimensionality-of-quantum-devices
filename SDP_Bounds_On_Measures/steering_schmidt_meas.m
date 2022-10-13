function [Schmidt_Meas,Unsteer,Steer] = steering_schmidt_meas(InputAssem)
%This function takes a steering assemblage, and gives an upper bound on 
%the steering schmidt measure, by maximising the weight on the unsteerable
%part (greedy algorithm) and assigning maximal dimension to the rest.

%Note: If the assemblage is a qubit assemblage, this code gives the TRUE
%value (not an upper bound).

%InputAssem should be a steering ensemble of the form 
%MeasSet(output,input,:,:) is a dim x dim matrix. (so a 4d array). 
% It returns: 

%Schmidt_Meas: An UPPER BOUND on the steering schmidt measure. 
%Unsteer: a 4d array (out,in,:,:) corresponding to the unsteerable part.
%steer: a 4d array (out,in,:,:) corresponding to the steerable part.

%NOTE: 
%If the assemblage states are qubit states, then this function
%gives the Schmidt measure of the assemblage; otherwise,
%it gives an upper bound on the true value.



[ANum,XNum,Dim1,Dim2]=size(InputAssem);
if Dim1~=Dim2
    error("Not square!")
else
    Dim=Dim1;
end
Unsteer=zeros(ANum,XNum,Dim,Dim);
Steer=zeros(ANum,XNum,Dim,Dim);
cvx_begin sdp quiet
variable LHS(Dim,Dim,ANum^XNum) Hermitian semidefinite %Note the reverse of the other arrays.
variable rhoU(Dim,Dim) Hermitian semidefinite
minimise 1-trace(rhoU)
subject to
IdCon=0;
for xcon=1:XNum
    for acon=1:ANum
    Constraint=0.5*(squeeze(InputAssem(acon,xcon,:,:))+conj(transpose(squeeze(InputAssem(acon,xcon,:,:)))));
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
rhoU==IdCon
cvx_end

if Dim>2
   disp('Only an upper bound!')
end
Schmidt_Meas=log2(Dim)*(1-trace(rhoU));
for xcon= 1:XNum
    for acon= 1:ANum
        Steer(acon,xcon,:,:)=squeeze(InputAssem(acon,xcon,:,:));
        for lcon= 1:ANum^XNum
        Achoices=dec2base(lcon-1,ANum,XNum);
        if Achoices(xcon)==num2str(acon-1)
            Steer(acon,xcon,:,:)=squeeze(Steer(acon,xcon,:,:))-squeeze(LHS(:,:,lcon));
            Unsteer(acon,xcon,:,:)=squeeze(Unsteer(acon,xcon,:,:))+squeeze(LHS(:,:,lcon));
        end
        end
    end
end