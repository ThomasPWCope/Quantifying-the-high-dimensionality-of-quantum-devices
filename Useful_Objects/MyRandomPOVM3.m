function [POVM] = MyRandomPOVM3(outcomes)
%This function creates a random qubit POVM; it is completely heuristic, in
%that it does not follow any well justified distribution. It is, however,
%useful, it that sets of POVMs created with this tend to highlight
%dimension measure properties of POVMs.

%Inputs: 
%outcomes - the number of desired outcomes for the POVM.
%Outputs:
%POVM: a (2,2,outcomes) array giving a single qubit POVM.

POVM=zeros(2,2,outcomes);
for i=1:outcomes
    randmultiplier=0.8+0.4*rand();
    Bloch=rand(1,3);
    Bloch=Bloch/norm(Bloch);
    Part1=randmultiplier*(eye(2)+Bloch(1)*Pauli(1)+Bloch(2)*Pauli(2)+Bloch(3)*Pauli(3));
    randmultiplier=0.8+0.4*rand();
    Bloch=rand(1,3);
    Bloch=Bloch/norm(Bloch);
    Part2=randmultiplier*(eye(2)+Bloch(1)*Pauli(1)+Bloch(2)*Pauli(2)+Bloch(3)*Pauli(3));
    POVM(:,:,i)=0.99*Part1+0.01*Part2;
end
POVMSum=sum(POVM,3);
for i=1:outcomes
    POVM(:,:,i)=(POVMSum^(-1/2))*POVM(:,:,i)*(POVMSum^(-1/2));
    POVM(:,:,i)=1/2*(ctranspose(POVM(:,:,i))+POVM(:,:,i));
end

    