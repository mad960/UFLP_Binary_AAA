clearvars; 
clc;

Delta = 1.5;
Ap = 0.5;
e = 0.3;

N=40;
transFuncIndex=7; % 1-10  
maxFEs = 80000;
MaxHesaplama=round(maxFEs/N); % Maximum numbef of iterations
load('Cap_Problems.mat');
probNo = size(Problems, 1);

for capNo=1:15
capName = Problems{capNo,1};
bestFitnessValue = Problems{capNo, 2};
matrix = Problems{capNo,3};
initializationCost = Problems{capNo,4};

LB = -8;
UB = 8;
D = size(matrix,1);
 
runtime = 5;
    for run=1: runtime
        [ObjMin, BestColony] = AAA(MaxHesaplama, LB, UB, N, D, Delta, Ap, e, matrix, initializationCost,transFuncIndex);
        BestBinaryColony = CreateVectorByTransferFunction(BestColony, transFuncIndex);
        fprintf('\nBest Binary Colony: %s\n',int2str(BestBinaryColony));
    end
end

