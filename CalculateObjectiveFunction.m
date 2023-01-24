function Costs = CalculateObjectiveFunction(Pop, matrix, initializationCost,transFuncIndex)

    binPop=zeros(size(Pop,1),size(Pop,2));
    Costs=zeros(size(Pop,1),1);
    numberOfFacility=size(Pop,2); %%% also dimension of the problem
	numberOfCustomer=size(matrix,2);

for j=1: size(Pop,1)
    binaryVector=CreateVectorByTransferFunction(Pop(j,:),transFuncIndex);
    if sum(binaryVector)==0
        binaryVector(randi(numberOfFacility))=1;
    end
    binPop(j,:)=binaryVector;
	totalCost=0;
	openIndices= binaryVector==1;

	for i=1: numberOfFacility   
		totalCost= totalCost+ binaryVector(i)*initializationCost(i);
	end

	for i=1: numberOfCustomer %
		minVal= min(matrix(openIndices,i)); %% min degere sahip müþteri id si al.
		totalCost= totalCost+  minVal;
    end
    Costs(j)=totalCost;
end

end

