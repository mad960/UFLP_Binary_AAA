function [ObjMin, BestColony] = AAA(MaxHesaplama, LB, UB, N, D, Delta, Ap, e, matrix, initializationCost,transFuncIndex)

    Colonies = rand(N,D)*(UB-LB)+LB;
    ObjectiveFun = ones(N,1);
    for i=1:N
        ObjectiveFun(i) = CalculateObjectiveFunction(Colonies(i,:), matrix, initializationCost,transFuncIndex); 
    end
    
    G = ones(N,1);
    G = Greatness(G, ObjectiveFun);
    
    StarvationV = zeros (N,1);

    t=0;
    while (t<MaxHesaplama)

        FrictionSurfaces = 2*pi*(3*G/4*pi).^(2/3);
        Energy = NormalizeVector(G);
        
        for i=1: N
            StarvationMarker = 1;
            currentObjFunValue = ObjectiveFun(i);
            
            while(Energy(i)>0)
                j = find(ObjectiveFun == min(ObjectiveFun));          
                RandCells = randperm(D,D);               
                if (D==1)
                    currentCells = [Colonies(i,RandCells(1))];
                    Colonies(i,RandCells(1)) = Colonies(i,RandCells(1)) + ( Colonies(j(1),RandCells(1)) -  Colonies(i,RandCells(1))) * (Delta - FrictionSurfaces(i)) * (2*rand(1)-1);
                elseif (D==2)
                    currentCells = [Colonies(i,RandCells(1)), Colonies(i,RandCells(2))];
                    Colonies(i,RandCells(1)) = Colonies(i,RandCells(1)) + ( Colonies(j(1),RandCells(1)) -  Colonies(i,RandCells(1))) * (Delta - FrictionSurfaces(i)) * (2*rand(1)-1);
                    Colonies(i,RandCells(2)) = Colonies(i,RandCells(2)) + ( Colonies(j(1),RandCells(2)) -  Colonies(i,RandCells(2))) * (Delta - FrictionSurfaces(i)) * (2*rand(1)-1);
                else
                    currentCells = [Colonies(i,RandCells(1)), Colonies(i,RandCells(2)), Colonies(i,RandCells(3))];
                    Colonies(i,RandCells(1)) = Colonies(i,RandCells(1)) + ( Colonies(j(1),RandCells(1)) -  Colonies(i,RandCells(1))) * (Delta - FrictionSurfaces(i)) * (2*rand(1)-1);
                    Colonies(i,RandCells(2)) = Colonies(i,RandCells(2)) + ( Colonies(j(1),RandCells(2)) -  Colonies(i,RandCells(2))) * (Delta - FrictionSurfaces(i)) * (2*rand(1)-1);
                    Colonies(i,RandCells(3)) = Colonies(i,RandCells(3)) + ( Colonies(j(1),RandCells(3)) -  Colonies(i,RandCells(3))) * (Delta - FrictionSurfaces(i)) * (2*rand(1)-1);
                end
                
                ObjectiveFun(i) = CalculateObjectiveFunction(Colonies(i,:), matrix, initializationCost,transFuncIndex);                
                t = t+1;
                Energy(i) = Energy(i) - (e/2);
                
                if(ObjectiveFun(i)<currentObjFunValue)
                    StarvationMarker = 0;
                else
                    if (D == 1)
                        Colonies(i,RandCells(1)) = currentCells(1);
                    elseif (D == 2)
                        Colonies(i,RandCells(1)) = currentCells(1);
                        Colonies(i,RandCells(2)) = currentCells(2);
                    else
                        Colonies(i,RandCells(1)) = currentCells(1);
                        Colonies(i,RandCells(2)) = currentCells(2);
                        Colonies(i,RandCells(3)) = currentCells(3);
                    end                    
                    ObjectiveFun(i) = currentObjFunValue;
                    Energy(i) = Energy(i) - (e/2);
                     
                end  
                fprintf('Time = %d     Iteration = %d     ObjVal = %d     Colony = %d     ObjMin = %d     Best Colony = %d\n', t, i, ObjectiveFun(i) , find(ObjectiveFun==ObjectiveFun(i)), min(ObjectiveFun), find(ObjectiveFun==min(ObjectiveFun)));
            end
            if(StarvationMarker == 1)
                StarvationV(i) = StarvationV(i) + 1;
            end           
        end
        G = Greatness(G, ObjectiveFun);
        cell = randperm(D,1);
        indexMax = G ==max(G);
        indexMin = G ==min(G);
        Colonies(indexMin,cell) = Colonies(indexMax,cell);
        if(rand < Ap)
            StarvationV(i) = max(StarvationV) + (max(G) - max(StarvationV))*rand;
        end                              
    end
    fprintf(' \n');
    ObjMin = min(ObjectiveFun);
    colonyPosition= find(ObjectiveFun==ObjMin);
    fprintf('\nTime = %d  ObjMin = %d  Colony = %d\n',t, ObjMin , colonyPosition);
    
    for c=1:D
        BestColony(c) = Colonies(colonyPosition(1),c);
    end
    
end

%%Buyukluk Fonksiyonu
function GV = Greatness(GV, ObjFun)
    ObjFun = 1 - NormalizeVector(ObjFun);
        for i=1:length(GV)
            Ks = abs(GV(i)/2);
            U = ObjFun(i) / (Ks + ObjFun(i));
            GV(i) = GV(i) + U * GV(i);
        end
end

%%Normalize Fonksiyonu
function normalizedV = NormalizeVector (V)
    temp = ones(size(V));
        for i=1:length(V)
            if (max(V)-min(V) == 0)
                temp(i)=V(i)-min(V);
            end
            if (max(V)-min(V) ~= 0)
                temp(i)=(V(i)-min(V))/(max(V)-min(V));
            end
        end
    normalizedV = temp;
end

