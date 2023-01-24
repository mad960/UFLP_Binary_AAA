    function binaryVector=CreateVectorByTransferFunction(continuousVector,transFuncIndex)
        
        
        if transFuncIndex==1
            binaryVector=mod(round(mod(continuousVector,2)),2);
            return;
        end
        
        if transFuncIndex==2
            binaryVector= mod( floor(continuousVector),2);
            return;
        end
        
        binaryVector=zeros(size(continuousVector));
        
        for j=1: numel(continuousVector)
            
            x=continuousVector(j);
            
            if transFuncIndex==3
                s=1/(1+exp(-2*x)); %S1 transfer function
                
            else  if transFuncIndex==4
                    s=1/(1+exp(-x));   %S2 transfer function
                    
                else  if transFuncIndex==5
                        s=1/(1+exp(-x/2)); %S3 transfer function
                        
                    else  if transFuncIndex==6
                            s=1/(1+exp(-x/3));  %S4 transfer function
                            
                        else if transFuncIndex==7
                                s=abs(erf(((sqrt(pi)/2)*x))); %V1 transfer function
                                
                            else  if transFuncIndex==8
                                    s=abs(tanh(x)); %V2 transfer function
                                    
                                else  if transFuncIndex==9
                                        s=abs(x/sqrt((1+x^2))); %V3 transfer function
                                        
                                    else if transFuncIndex==10
                                            s=abs((2/pi)*atan((pi/2)*x)); %V4 transfer function (VPSO)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            if transFuncIndex<=6 %S-shaped transfer functions
                if rand<s % Equation (4) and (8)
                    binaryVector(1,j)=1;
                else
                    binaryVector(1,j)=0;
                end
                
            else if transFuncIndex>6 && transFuncIndex<=10 %V-shaped transfer functions
                    if rand<s %Equation (10)
                        binaryVector(1,j)=~(mod(x,1)>0.5);
                    else
                        binaryVector(1,j)=(mod(x,1)>0.5);
                    end
                end
            end
        end
        
    end