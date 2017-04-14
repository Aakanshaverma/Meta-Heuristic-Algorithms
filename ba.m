%ABC Algorithm
%Plots a graph iterations vs. best solution

clc;
clear;
close all;



CostFunction=@(x) Sphere(x);       
nVar=5;            
VarSize=[1 nVar];  

VarMin=-10;        
VarMax= 10;         



MaxIt=1000;          
nScoutBee=30;                           

nSelectedSite=round(0.5*nScoutBee);     

nEliteSite=round(0.4*nSelectedSite);   
nSelectedSiteBee=round(0.5*nScoutBee);  
nEliteSiteBee=2*nSelectedSiteBee;       

r=0.1*(VarMax-VarMin);	

rdamp=0.95;            


empty_bee.Position=[];
empty_bee.Cost=[];


bee=repmat(empty_bee,nScoutBee,1);

for i=1:nScoutBee
    bee(i).Position=unifrnd(VarMin,VarMax,VarSize);
    bee(i).Cost=CostFunction(bee(i).Position);
end


[~, SortOrder]=sort([bee.Cost]);
bee=bee(SortOrder);

BestSol=bee(1);


BestCost=zeros(MaxIt,1);



for it=1:MaxIt
   
    for i=1:nEliteSite
        
        bestnewbee.Cost=inf;
        
        for j=1:nEliteSiteBee
            newbee.Position=PerformBeeDance(bee(i).Position,r);
            newbee.Cost=CostFunction(newbee.Position);
            if newbee.Cost<bestnewbee.Cost
                bestnewbee=newbee;
            end
        end

        if bestnewbee.Cost<bee(i).Cost
            bee(i)=bestnewbee;
        end
        
    end
    
    for i=nEliteSite+1:nSelectedSite
        
        bestnewbee.Cost=inf;
        
        for j=1:nSelectedSiteBee
            newbee.Position=PerformBeeDance(bee(i).Position,r);
            newbee.Cost=CostFunction(newbee.Position);
            if newbee.Cost<bestnewbee.Cost
                bestnewbee=newbee;
            end
        end

        if bestnewbee.Cost<bee(i).Cost
            bee(i)=bestnewbee;
        end
        
    end
    

    for i=nSelectedSite+1:nScoutBee
        bee(i).Position=unifrnd(VarMin,VarMax,VarSize);
        bee(i).Cost=CostFunction(bee(i).Position);
    end
    
    [~, SortOrder]=sort([bee.Cost]);
    bee=bee(SortOrder);
    
    
    BestSol=bee(1);
    
  
    BestCost(it)=BestSol.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    
    r=r*rdamp;
    
end



figure;
plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');

