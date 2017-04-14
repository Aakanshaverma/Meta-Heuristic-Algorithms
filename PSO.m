%PSO Algorithm



function [MinCost] = PSO(ProblemFunction, DisplayFlag)
if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end

[OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction);

OPTIONS.Keep = 2; 
OPTIONS.neighbors = 0; 
OPTIONS.w = 0.3; 
OPTIONS.c1 = 1; 
OPTIONS.c2 = 1; 
OPTIONS.c3 = 1; 

vel = zeros(OPTIONS.popsize, OPTIONS.numVar); 
pbest = Population; 
nbest = Population; 
gbest = Population(1); 


for GenIndex = 1 : OPTIONS.Maxgen
    if ~OPTIONS.OrderDependent
        
        for i = 1 : OPTIONS.popsize
            [chrom, indices] = sort(Population(i).chrom);
            Population(i).chrom = chrom;
            VelTemp = vel(i, :);
            for j = 1 : OPTIONS.numVar
                vel(i, j) = VelTemp(indices(j));
            end
        end
    end
    
    if Population(1).cost < gbest.cost
        gbest = Population(1);
    end
    
    for i = 1 : OPTIONS.popsize 
        
        if Population(i).cost < pbest(i).cost
            pbest(i) = Population(i);
        end
        
        Distance = zeros(OPTIONS.popsize, 1);
        for j = 1 : OPTIONS.popsize 
            Distance(j) = norm(Population(i).chrom-Population(j).chrom);
        end
        [Distance, indices] = sort(Distance);
        nbest(i).cost = inf;
        for j = 2 : OPTIONS.neighbors+1
            nindex = indices(j);
            if Population(nindex).cost < nbest(i).cost
                nbest(i) = Population(nindex);
            end
        end
    end
    
    for i = OPTIONS.Keep+1 : OPTIONS.popsize
        r = rand(3, OPTIONS.numVar);
        x = Population(i).chrom;
        deltaVpersonal = OPTIONS.c1 * r(1,:) .* (pbest(i).chrom - x);
        deltaVswarm = OPTIONS.c2 * r(2,:) .* (gbest.chrom - x);
        deltaVneighborhood = OPTIONS.c3 * r(3,:) .* (nbest(i).chrom - x);
        vel(i,:) = OPTIONS.w * vel(i,:) + deltaVpersonal + deltaVswarm + deltaVneighborhood;
        Population(i).chrom = x + vel(i,:);
    end 
     
    Population = ClearDups(Population, MaxParValue, MinParValue);
    
    Population = FeasibleFunction(OPTIONS, Population);
    
    Population = CostFunction(OPTIONS, Population);
    
    Population = PopSort(Population);
    
    [AverageCost, nLegal] = ComputeAveCost(Population);
    
    MinCost = [MinCost Population(1).cost];
    AvgCost = [AvgCost AverageCost];
    if DisplayFlag
        disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ',...
            num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
    end
end
Conclude(DisplayFlag, OPTIONS, Population, nLegal, MinCost);

return;
