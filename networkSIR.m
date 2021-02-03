function [SVec, IVec, RVec, tracking] = networkSIR(n, a, b, infectivePeriod, A, tf, I0)

hwait = waitbar(0,'Please wait. Generating Networked Model');

%set population size for the networks
popSize = n;

R = 0;              %initial recovered population
S = popSize - I0;    %susceptible population
I = I0;

IVec = [I];
RVec = [R];
SVec = [S];

%tracking who is S, I and R
%  in row 2, 0=susceptible, 1=infected, 2=recovered
tracking = [1:popSize];
tracking(2,:) = zeros;
tracking(3,:) = Inf;
tracking(7,:) = zeros;

%randomly assign initially infected
for i = 1:I

    %random number generator
    r = randi(popSize, 1);
    patient0 = r;

    %setting infection status and time of infection
    tracking(2,r) = 1;
    tracking(3,r) = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 2:tf
    waitbar(t/tf,hwait,sprintf('Please wait. Generating networked model\n%.1f %%',(t/tf)*100));
    S = 0;
    I = 0;
    R = 0;
    for p = tracking(1,:)
        neighbours = neighbors(graph(A),p)';
    
        for i = 1:b
            if not(isempty(neighbours))
                ind = randperm(length(neighbours), 1);
                visits(i) = neighbours(ind);                
                neighbours(neighbours==visits(i)) = [];
            end
        end
        if tracking(2,p) == 1 && (t - tracking(3,p)) > infectivePeriod
               tracking(2,p) = 2; 
               R = R+1;
               I = I-1;
        end   
        for visit = visits
            if tracking(2,p) == 0 && tracking(2, visit) == 1
               r = rand;
               if r <= a
                   tracking(2,p) = 1;
                   tracking(3,p) = t;
                   I = I+1;
                   S = S-1;
                   tracking(7,visit) = tracking(7,visit) + 1;
               end
            end
            
            if tracking(2,p) == 1 && tracking(2, visit) == 0
               r = rand;
               if r <= a
                  tracking(2, visit) = 1;
                  tracking(3,visit) = t;
                  I = I+1;
                  S = S-1;
                  tracking(7,p) = tracking(7,p) + 1;
               end
            end 
            
            
        end
         
    end
    
    IVec(t) = IVec(t-1)+I;
    RVec(t) = RVec(t-1)+R;
    SVec(t) = SVec(t-1)+S;
    
    
end

for n = 1:popSize
    tracking(4,n) = local_clustering_coefficient(A,n);
end
for n = 1:popSize
    tracking(5,n) = degree(graph(A),n);
end
for n = 1:popSize
    [TR,D] = shortestpathtree(graph(A), n, patient0);
    tracking(6,n) = D;
end


delete(hwait);
end
