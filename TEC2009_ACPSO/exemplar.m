function [pbestf] = exemplar(x,pbest,j)
    % j as particle index
    % Learning Probability PCi
    [nop dim]= size(x); dim = dim - 1;
%     pbestf = zeros(1,dim);
pbestf = [];
    di = 1:nop; % discretized particle index
    Pc = 0.05+0.45*((exp(10*(di-1)/(nop-1))-1)/(exp(10)-1));
    % ======================================
    % calculatind exepmlar Pbest Function
    for k=1:dim
        if rand>Pc(j)
            pbestf(k) = pbest(j,k);
        else
            pind = di; %insert raw particles index in particle index vector
            pind = pind(pind ~= j);
            f1ind = pind(randi([1 length(pind)]));
            pind = pind(pind ~= f1ind);
            f2ind = pind(randi([1 length(pind)]));
            if pbest(f1ind,end)<pbest(f2ind,end)
                pbestf(k) = x(f1ind,k);
            else
                pbestf(k) = x(f2ind,k);
            end
        end
    end
    if pbestf == pbest(j,1:end-1)
        Dind = randi([1 dim]);
        Pind = randi([1 nop]);
        pbestf(Dind) = pbest(Pind,Dind);
    end
end