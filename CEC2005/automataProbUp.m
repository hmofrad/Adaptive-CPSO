function p = automataProbUp(p,imp_tag,swarmInd,swarmDim,swarmTable,alpha,beta)
    global action
    [swarmNum dim] = size(swarmTable);
    r = length(action);
    swarm = 1:swarmNum;
    swarmNind = swarm ~=swarmInd;
% nop = 30;
%     if (imp_tag >= nop/2)
    if (imp_tag==1)
       p(swarmInd,swarmDim) = p(swarmInd,swarmDim) + alpha.*(1 - p(swarmInd,swarmDim));   % desired action
       p(swarmNind,swarmDim) = (1-alpha).*p(swarmNind,swarmDim);
    else
        p(swarmInd,swarmDim) = (1 - beta).*p(swarmInd,swarmDim);   % non-desired action
        p(swarmNind,swarmDim) = (beta/(r-1))+(1-beta).*p(swarmNind,swarmDim);
    end
end