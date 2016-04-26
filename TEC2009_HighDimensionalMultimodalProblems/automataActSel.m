% Multi dimentional Action Selection
function swarmTable = automataActSel(p) 
    global action
    [swarmNum dim] = size(p);
    for i =1:dim
        rw = randsrc(100,1,[action; p(:,i)']); % roulette wheel
        act = rw(randi(length(rw)));  % selected action
        swarmTable(act,i) = 1;
        swarmTable(act ~= action,i) = 0;
    end
end