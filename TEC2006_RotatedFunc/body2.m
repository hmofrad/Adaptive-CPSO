clear,clc
%CLPSO
swarmNum = [3 6];nop = 40;dim=30;
for i=1:2
    for j=9:14
        for k=1:25
            if i == 1
                [gh fh]=icpsoh(j,nop,dim,swarmNum(i));
                G(j-8,k) =gh(end);
                F(j-8,k) =fh(end);
            else
                [gh fh]=icpsoh(j,nop,dim,swarmNum(i));
                G(j-8+6,k) =gh(end);
                F(j-8+6,k) =fh(end);
            end
            i,j,k
        end
    end
end
save('icpsohR30D.mat')
