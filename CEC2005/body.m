clear,clc
alpha = 0.1;beta = 0.1;   %LRP
swarmNum = [3 6];nop = 30;dim=30;
for i=1:2
    for j=15:25
        for k=1:25
            if i == 1
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-14,k) =gh(end);
                F(j-14,k) =fh(end);
            else
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-14+11,k) =gh(end);
                F(j-14+11,k) =fh(end);
            end
            swarmNum(i),j,k
        end
    end
end
save('ACPSOCF3RP.mat')
%   ==================================
clear,clc
alpha = 0.1;beta = 0.01;   %LReP
swarmNum = [3 6];nop = 30;dim=30;
for i=1:2
    for j=15:25
        for k=1:25
            if i == 1
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-14,k) =gh(end);
                F(j-14,k) =fh(end);
            else
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-14+11,k) =gh(end);
                F(j-14+11,k) =fh(end);
            end
            swarmNum(i),j,k,alpha
        end
    end
end
save('ACPSOCF3ReP.mat')
%   ==================================
clear,clc
alpha = 0.1;beta = 0.0;   %LRI
swarmNum = [3 6];nop = 30;dim=30;
for i=1:2
    for j=15:25
        for k=1:25
            if i == 1
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-14,k) =gh(end);
                F(j-14,k) =fh(end);
            else
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-14+11,k) =gh(end);
                F(j-14+11,k) =fh(end);
            end
            swarmNum(i),j,k,beta
        end
    end
end
save('ACPSOCF3RI.mat')