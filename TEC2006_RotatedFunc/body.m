clear,clc
alpha = 0.1;beta = 0.1;   %LRP
swarmNum = [3 6];nop = 10;dim=10;
for i=1:2
    for j=9:14
        for k=1:25
            if i == 1
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-8,k) =gh(end);
                F(j-8,k) =fh(end);
            else
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-8+6,k) =gh(end);
                F(j-8+6,k) =fh(end);
            end
            swarmNum(i),i,j
        end
    end
end
save('ACPSORF1RP.mat')
% =====================
clear,clc
alpha = 0.1;beta = 0.01;   %LReP
swarmNum = [3 6];nop = 10;dim=10;
for i=1:2
    for j=9:14
        for k=1:25
            if i == 1
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-8,k) =gh(end);
                F(j-8,k) =fh(end);
            else
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-8+6,k) =gh(end);
                F(j-8+6,k) =fh(end);
            end
            swarmNum(i),i,j,alpha
        end
    end
end
save('ACPSORF1ReP.mat')
% =====================
clear,clc
alpha = 0.1;beta = 0.0;   %LRI
swarmNum = [3 6];nop = 10;dim=10;
for i=1:2
    for j=9:14
        for k=1:25
            if i == 1
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-8,k) =gh(end);
                F(j-8,k) =fh(end);
            else
                [gh fh]=ACPSO(j,nop,dim,swarmNum(i),alpha,beta);
                G(j-8+6,k) =gh(end);
                F(j-8+6,k) =fh(end);
            end
            swarmNum(i),i,j,beta
        end
    end
end
save('ACPSORF1RI.mat')
% =====================