clear,clc
alpha = 0.1;beta = 0.0;   %LRP
swarmNum = 3;nop = 10;dim=10;
for i=1:1
    for j=1:10
        [gh fh]=ACPSO(i,nop,dim,swarmNum,alpha,beta);
        G3(i,j) =gh(end);
        F3(i,j) =fh(end);
        3,i,j
    end
end
swarmNum = 6;
for i=1:1
    for j=1:10
        [gh fh]=ACPSO(i,nop,dim,swarmNum,alpha,beta);
        G6(i,j) =gh(end);
        F6(i,j) =fh(end);
        6,i,j
    end
end
% save('F1LRP.mat')
% % =====================
% clear,clc
% alpha = 0.1;beta = 0.01;   %LReP
% swarmNum = 3;nop = 10;dim=10;
% for i=1:1
%     for j=1:25
%         [gh fh]=ACPSO(i,nop,dim,swarmNum,alpha,beta);
%         G3(i,j) =gh(end);
%         F3(i,j) =fh(end);
%         3,i,j,alpha
%     end
% end
% swarmNum = 6;
% for i=1:1
%     for j=1:25
%         [gh fh]=ACPSO(i,nop,dim,swarmNum,alpha,beta);
%         G6(i,j) =gh(end);
%         F6(i,j) =fh(end);
%         6,i,j,alpha
%     end
% end
% save('F1LReP.mat')
% % =====================
% clear,clc
% alpha = 0.1;beta = 0.0;   %LRI
% swarmNum = 3;nop = 40;dim=30;
% for i=1:1
%     for j=1:25
%         [gh fh]=ACPSO(i,nop,dim,swarmNum,alpha,beta);
%         G3(i,j) =gh(end);
%         F3(i,j) =fh(end);
%         3,i,j,beta
%     end
% end
% swarmNum = 6;
% for i=1:1
%     for j=1:25
%         [gh fh]=ACPSO(i,nop,dim,swarmNum,alpha,beta);
%         G6(i,j) =gh(end);
%         F6(i,j) =fh(end);
%         6,i,j,beta
%     end
% end
% save('F1LRI.mat')