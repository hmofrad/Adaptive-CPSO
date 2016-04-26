%   Randomized CPSO-H
function [gbest_hist fe_hist] = rcpsoh(f,nop,dim,swarmNum)
%     clear,clc
%     f= 9;nop = 10;dim=10;swarmNum=3;
    endgen = 10^4;    % maximum Generation
    if dim == 10, max_fe = 3*10^4; end
    if dim == 30, max_fe = 3*10^5; end

    global  orthm1 orthm2 swarm1 swarm2 K1 K2 K1len K2len s
    if f == 9;  Ub = 32.768; end % Rotated Ackley
    if f == 10; Ub = 600;    end % Rotated Griewank
    if f == 11; Ub = 0.5;    end % Rotated Weierstrass
    if f == 12; Ub = 5.12;   end % Rotated Rastrigin
    if f == 13; Ub = 5.12;   end % Rotated Rastrigin_noncont
    if f == 14; Ub = 500;    end % Rotated Schewfel
    % =================  Rotation =======================
    s = swarmNum;
    K1 = mod(dim,swarmNum); K1len = ceil(dim/swarmNum);    orthm1 = [];
    K2 = swarmNum - K1;     K2len = floor((dim/swarmNum)); orthm2 = [];
    %  define rotation
    K1 = mod(dim,swarmNum); K1len = ceil(dim/swarmNum);    orthm1 = [];
    K2 = swarmNum - K1;     K2len = floor((dim/swarmNum)); orthm2 = [];
    % =========================
    rc=1:dim;swarm1=[]; swarm2=[];
    if K1 == 0
        for i=1:K2
           orthm2 = [orthm2 orthm_generator(K2len)];
        end
    else
        for i=1:K1
            orthm1 = [orthm1 orthm_generator(K1len)];
        end
        for i=1:K1
            swarm1(i,:)=rc((i*K1len-K1len+1):(i*K1len));
        end
        if K2len == 1
            orthm2 = [];
        else
            for i=1:K2
               orthm2 = [orthm2 orthm_generator(K2len)];
            end
        end
    end
    for i=1:K2
        swarm2(i,:)=rc((K1*K1len+i*K2len-K2len+1):(K1*K1len+i*K2len));
    end
    %  randomized swarm members
    rc=randperm(dim);swarm11=[]; swarm22=[];
    if K1 ~= 0
        for i=1:K1
            swarm11(i,:)=rc((i*K1len-K1len+1):(i*K1len));
        end
    end
    for i=1:K2
        swarm22(i,:)=rc((K1*K1len+i*K2len-K2len+1):(K1*K1len+i*K2len));
    end

    Lb = -Ub;    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    c1 = 1.49445; c2 = 1.49445;
    w_max = 0.9;  w_min = 0.4;
    [bst ind]=min(fit_func({f},x));    
    gbest = x(ind,:); gbest(end+1) = bst;
    pbest = x; %initialize Best Particle Position
    gbest_hist = []; fe_hist=[];
    % ========================
    % initilise Q Swarm parameters
    % ========================
    x1 = Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    spd1 = Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    x1(:,end+1)=0;
    for i=1:nop
        x1(i,end) = fit_func({f},x1(i,1:end-1));
    end
    pbest1 = x1; %initialise Best Particle Position
    [bst ind] = min(x1(:,end));
    gbest1 = x1(ind,:); % initialise global best position
    fe =0;i=1;
    while i<=endgen && fe<=max_fe
        w = w_max-(w_max-w_min)*i/endgen;
        if mod(i,2) == 1
            if K1 == 0
                for ii=1:K2
                    for k=1:nop
                        if fit_func({f},b(gbest(1:end-1),swarm22(ii,:),x(k,swarm22(ii,:)))) <...
                                fit_func({f},b(gbest(1:end-1),swarm22(ii,:),pbest(k,swarm22(ii,:))))
                            pbest(k,swarm22(ii,:)) = x(k,swarm22(ii,:));
                        end
                        if fit_func({f},b(gbest(1:end-1),swarm22(ii,:),pbest(k,swarm22(ii,:)))) < gbest(end)
                            gbest(swarm22(ii,:)) = pbest(k,swarm22(ii,:));
                            gbest(end)=fit_func({f},gbest(1:end-1));
                        end
                    end
                    spd(:,swarm22(ii,:)) = w.* spd(:,swarm22(ii,:))...
                        +c1.*rand(nop,K2len).*(pbest(:,swarm22(ii,:))-x(:,swarm22(ii,:)))...
                        +c2.*rand(nop,K2len).*(repmat(gbest(swarm22(ii,:)),nop,1)-x(:,swarm22(ii,:)));
                    spd(:,swarm22(ii,:))=(spd(:,swarm22(ii,:))>Vmax).*Vmax+(spd(:,swarm22(ii,:))<=Vmax).*spd(:,swarm22(ii,:)); 
                    spd(:,swarm22(ii,:))=(spd(:,swarm22(ii,:))<(Vmin)).*(Vmin)+(spd(:,swarm22(ii,:))>=(Vmin)).*spd(:,swarm22(ii,:));
                    x(:,swarm22(ii,:)) = x(:,swarm22(ii,:))+spd(:,swarm22(ii,:)); 
                end
            else
                for ii=1:K1
                    for k=1:nop
                        if fit_func({f},b(gbest(1:end-1),swarm11(ii,:),x(k,swarm11(ii,:)))) <...
                                fit_func({f},b(gbest(1:end-1),swarm11(ii,:),pbest(k,swarm11(ii,:))))
                            pbest(k,swarm11(ii,:)) = x(k,swarm11(ii,:));
                        end
                        if fit_func({f},b(gbest(1:end-1),swarm11(ii,:),pbest(k,swarm11(ii,:)))) < gbest(end)
                            gbest(swarm11(ii,:)) = pbest(k,swarm11(ii,:));
                            gbest(end)=fit_func({f},gbest(1:end-1));
                        end
                    end
                    spd(:,swarm11(ii,:)) = w.* spd(:,swarm11(ii,:))...
                        +c1.*rand(nop,K1len).*(pbest(:,swarm11(ii,:))-x(:,swarm11(ii,:)))...
                        +c2.*rand(nop,K1len).*(repmat(gbest(swarm11(ii,:)),nop,1)-x(:,swarm11(ii,:)));
                    spd(:,swarm11(ii,:))=(spd(:,swarm11(ii,:))>Vmax).*Vmax+(spd(:,swarm11(ii,:))<=Vmax).*spd(:,swarm11(ii,:)); 
                    spd(:,swarm11(ii,:))=(spd(:,swarm11(ii,:))<(Vmin)).*(Vmin)+(spd(:,swarm11(ii,:))>=(Vmin)).*spd(:,swarm11(ii,:));
                    x(:,swarm11(ii,:)) = x(:,swarm11(ii,:))+spd(:,swarm11(ii,:)); 
                end
                for ii=1:K2
                    for k=1:nop
                        if fit_func({f},b(gbest(1:end-1),swarm22(ii,:),x(k,swarm22(ii,:)))) <...
                                fit_func({f},b(gbest(1:end-1),swarm22(ii,:),pbest(k,swarm22(ii,:))))
                            pbest(k,swarm22(ii,:)) = x(k,swarm22(ii,:));
                        end
                        if fit_func({f},b(gbest(1:end-1),swarm22(ii,:),pbest(k,swarm22(ii,:)))) < gbest(end)
                            gbest(swarm22(ii,:)) = pbest(k,swarm22(ii,:));
                            gbest(end)=fit_func({f},gbest(1:end-1));
                        end
                    end
                    spd(:,swarm22(ii,:)) = w.* spd(:,swarm22(ii,:))...
                        +c1.*rand(nop,K2len).*(pbest(:,swarm22(ii,:))-x(:,swarm22(ii,:)))...
                        +c2.*rand(nop,K2len).*(repmat(gbest(swarm22(ii,:)),nop,1)-x(:,swarm22(ii,:)));
                    spd(:,swarm22(ii,:))=(spd(:,swarm22(ii,:))>Vmax).*Vmax+(spd(:,swarm22(ii,:))<=Vmax).*spd(:,swarm22(ii,:)); 
                    spd(:,swarm22(ii,:))=(spd(:,swarm22(ii,:))<(Vmin)).*(Vmin)+(spd(:,swarm22(ii,:))>=(Vmin)).*spd(:,swarm22(ii,:));
                    x(:,swarm22(ii,:)) = x(:,swarm22(ii,:))+spd(:,swarm22(ii,:)); 
                end
            end
            fe = fe + swarmNum*nop;
            rc=randperm(round(nop/2));   ind=rc(1);
            if x1(ind,end) == gbest1(end)
                ind=rc(2);
            end
            x1(ind,:) =  gbest;
        else
            % Q SWARM
            for j=1:nop
                if (x1(j,end) < pbest1(j,end))
                    pbest1(j,:) = x1(j,:);
                end
                if (pbest1(j,end) < gbest1(end))
                    gbest1 = pbest1(j,:);
                end
            end
            spd1 = w.* spd1+c1.*rand(nop,dim).*(pbest1(:,1:end-1)-x1(:,1:end-1))...
                               +c2.*rand(nop,dim).*(repmat(gbest1(1:end-1),nop,1) - x1(:,1:end-1));
            spd1=(spd1>Vmax).*Vmax+(spd1<=Vmax).*spd1; 
            spd1=(spd1<(Vmin)).*(Vmin)+(spd1>=(Vmin)).*spd1;
            x1(:,1:end-1) = x1(:,1:end-1)+spd1;
            x1(:,end) = fit_func({f},x1(:,1:end-1));
            fe = fe+nop;
        % ==== Interval THREE ====
            if K1 == 0
                for ii=1:K2
                    rc=randperm(round(nop/2));   ind=rc(1);
                    if x(ind,swarm22(ii,:)) == gbest(swarm22(ii,:))
                        ind=rc(2);
                    end
                    x(ind,swarm22(ii,:)) = gbest1(swarm22(ii,:));
                end
            else
                for ii=1:K1
                    rc=randperm(round(nop/2));   ind=rc(1);
                    if x(ind,swarm11(ii,:)) == gbest(swarm11(ii,:))
                        ind=rc(2);
                    end
                    x(ind,swarm11(ii,:)) = gbest1(swarm11(ii,:));
                end
                for ii=1:K2
                    rc=randperm(round(nop/2));   ind=rc(1);
                    if x(ind,swarm22(ii,:)) == gbest(swarm22(ii,:))
                        ind=rc(2);
                    end
                    x(ind,swarm22(ii,:)) = gbest1(swarm22(ii,:));
                end
            end
        end
        gbest_hist = [gbest_hist gbest(end)]; fe_hist = [fe_hist fe];
        if mod(i,500)==0,fprintf('fun=%u,Gene=%u,Fit_Eval=%u,Gbest=%e\n',f,i,fe,gbest(end)),end
        if ((i>1300 && gbest_hist(length(gbest_hist)-500) == gbest(end)) || gbest(end) ==0),break,end
        i=i+1;
    end
end