% function [x_std  gbest_fit gbest_hist] = pso(f,bnd,dim,nop,endgen)
clear,clc
f= 12;nop = 30;dim=30;swarmNum=6;alpha=0.1;beta=0.1;
%     [gh fe]=ACPSO(15,30,30,3,0.01,0.01);
    endgen = 10^4;    % maximum Generation
    max_fe = 10^6; % current fitness evaluation
    global initial_flag action 
    initial_flag = 0;
    if f == 1; Bounds=[-100,100]; f_bias=-450;end
    if f == 2; Bounds=[-100,100]; f_bias=-450;end
    if f == 3; Bounds=[-100,100]; f_bias=-450;end
    if f == 4; Bounds=[-100,100]; f_bias=-450;end
    if f == 5; Bounds=[-100,100]; f_bias=-310;end
    if f == 6; Bounds=[-100,100]; f_bias=+390;end
    if f == 7; Bounds=[0,600];    f_bias=-180;end
    if f == 8; Bounds=[-32, 32];  f_bias=-140;end
    if f == 9; Bounds=[-5,5];     f_bias=-330;end
    if f == 10;Bounds=[-5,5];     f_bias=-330;end
    if f == 11;Bounds=[-0.5,0.5]; f_bias=+90;end
    if f == 12;Bounds=[-100,100]; f_bias=-460;end
    if f == 13,Bounds=[-3,1];     f_bias=-130;end
    if f == 14;Bounds=[-100,100]; f_bias=-300;end
    if f == 15;Bounds=[-5,5];     f_bias= 120;end
    if f == 16;Bounds=[-5,5];     f_bias= 120;end
    if f == 17;Bounds=[-5,5];	  f_bias= 120;end
    if f == 18;Bounds=[-5,5];	  f_bias=10;end
    if f == 19;Bounds=[-5,5];	  f_bias=10;end
    if f == 20;Bounds=[-5,5];     f_bias=10;end
    if f == 21;Bounds=[-5,5];     f_bias=360;end
    if f == 22;Bounds=[-5,5];	  f_bias=360;end
    if f == 23;Bounds=[-5,5]; 	  f_bias=360;end
    if f == 24;Bounds=[-5,5];	  f_bias=260;end
    if f == 25;Bounds=[-2,5]; 	  f_bias=260;end
    w = 0.72;
    Ub = Bounds(2);
    Lb = -Ub;    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity    nop = 10; dim = 9;
    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    x(:,end+1)=0;
    c1 = 1.49445;    c2 = 1.49445;
    x(:,end) = benchmark_func(x(:,1:end-1),f)-f_bias;
    pbest = x; %initialize Best Particle Position
    [mn ind] = min(x(:,end));
    gbest = x(ind,:);
    gbest_hist = [];
    fe = 0;
    for i=1:endgen
        for j=1:nop
            if (x(j,end)<pbest(j,end))
                pbest(j,:) = x(j,:);
            end
            if (pbest(j,end)<gbest(end))
                gbest = pbest(j,:);
            end 
            spd(j,:) = w.* spd(j,:)+c1.*rand(1,dim).*(pbest(j,1:end-1)-x(j,1:end-1))...
                                   +c2.*rand(1,dim).*(gbest(1:end-1) - x(j,1:end-1));
%             ind = find(abs(spd>Vmax));
%             spd(ind) = Vmax.*sign(spd(ind));
            spd=(spd>Vmax).*Vmax+(spd<=Vmax).*spd; 
            spd=(spd<(Vmin)).*(Vmin)+(spd>=(Vmin)).*spd;

            x(j,1:end-1) = x(j,1:end-1)+spd(j,:);
        end
        fe = fe+nop;
%         x(:,end) = abs(benchmark_func(x(:,1:end-1),f)-f_bias);
        x(:,end) = benchmark_func(x(:,1:end-1),f)-f_bias;
%         if mod(i,200) == 0,fprintf('iteration=%u,gbest=%e\n',i,gbest(end));end
        gbest_hist = [gbest_hist gbest(end)];
        fprintf('fun=%u,Gene=%u,Fit_Eval=%u,Gbest=%e\n',f,i,fe,gbest(end))
    end
% end