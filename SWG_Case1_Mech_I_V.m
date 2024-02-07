
Par=[0.209782211319144;4.50613392290133e-06;2.38108020743778e-05;1.27393787122488e-06;0.00854743569403738;1.37641724324446e-05;1.50182377649022e-05;0.0633574675792873]

lst_ODE_model(Par,1)

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)

return

function z=lst_ODE_model(Par,g)

Par(:)=abs(Par(:));

ok=0;


load('Dados01') % Read data as Ta(time) and Ya(Grass, sheep, wolves)

Tspan=50:1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode45(@ODE_SWG_wC2,Tspan,[Ya(Tspan(1),:) ],options,Par,[0 0 0], +inf);
if t(end)<Tspan(end)
    ok=ok+1e12;
end

if g==1
    newcolors = [0.25 0.80 0.54;0.7 0.7 0.7;0.83 0.14 0.14];
    
    
    figure(1)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    
    hold on
    plot(Ta,Ya,'--', 'LineWidth', 1.5)
    hold off
    legend('grass','sheep','wolves')
    
end
dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst1=sum(sum(lst));
%%
load('Dados02') % Read data as Ta(time) and Ya(Grass, sheep, wolves)
Tspan=1051:2000; % Daraset 'Dados02' was created by shifting the parameters 
                 % from 0-999, letting the system reach a steady-state, and
                 % shifting the parameters back to the original values at
                 % 1000 and let the system reach its normal steady state

[t,y] = ode45(@ODE_SWG_wC2,Tspan,[Ya(Tspan(1),:)],options,Par,[0 0 0],+inf);
if t(end)<Tspan(end)
    ok=ok+1e12;
end
Yss=y(end,:); %Steady-state of the system

if g==1
    figure(2)
    colororder(newcolors)
    plot(Ta(1000:end),Ya(1000:end,:),'--', 'LineWidth', 1.5)
    hold on
    plot(t,y, 'LineWidth', 1.5)
    hold off
    legend('grass','sheep','wolves')
    xlim([1000 2000])
        
end

dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst2=sum(sum(lst));

%%
load('DadosConGrass2')
load('DadosConGrass2s')
Tspan=500:2500;

[t,y] = ode45(@ODE_SWG_wC2,Tspan,[Yss],options,Par,[0.020 0 0],1000);
if t(end)<Tspan(end)
    ok=ok+1e12;
end

if g==1
    figure(4)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:2500),Ya(500:2500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:2500),Ya(500:2500,:)+3*Yas(500:2500,:),'--', 'LineWidth', 1)
    plot(Ta(500:2500),Ya(500:2500,:)-3*Yas(500:2500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 2500])
    
end

dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst4=sum(sum(lst));

%%
load('DadosConSheep2')
load('DadosConSheep2s')
Tspan=500:2500;

[t,y] = ode45(@ODE_SWG_wC2,Tspan,[Yss],options,Par,[0 0.020 0],1000);
if t(end)<Tspan(end)
    ok=ok+1e12;
end

if g==1
    figure(5)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:2500),Ya(500:2500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:2500),Ya(500:2500,:)+3*Yas(500:2500,:),'--', 'LineWidth', 1)
    plot(Ta(500:2500),Ya(500:2500,:)-3*Yas(500:2500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 2500])
    
end

dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst5=sum(sum(lst));

%%
load('DadosConWolves1.5.mat')
load('DadosConWolves1.5s.mat')
Tspan=500:2500;

[t,y] = ode45(@ODE_SWG_wC2,Tspan,[Yss],options,Par,[0 0 0.015],1000);
if t(end)<Tspan(end)
    ok=ok+1e12;
end

if g==1
    figure(6)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:2500),Ya(500:2500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:2500),Ya(500:2500,:)+3*Yas(500:2500,:),'--', 'LineWidth', 1)
    plot(Ta(500:2500),Ya(500:2500,:)-3*Yas(500:2500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 2500])
    
end

dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst6=sum(sum(lst));
%%

if g==1 
    disp([lst1 lst2 lst4 lst5 lst6 ])
end

z=lst1+lst2+lst4+lst5+lst6+ok;
end
function dxdt=ODE_SWG_wC2(t,x,p,uu,tu) %This is the ODE function

u=[0 0 0];
if t>tu u=uu; end

dxdt=[p(1)*x(1)-p(2)*x(1)^2-p(3)*x(1)*x(2)
    p(4)*x(1)*x(2)-p(5)*x(2)-p(6)*x(2)*x(3)
    p(7)*x(2)*x(3)-p(8)*x(3)];

dxdt=dxdt+[-u(1)*x(1)
           -u(2)*x(2)
           -u(3)*x(3)];
end
