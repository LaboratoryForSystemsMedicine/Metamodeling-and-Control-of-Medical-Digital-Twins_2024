%S-system
Par=[2069.42893661126;6.76828521866173e-05;1.99765759524512e-07;6.26895926429618e-06;2.73072268196848e-05;0.0175520277401400;-0.0537200338692112;1.04354163421149;0.910081754956121;0.259710658582638;0.129648272522659;0.269749920753349;0.131384586573752;0.899446827103169;1.45427430908421;1.53480975269071;0.716155230053845;-0.260292645375453;-0.0928360230258138;-0.123390351722462;-0.0845052589143311;0.242013107221230;1.01198019045722;1.04445089043740];

lst_ODE_model(Par,1)

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)


return

function z=lst_ODE_model(Par,g)

Par(1:6)=abs(Par(1:6));
ok=0;

load('Dados01')
Tspan=50:1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:) ],options,Par);
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
load('Dados02')
Tspan=1051:2000;

[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:)],options,Par);
if t(end)<Tspan(end)
    ok=ok+1e12;
end

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
     
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,[0.020 0 0]);
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
     
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,[0 0.020 0]);
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
     
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,[0 0 0.015]);
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


if g==1 disp([lst1 lst2  lst4 lst5 lst6]);end

z=lst1+lst2+lst4+lst5+lst6+ok;
end
function dxdt=ODE_SWG(t,x,p) %This is the ODE function

Pa=reshape(p,[6 4]);
Pr=nan(6,4);
Pr(:,1)=Pa(:,1);
Pr(:,2)=x(1).^Pa(:,2);
Pr(:,3)=x(2).^Pa(:,3);
Pr(:,4)=x(3).^Pa(:,4);

Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)
      Pr(3)-Pr(4)
      Pr(5)-Pr(6)];

end
function dxdt=ODE_SWG_wC2(t,x,p,uu) %This is the ODE function

u=[0 0 0];
if t>1000 u=uu; end

Pa=reshape(p,[6 4]);
Pr=nan(6,4);
Pr(:,1)=Pa(:,1);
Pr(:,2)=x(1).^Pa(:,2);
Pr(:,3)=x(2).^Pa(:,3);
Pr(:,4)=x(3).^Pa(:,4);

Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)-u(1)*x(1)
      Pr(3)-Pr(4)-u(2)*x(2)
      Pr(5)-Pr(6)-u(3)*x(3)];

end