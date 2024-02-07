
Par=[1101.56066143358;2.95182065883021e-05;9.72853126647962e-07;1.91385037901842;7.84362600961635e-07;1.26570409541109e-05;0.0613654242587333;-0.0611477933738020;1.08667353832568;0.776810343311747;-0.486323213123282;0.0130619420299988;0.0629454415177019;0.0262145077893222;0.150276670828555;0.891745745356395;1.36118282775825;0.830813549752694;1.66332965702456;0.920665169816947;-0.0699308359026207;0.000705888820201751;-0.0418167654121998;0.0299791926048226;0.276134902900581;0.628362758945057;0.972602707658072;1.02455180137201]


lst_ODE_model(Par,1)
drawnow; pause(.1)

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)

return

function z=lst_ODE_model(Par,g)
Par(1:7)=abs(Par(1:7));
Par(8:end)=Par(8:end).*(abs(Par(8:end))>0.03);
ok=0;lst1=0;lst2=0;

load('Dados01')
Tspan=50:1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:) ],options,Par);
if t(end)<Tspan(end)
    ok=ok+1e12
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
    ok=ok+1e12
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
    ok=ok+1e12
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
    ok=ok+1e12
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
    ok=ok+1e12
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

lstl1=(reshape(Par,[7 4]));
lstl1=((lstl1(:,2:4))*1000).^2; 
lstl1(isnan(lstl1))=0;
lstl1=sum(sum(lstl1))*10;

if g==1 disp([lst1 lst2 lst4 lst5 lst6 lstl1 ok]);end

z=lst1+lst2+lstl1+lst4+lst5+lst6+ok;
end
function dxdt=ODE_SWG(t,x,p) %This is the ODE function

K=255*255;


Pn=reshape(1:28,[7 4]);
Pr(:,1)=p(Pn(:,1));
Pr(:,2:4)=repmat([x(1) x(2) x(3)],7,1).^p(Pn(:,2:4));
Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)
      Pr(3)-Pr(4)-Pr(5)
      Pr(6)-Pr(7)];

end
function dxdt=ODE_SWG_wC2(t,x,p,uu) %This is the ODE function

K=255*255;

u=[0 0 0];
if t>1000 u=uu; end

Pn=reshape(1:28,[7 4]);
Pr(:,1)=p(Pn(:,1));
Pr(:,2:4)=repmat([x(1) x(2) x(3)],7,1).^p(Pn(:,2:4));
Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)-u(1)*x(1)
      Pr(3)-Pr(4)-Pr(5)-u(2)*x(2)
      Pr(6)-Pr(7)-u(3)*x(3)];

end