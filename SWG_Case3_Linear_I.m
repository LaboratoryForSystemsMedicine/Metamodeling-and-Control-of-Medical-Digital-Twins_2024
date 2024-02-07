
Par=[-0.0531680410249922;-0.0464706312973534;-0.00785451223377784;0.000133489352694774;-0.270427611528254;-0.00479917470584026;0.333916097071607;-0.227361444830694;0.0175969929695567]

lst_ODE_model(Par,1)
drawnow; pause(2)

% options = optimset('Display','iter','MaxFunEvals',3000,'MaxIter',3000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)

return

function z=lst_ODE_model(Par,g)

ok=0;lst1=0;lst2=0;lstl1=0;

load('Dados01')
SSa=mean(Ya(900:1000,:));


Tspan=50:1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:) ],options,Par,SSa);
if t(end)<Tspan(end)
    ok=1e12
end
Yss=y(end,:);
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
NP=sum(sum(~isnan(lst)))/3;
lst(isnan(lst))=0;
lst1=sum(sum(lst));
%%
load('Dados02')

SSa=mean(Ya(1900:2000,:));

Tspan=1051:2000;
     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:)],options,Par,SSa);
if t(end)<Tspan(end)
    ok=1e12
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
%     y(end,:)
end
dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;

lst(isnan(lst))=0;
lst2=sum(sum(lst));
%%


if g==1 disp([lst1 lstl1 ok]);end

z=lst1+lstl1+ok;

end

function dxdt=ODE_SWG(t,x,p,ss) %This is the ODE function
K=255*255;

Pr=reshape(p,3,3);

dxdt=Pr*(x-ss');

end
function dxdt=ODE_SWG_wC(t,x,p,ss) %This is the ODE function
K=255*255;

u=0;
if t>1500 u=.01; end

Pr=reshape(p,3,3);

dxdt=Pr*(x-ss');

dxdt(3)=dxdt(3)-u*x(3);

end
