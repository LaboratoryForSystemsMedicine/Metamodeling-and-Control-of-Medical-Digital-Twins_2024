
Par=[0.00640830930585624;0.170690546476125;0.00278321488398707;0.0482121163635228;1.02999417417886;0.180583726676300;0.181424298656094;-0.269876072057164;-0.233033986220620;0.496038664790223;-0.201049956929962;0.504936421642930;0.111207631243911;-0.172672424662939;0.926126513433303;0.203163808160835;-0.406423140365174;-0.287521472086603;0.167563882016180;-0.582449263855676;-0.192636122480694;-0.0413164153593549;-0.495779960824532;0.652922322995419];

lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)

return

function z=lst_ODE_model(Par,g)
load('MetPw_TrainingDatasets')
ok=0;
Par(1:4)=abs(Par(1:4));

output=[0:1:length(temp1)-1]';
output=[output temp1];

Tspan=100:100000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@odeEvents);
tic
[t,y] = ode15s(@MetPathODE,Tspan,output(Tspan(1)+1,2:6),options,Par);
if or(t(end)<Tspan(end),any(imag(y)>0))
    ok=ok+1e12;
end

if g==1
    figure(10)
    newcolors = [0.83 0.14 0.14
        1.00 0.54 0.00
        0.47 0.25 0.80
        0.25 0.80 0.54
        0.54 0.54 0.54];
    
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(output(:,1),output(:,2:end),'--', 'LineWidth', 1.5)
    hold off
    legend('S','P','Q','R','T')
    title('data mean simu output no dilution')
    
end

dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst1=sum(sum(lst));
%%

output=[0:1:length(temp2)-1]';
output=[output temp2];

Tspan=100:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE_flow,Tspan,output(Tspan(1)+1,2:6),options,Par,[1.0 0.0005]);
if t(end)<Tspan(end)
    ok=ok+1e12;
end

if g==1
    figure(11)
    newcolors = [0.83 0.14 0.14
        1.00 0.54 0.00
        0.47 0.25 0.80
        0.25 0.80 0.54
        0.54 0.54 0.54];
    
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(output(:,1),output(:,2:end),'--', 'LineWidth', 1.5)
    hold off
    title('simu output dilution feedmean=1.0')
    
end

dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst2=sum(sum(lst));
if any(y<0,'All')
    lst2=lst2*10;
    %         'negative'
end

%%
output=[0:1:length(temp3)-1]';
output=[output temp3];

Tspan=100:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE_flow,Tspan,output(Tspan(1)+1,2:6),options,Par,[0.2 0.0005]);
if t(end)<Tspan(end)
    ok=ok+1e12;
end

if g==1
    figure(12)
    newcolors = [0.83 0.14 0.14
        1.00 0.54 0.00
        0.47 0.25 0.80
        0.25 0.80 0.54
        0.54 0.54 0.54];
    
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(output(:,1),output(:,2:end),'--', 'LineWidth', 1.5)
    hold off
    title('simu output dilution feedmean=0.2')
   
end

dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst3=sum(sum(lst));
if any(y<0,'All')
    lst3=lst3*10;
    %         'negative'
end
%%
if g==1
    [lst1 lst2 lst3]
end

z=lst1+lst2+lst3+ok;
end

function dxdt=MetPathODE(t,x,p) 

MS=[-1,0,0,0;1,-1,0,-1;0,1,-1,0;0,0,1,0;0,0,0,1];
F=[p(1)*x(1)^p(5)*x(2)^p(9)*x(3)^p(13)*x(4)^p(17)*x(5)^p(21)  %1
    p(2)*x(1)^p(6)*x(2)^p(10)*x(3)^p(14)*x(4)^p(18)*x(5)^p(22)  %2
    p(3)*x(1)^p(7)*x(2)^p(11)*x(3)^p(15)*x(4)^p(19)*x(5)^p(23)  %3
    p(4)*x(1)^p(8)*x(2)^p(12)*x(3)^p(16)*x(4)^p(20)*x(5)^p(24)];  %4

dxdt=MS*F;

end
function dxdt=MetPathODE_flow(t,x,p,flow) 
in=flow(1);
dil=flow(2);
MS=[-1,0,0,0;1,-1,0,-1;0,1,-1,0;0,0,1,0;0,0,0,1];
F=[p(1)*x(1)^p(5)*x(2)^p(9)*x(3)^p(13)*x(4)^p(17)*x(5)^p(21)  %1
    p(2)*x(1)^p(6)*x(2)^p(10)*x(3)^p(14)*x(4)^p(18)*x(5)^p(22)  %2
    p(3)*x(1)^p(7)*x(2)^p(11)*x(3)^p(15)*x(4)^p(19)*x(5)^p(23)  %3
    p(4)*x(1)^p(8)*x(2)^p(12)*x(3)^p(16)*x(4)^p(20)*x(5)^p(24)];  %4

dxdt=MS*F;

dxdt(1)=dxdt(1)+in;
dxdt=dxdt-x*dil;

end
function [value,isterminal,direction] = odeEvents(t,y,Par,SSa)
x=toc-10;
if x<0
    x=0;
end

value = x;     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;
end
