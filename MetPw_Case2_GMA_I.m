
Par=[0.00270595874395440;0.0915455624778722;0.000915238403295641;0.00715902548238057;1.30417107348024;-0.0434933173786568;0.384105758331607;-0.0668968875758788;-0.254866528554665;0.547471900239930;-0.213151604687269;0.435700026116517;-0.0630063459019779;0.0573501028880305;0.826247438230893;0.341883111848523;-0.529463464045609;-0.174500728507968;0.151467465895348;-0.605086547526893;-0.112718203207009;-0.0960564150272877;-0.491284919860753;0.550286514355907]


lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)


return

function z=lst_ODE_model(Par,g)
ok=0; 
Par(1:4)=abs(Par(1:4));

load('Met_Pathwayv2_S80k_P20k_Q20k_noDil')

Tspan=0:100000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@odeEvents);
tic
[t,y] = ode15s(@MetPathODE,Tspan,[80000 20000 20000 10 10],options,Par);
if or(t(end)<Tspan(end),any(imag(y)>0))
    ok=ok+1e12
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

end

dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst1=sum(sum(lst));
%%

load('Met_Pathwayv2_S80k_P20k_Q20k_Dil0005In1.mat')

Tspan=0:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE_flow,Tspan,[80000 20000 20000 10 10],options,Par,[1 0.0005]);
if t(end)<Tspan(end)
    ok=ok+1e12
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
   
end

dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst2=sum(sum(lst));
if any(y<0,'All')
        lst2=lst2*10;
%         'negativos'
    end
%%

%%
if g==1
    [lst1 lst2 ]
end

z=lst1+lst2+ok;
end

function dxdt=MetPathODE(t,x,p) %This is the ODE function

MS=[-1,0,0,0;1,-1,0,-1;0,1,-1,0;0,0,1,0;0,0,0,1];
F=[p(1)*x(1)^p(5)*x(2)^p(9)*x(3)^p(13)*x(4)^p(17)*x(5)^p(21)  %1
   p(2)*x(1)^p(6)*x(2)^p(10)*x(3)^p(14)*x(4)^p(18)*x(5)^p(22)  %2
   p(3)*x(1)^p(7)*x(2)^p(11)*x(3)^p(15)*x(4)^p(19)*x(5)^p(23)  %3
   p(4)*x(1)^p(8)*x(2)^p(12)*x(3)^p(16)*x(4)^p(20)*x(5)^p(24)];  %4

dxdt=MS*F;

end
function dxdt=MetPathODE_flow(t,x,p,flow) %This is the ODE function
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
