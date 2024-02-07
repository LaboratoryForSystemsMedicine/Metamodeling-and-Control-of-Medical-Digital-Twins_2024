
Par=[2.73617416822943;7598.11912789107;25277.8384344739;1376.04738153643;3.88000690205551;7583.00127365406;48974.4688480864;2.70482475407282;9193.24733643968;19713.1798812627;2.28282098395244;3523.17614635357;10.0566178364343;2484.69090420343;1507.43652703624;1.70132123400989e+18;15618.2802962940]

old_lst=lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% new_lst=lst_ODE_model(new_par,1)


return

function z=lst_ODE_model(Par,g)
ok=0; 
Par=abs(Par);

load('Met_Pathwayv2_S80k_P20k_Q20k_noDil')

Tspan=0:100000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE,Tspan,[80000 20000 20000 10 10],options,Par);
if t(end)<Tspan(end)
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
%%

z=lst1+lst2+ok;
end

function dxdt=MetPathODE(t,x,p) %This is the ODE function

FA=p(1)*x(1)/p(2)/(1+x(1)/p(2)+x(2)/p(3)+x(4)/p(4));
FE=p(5)*x(2)/p(6)/(1+x(2)/p(6)+x(3)/p(7));
FI=p(8)*x(3)/p(9)/(1+x(3)/p(9)+x(4)/p(10));
FO=(p(11)*x(2)/p(12)+p(13)*(x(4)/p(14))*(x(2)/p(15)))/(1+x(2)/p(12)+x(5)/p(16)+(x(4)/p(14))*(1+x(2)/p(15)+x(5)/p(17)));

dxdt=[-FA        %1 S
       FA-FE-FO  %2 P
       FE-FI     %3 Q
       FI    %4 T
       FO    %5 R
       ];

end
function dxdt=MetPathODE_flow(t,x,p,flow) %This is the ODE function
in=flow(1);
dil=flow(2);

FA=p(1)*x(1)/p(2)/(1+x(1)/p(2)+x(2)/p(3)+x(4)/p(4));
FE=p(5)*x(2)/p(6)/(1+x(2)/p(6)+x(3)/p(7));
FI=p(8)*x(3)/p(9)/(1+x(3)/p(9)+x(4)/p(10));
FO=(p(11)*x(2)/p(12)+p(13)*(x(4)/p(14))*(x(2)/p(15)))/(1+x(2)/p(12)+x(5)/p(16)+(x(4)/p(14))*(1+x(2)/p(15)+x(5)/p(17)));

dxdt=[in-FA-x(1)*dil        %1 S
       FA-FE-FO-x(2)*dil  %2 P
       FE-FI-x(3)*dil     %3 Q
       FI-x(4)*dil    %4 R
       FO-x(5)*dil    %5 T
       ];

end
