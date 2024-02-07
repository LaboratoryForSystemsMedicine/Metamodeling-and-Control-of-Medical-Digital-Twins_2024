
Par=[2.64446640558550;6914.65232751257;24119.1733356704;1291.78780305338;3.34356194466207;3594.79818496842;17601.6585503908;2.57390072789384;8489.16090388438;21221.0908879279;0.302156187324374;1568.14536443261;8.97847898824731;351.513289087267;2323.27749676967;2.61913653349181e+26;123792.255792294];

lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% new_lst=lst_ODE_model(new_par,1)


return

function z=lst_ODE_model(Par,g)
load('MetPw_TrainingDatasets')
ok=0; 
Par=abs(Par);

output=[0:1:length(temp1)-1]';
output=[output temp1];

Tspan=100:100000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE,Tspan,output(Tspan(1)+1,2:6),options,Par);
if t(end)<Tspan(end)
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
[t,y] = ode15s(@MetPathODE_flow,Tspan,output(Tspan(1)+1,2:6),options,Par,[1 0.0005]);
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
%%

if g==1
    [lst1 lst2 lst3 ok]
end

z=lst1+lst2+lst3+ok;
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

dxdt=[in-FA-x(1)*dil    %1 S
         FA-FE-FO-x(2)*dil    %2 P
         FE-FI-x(3)*dil    %3 Q
         FI-x(4)*dil    %4 R
         FO-x(5)*dil    %5 T
         ];

end
