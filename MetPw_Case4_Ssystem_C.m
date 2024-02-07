
Par=[0;0.00601188207042535;0.00867385572855826;0.189775327323142;0.135058475430618;0.00583493958898133;0.00441214421772229;0;0.0172260510283459;0;0;1.09825793360593;1.12900331250909;0.0700490576000509;0.0590603780112077;0.0267535325540021;0.0334965297753868;0;-0.0197529233980140;0;0;-0.254893939865513;-0.245865630139353;0.221529998084442;0.624336634260002;-0.157391933118318;-0.160087074419801;70.8418459970096;0.467669514426608;0;0;0.118556873096018;0.234615056294142;0.328176022831244;-0.252379519272501;1.01173845445508;0.994057067284675;0;0.0807042789827792;71.3844796244336;0;-0.634473321392564;-0.647582803803906;0.124236195363149;0.0939037702307317;0.0705188219153652;0.0964523243509500;0;-0.327067115730190;0;0;-0.0479493057401740;-0.0845031885001927;-0.310713490882139;-0.233169377753989;-0.417945346572671;-0.404674308687365;0;0.390083412264538;0];

lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)


return

function z=lst_ODE_model(Par,g)
load('MetPw_TrainingDatasets')

ok=0; lst1=0; lst2=0; lst3=0;
Par(1:10)=abs(Par(1:10));

output=[0:1:length(temp1)-1]';
output=[output temp1];

Tspan=100:100000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@odeEvents);
tic
[t,y] = ode15s(@MetPathODE,Tspan,output(Tspan(1)+1,2:6),options,Par);
if or(t(end)<Tspan(end),any(imag(y)>0))
    ok=ok+1e16
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

if t(end)>=Tspan(end)
dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst1=sum(sum(lst));
end

%%

output=[0:1:length(temp2)-1]';
output=[output temp2];

Tspan=100:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@odeEvents);
[t,y] = ode15s(@MetPathODE_flow,Tspan,output(Tspan(1)+1,2:6),options,Par,[1.0 0.0005]);
if or(t(end)<Tspan(end),any(imag(y)>0))
    ok=ok+1e16
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
 if t(end)>=Tspan(end)
dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst2=sum(sum(lst));
 end
 
%%

output=[0:1:length(temp3)-1]';
output=[output temp3];

Tspan=100:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@odeEvents);
[t,y] = ode15s(@MetPathODE_flow,Tspan,output(Tspan(1)+1,2:6),options,Par,[.2 0.0005]);
if or(t(end)<Tspan(end),any(imag(y)>0))
    ok=ok+1e16
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
    title('simu_output_dilution_feed_mean=0.2_nohdr.dat')
         
end
 if t(end)>=Tspan(end)
dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst3=sum(sum(lst));
 end

%%

if g==1
    [lst1 lst2 lst3 ok]
end

z=lst1+lst2+lst3+ok;
end

function dxdt=MetPathODE(t,x,p) %This is the ODE function

pp=p';

PL=[pp(1:10)' x(1).^pp(11:20)' x(2).^pp(21:30)' x(3).^pp(31:40)' x(4).^pp(41:50)' x(5).^pp(51:60)'];
F=prod(PL,2);  

dxdt=[F(1)-F(2)
      F(3)-F(4)
      F(5)-F(6)
      F(7)-F(8)
      F(9)-F(10)];

end
function dxdt=MetPathODE_flow(t,x,p,flow) %This is the ODE function
in=flow(1);
dil=flow(2);

pp=reshape(p,10,6);


PL=[pp(1:10)' x(1).^pp(11:20)' x(2).^pp(21:30)' x(3).^pp(31:40)' x(4).^pp(41:50)' x(5).^pp(51:60)'];
F=prod(PL,2);  

dxdt=[in+F(1)-F(2)-x(1)*dil
         F(3)-F(4)-x(2)*dil
         F(5)-F(6)-x(3)*dil
         F(7)-F(8)-x(4)*dil
         F(9)-F(10)-x(5)*dil];

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