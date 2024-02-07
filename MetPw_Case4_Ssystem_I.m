

Par=[2.01000030504152e-07;6.17840285506539e-06;1.10639139995044e-06;0.000110159137086306;0.113546414312275;0.00453182677331608;0.00787946513524819;4184.60935960298;9.72354096343329e-06;7.45235374915438e-07;0.00201466862736885;1.35806851827259;1.47125608715638;-0.0413332815214503;0.125323151373965;-0.187986414879318;-0.133320369190795;-1.14194906409362e-05;0.0898636908657002;0.00309915005096666;0.00156599346524686;0.000988228141502511;-0.365654395951979;0.742112769922002;1.06197306096321;-0.434502556076263;-0.122077198374000;-0.360712897676352;0.781793251576705;-0.00229359872198341;-2.75543128480960e-06;-0.0728988179221552;0.307587158960334;0.279130397969031;-1.02311605949548;1.24529297006206;0.926938704087483;-1.31534844351156;0.0947342714712098;0.000277660614865001;-0.0157963545891748;-0.464600830019000;-0.546651910705414;0.110322407091817;0.249606169665868;-0.176750216887226;-0.160541997233235;0.114997474041835;0.315865039695119;-0.000108308488959499;0.0149524576062521;0.211847079455893;0.124749438313031;0.103162100471898;-0.169979149549653;0.000994970637311251;0.00130141249825413;-0.275524539554565;0.195649325612846;0.00364796436815905]



lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)



return

function z=lst_ODE_model(Par,g)
ok=0; 
Par(1:10)=abs(Par(1:10));

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
    ok=ok+2e12
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

Tspan=0:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE_flow,Tspan,[80000 20000 20000 10 10],options,Par,[.1 0.0005]);
if t(end)<Tspan(end)
    ok=ok+3e12;
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
    ylim([0 3000])
    if any(imag(y),'All')
        ok=ok+1e13;
    end
     if any(std(y(20000:end,:))>2,'All')
         ok=ok+1e11;
         %[12 std(y(20000:end,:))]
     end
     
    
end
%%

Tspan=0:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE_flow,Tspan,[80000 20000 20000 10 10],options,Par,[.5 0.0005]);
if t(end)<Tspan(end)
    ok=ok+4e12;
end
if g==1
    figure(13)
    newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54
             0.54 0.54 0.54];
         
         colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    ylim([0 3000])
   if any(imag(y),'All')
        ok=ok+1e13;
    end
    if any(std(y(20000:end,:))>2,'All')
         ok=ok+1e11;
         %[13 std(y(20000:end,:))]
     end
end


%%
pp=reshape(Par,10,6);
lstp=abs(pp(:,2:6));


lstp=sum(sum(lstp));
lstp=lstp/10^(floor(log10(lstp)))*10^(floor(log10(lst1))-1);

if g==1
    [lst1 lst2 lstp ok]
    
end

z=lst1+lst2+lstp+ok;
end

function dxdt=MetPathODE(t,x,p) %This is the ODE function

pp=reshape(p,10,6);
pp(abs(pp)<1e-3)=0;
pp(:,1)=p(1:10);

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
pp(abs(pp)<1e-3)=0;
pp(:,1)=p(1:10);

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