
Par=[-9.08325271900039e-06;1.89060108857412e-06;1.11281444410269e-08;6.51441077566018e-07;3.41494288640205e-06;-5.52492957539392e-05;-0.000455120317952113;0.000164908382137680;-3.35809852601176e-05;0.000363360499075473;-2.24346928582219e-05;-1.98018619461469e-05;-8.67855816532674e-05;9.03965567725307e-05;5.37407514099764e-05;1.52823939247211e-05;1.53027725528856e-06;3.01593755647163e-06;-2.93008121792157e-06;5.42467581908301e-06;-3.40746754247374e-06;-2.56618061947868e-06;-2.47120370437110e-06;6.00802844781644e-07;-4.99217091464651e-06]
  


lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)


return

function z=lst_ODE_model(Par,g)
ok=0; lst1=0; lsti=0; lstx=0;
ini=0;

load('Met_Pathwayv2_S80k_P20k_Q20k_noDil');

Tspan=ini:100000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@odeEvents);
tic

[t,y] = ode15s(@MetPathODE,Tspan,output(ini+1,2:end),options,Par); 
if or(t(end)<Tspan(end),any(imag(y)>0))
    ok=ok+1e14;
end    
   if ok==0 
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
    if any(y<0,'All')
        lst1=lst1*10;
        'negativos'
    end
   end
   %%

load('Met_Pathwayv2_S80k_P20k_Q20k_Dil0005In1.mat');


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
    assignin('base','output',output(1:2000:end,:))
assignin('base','ty',[t(1:500:end) y(1:500:end,:)])   
end

dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst2=sum(sum(lst));
if any(y<0,'All')
        lst2=lst2*10;
        'negativos'
end

%%

Tspan=0:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE_flow,Tspan,[80000 20000 20000 10 10],options,Par,[.1 0.0005]);
if t(end)<Tspan(end)
    ok=ok+3e12
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
    ylim([-1000 4000])
    
end
if any(y<0,'All')
        lstx=lstx+sum(sum(((y<0).*y).^2));
%         'negativos .1'
end
%%

Tspan=0:50000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,y] = ode15s(@MetPathODE_flow,Tspan,[80000 20000 20000 10 10],options,Par,[.5 0.0005]);
if t(end)<Tspan(end)
    ok=ok+4e12
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
    
end
if any(y<0,'All')
        lstx=lstx+sum(sum(((y<0).*y).^2));
%         'negativos .5'
end
%%

%%


if g==1
    [lst1 lst2 lstx*10]
end

z=lst1+lst2+lstx*10+ok;
end

function dxdt=MetPathODE(t,x,p) %This is the ODE function
ss=[0.0000   0.0000   0.0000    4.3226    7.6794]*1e4;
Pr=reshape(p(1:25),5,5);

X=(x-ss');

dxdt=Pr*X;

end
function dxdt=MetPathODE_flow(t,x,p,flow) %This is the ODE function
in=flow(1);
dil=flow(2);
ss=[0.0000   0.0000   0.0000    4.3226    7.6794]*1e4;
Pr=reshape(p(1:25),5,5);

X=(x-ss');

dxdt=Pr*X;

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
