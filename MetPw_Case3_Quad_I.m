warning off

Par=[1.08788006796825e-24;1.00016654989399e-24;3.47231170418797e-06;1.29868502414507e-26;4.88624972143440e-28;-4.17483717633520e-20;-0.000791646968960071;3.78428102939598e-26;-9.00341789880448e-21;0.000920850502743092;-3.10276750798453e-33;-2.43216934626296e-22;-0.000100995430466242;8.10605794849896e-05;3.29753226981123e-05;3.32055826683392e-05;-3.74369122962240e-27;-1.57088562256977e-22;-1.23929560944211e-06;1.54854867232725e-31;8.49951215943475e-07;8.24044658911237e-35;2.76438158077051e-28;-1.87921664916293e-17;-1.01351888882252e-05;5.18957853424483e-11;-2.33030566640172e-11;-2.83415912729087e-11;8.08965469008923e-29;-7.33513250917091e-29;7.06171037116604e-30;-3.10867604098722e-09;-4.44682257982137e-26;-1.17434518496350e-09;-9.79590804670437e-19;4.18128223724650e-32;4.03106253764490e-29;-9.38190298181348e-11;4.32700957890714e-18;7.96798618660199e-22;1.04014010916323e-10;-1.45095713284577e-23;-2.62783932452203e-18;1.66888032263158e-25;-6.36187656343929e-11;1.80243047777822e-10;-1.00770051853952e-10;4.90287105811815e-15;-3.36621344192835e-11;-7.92695221695255e-21;2.99372305324953e-09;4.08650203058687e-08;8.97304813167621e-09;3.79900167608768e-09;-4.60274717316329e-08;-4.61161502479196e-09;2.04003762155830e-21;3.67005001655732e-20;9.22647271710370e-27;5.41245271614902e-09;6.62740438420717e-20;-8.24986015792764e-29;6.31945657537015e-30;7.70002682150748e-25;-2.27745804188908e-18;-1.53811983418764e-33;3.53880919495865e-10;2.76374769461554e-27;6.30743783184525e-20;-2.12467596175532e-21;-1.39816386181444e-22;8.64307813027444e-23;1.71164211748089e-09;5.75155505896683e-16;-8.56387303568856e-22;2.50515900937401e-10;-1.10728769602253e-23;-4.09998163241263e-31;-1.23377781252218e-25;2.60938923421032e-10;-7.59838727451868e-11;-8.27412782471277e-31;6.42159530356465e-11;-1.31531900719341e-10;6.03414197708929e-32;1.43022984343540e-25;7.66150647365070e-26;1.06531502407512e-10;-1.00816907770403e-21;-1.25650926429669e-31;3.33895167429608e-10;-2.57285662745157e-11;4.69493444858728e-31;-2.06524814577257e-28;-1.66213638727884e-22;-9.85277392531334e-24;3.83664011175963e-11;-1.74334240517491e-11;6.86752395286793e-19;-1.17739832588341e-10]
  



lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)
% save('new_par','new_par')



return

function z=lst_ODE_model(Par,g)
ok=0; lst1=0; lsti=0;lstx=0;


ini=0;

load('Met_Pathwayv2_S80k_P20k_Q20k_noDil')

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
%         'negativos'
    end
   end
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

if g==1
    [lst1 lst2 lstx*10]
end

z=lst1+lst2+lstx*10+ok;
end

function dxdt=MetPathODE(t,x,p) %This is the ODE function
ss=[0.0000   0.0000   0.0000    4.3226    7.6794]*1e4;
Pr=reshape(p(1:25),5,5);
H=reshape(p(26:100),5,15);

X=(x-ss');
X2=[X(1)*X(1)
    X(1)*X(2)
    X(1)*X(3)
    X(1)*X(4)
    X(1)*X(5)
    X(2)*X(2)
    X(2)*X(3)
    X(2)*X(4)
    X(2)*X(5)
    X(3)*X(3)
    X(3)*X(4)
    X(3)*X(5)
    X(4)*X(4)
    X(4)*X(5)
    X(5)*X(5)];

dxdt=Pr*X+H*X2;

end
function dxdt=MetPathODE_flow(t,x,p,flow) %This is the ODE function
in=flow(1);
dil=flow(2);
ss=[0.0000   0.0000   0.0000    4.3226    7.6794]*1e4;
Pr=reshape(p(1:25),5,5);
H=reshape(p(26:100),5,15);

X=(x-ss');
X2=[X(1)*X(1)
    X(1)*X(2)
    X(1)*X(3)
    X(1)*X(4)
    X(1)*X(5)
    X(2)*X(2)
    X(2)*X(3)
    X(2)*X(4)
    X(2)*X(5)
    X(3)*X(3)
    X(3)*X(4)
    X(3)*X(5)
    X(4)*X(4)
    X(4)*X(5)
    X(5)*X(5)];

dxdt=Pr*X+H*X2;

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
