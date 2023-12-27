ParA=[1327/500; 7547; 22350; 1588];
ParE=[1340/500*1.5; 7794/2; 20830];
ParI=[1317/500; 17630/2; 21350];
ParO=[1295/500; 7053; 4459/500; 615; 8388/4; 24700; 98920];

Par=[ParA;ParE;ParI;ParO];
Par=[2.79052842151293;9013.83465468893;22058.6387511091;1605.32744700669;3.38171525558151;4137.70843365706;20170.7866670707;2.70925139764484;9187.69688973580;19520.6245375213;1.26514464608957;1620.27952571195;9.67336292372134;534.779106073429;2226.00287358831;5628.11935821274;79478.2456566545]
Par=[2.72446991166585;8262.43190770059;41063.3934559058;1511.57381315394;3.85384727615157;5776.06245080458;22948.8185329776;3.06091277266163;13425.5358639202;30725.7475085698;0.248637396749928;955.312260007463;9.92059737887939;211.456864654175;2877.04992313086;18076.8323513574;71312.7960658305]
Par1=[2.68726664436153;7800.86393669155;137402.445081227;1448.23578796063;3.82767700554393;5615.90360917634;22695.5387473684;2.78417200742962;10075.6416251211;21650.7779743723;0.980443061841845;2212.34113726344;9.05844733737339;486.025008302109;2389.20922994853;168920.352089300;59704.8588043371]
Par=[2.74927181906818;7975.12482188265;24994.0289640299;1438.45037696741;3.92996859517995;5966.17944288850;23129.9163732702;2.84236238984592;10749.4391326233;23331.8171232680;0.483280585160226;1517.58047814951;9.50552991268020;326.497690517203;2682.60258077200;167240.063195120;64316.9511371378]

Par=[2.73617416822943;7598.11912789107;25277.8384344739;1376.04738153643;3.88000690205551;7583.00127365406;48974.4688480864;2.70482475407282;9193.24733643968;19713.1798812627;2.28282098395244;3523.17614635357;10.0566178364343;2484.69090420343;1507.43652703624;1.70132123400989e+18;15618.2802962940]

old_lst=lst_ODE_model(Par,1)
drawnow

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% new_lst=lst_ODE_model(new_par,1)
% new_lst-old_lst

% clear Lsts Pars
% new_par=Par;
% for n=1:250
% options = optimset('Display','final','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,new_par,options,0);
% Lsts(n)=lst_ODE_model(new_par,1);
% Pars(:,n)=new_par;
% figure(100); plot(Lsts);title('MFA_ODE')
% drawnow
% save('MFA_ODE')
% end


return

function z=lst_ODE_model(Par,g)
ok=0; 
Par=abs(Par);

% load('Met_PathwayS80000noDil')
load('Met_Pathwayv2_S80000_P20000_Q20000_noDil')
temp=output(:,5);
output(:,5)=output(:,6);
output(:,6)=temp;

Tspan=0:100000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% [t,y] = ode15s(@MetPathODE,Tspan,[80000 0 0 0 0],options,Par);
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
% y(end,:)


end

dadosS=interp1(t,y,output(:,1));
lst=(dadosS-output(:,2:end)).^2;
lst(isnan(lst))=0;
lst1=sum(sum(lst));
%%

load('Met_Pathwayv2_S80000_P20000_Q20000_Dil0005In1.mat')
temp=output(:,5);
output(:,5)=output(:,6);
output(:,6)=temp;

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
    
    zz1=sum(y(1:50001,5)*0.0005)+sum(y(1:50001,4)*0.0005);
    zz2=sum(y(1:50001,1)*0.0005);
    zz3=zz2/zz1
    
    za1=sum(output(1:50001,5+1)*0.0005)+sum(output(1:50001,4+1)*0.0005);
    za2=sum(output(1:50001,1+1)*0.0005);
    za=za2/za1
    assignin('base','output',output(1:2000:end,:))
assignin('base','ty',[t(1:500:end) y(1:500:end,:)])
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
function dxdt=MetPathODE2(t,x,p) %This is the ODE function

FA=p(1)*x(1)/p(2)/(1+x(1)/p(2)+x(2)/p(3)+x(4)/p(4));
FE=p(5)*x(2)/p(6)/(1+x(2)/p(6)+x(3)/p(7));
FI=p(8)*x(3)/p(9)/(1+x(3)/p(9)+x(4)/p(10));
FO=(p(11)*x(2)/p(12)+p(13)*(x(4)/p(14))*(x(2)/p(15)))/(1+x(2)/p(12)+x(5)/p(16)+(x(4)/p(14))*(1+x(2)/p(15)+x(5)/p(17)));
in=.09*0; %normpdf(t,.5e5,10000)*500
d=0.0005*0;

dxdt=[-FA+in-d*x(1)        %1 S
       FA-FE-FO-d*x(2)  %2 P
       FE-FI-d*x(3)     %3 Q
       FI-d*x(4)    %4 R
       FO-d*x(5)   %5 T
       0 %6
       0 %7
       ];

end