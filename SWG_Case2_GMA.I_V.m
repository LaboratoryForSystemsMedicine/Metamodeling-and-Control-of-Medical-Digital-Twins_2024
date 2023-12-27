%original from v1
Par=[707.898797713976;2.74290709712722e-05;1.18449032179785e-06;7.63614372391576;2.47724294281499e-06;1.08954552559593e-05;0.129450032776724;-0.131322619148285;1.01002513278854;0.937497531445507;-0.940942480785047;-0.00313911396964099;0.0379856779311007;0.00216129790049945;0.0652059903326178;0.915235668006856;0.988680344035039;0.680196173029810;1.43529013066231;0.998693159011258;0.00973677896764006;0.138430231017208;-0.0624222334509402;0.0585716046376420;0.702043076249022;0.660600296309528;1.03520222002640;0.942070121122811]
%v2
% min kinetic orders
% Par=[897.448326832814;2.71567651594747e-05;1.09050868575415e-06;0.240327944250222;5.00456563021423e-06;9.79248912277158e-06;0.106665013844497;-0.0877196524423893;1.08748471266136;0.717370689177534;-0.293501309629295;-0.00275894701275943;0.0476499522536238;0.00397173286260035;0.227870019122627;0.916520994853719;1.38351082937825;-0.323252799556142;1.56688571521740;0.960343426046831;-0.0410880329175791;0.00888482598546845;-0.0359850222828553;0.00941487806988185;0.598081414083318;0.554784169095061;1.01837638168668;0.953293302987333]
% Par=[931.288439304380;2.70126470469888e-05;1.09741812893163e-06;0.0232677196945657;5.03674289783284e-06;1.09508802334730e-05;0.124098364477390;-0.0838120517317377;1.09041722943914;0.708896454125246;-0.0299108791296472;-0.00126320925924433;0.0470916115634871;-0.00512966549651448;0.230063778938654;0.932982842402985;1.39466882561427;-0.763124556698743;1.59805605851984;0.956485231496215;-0.0608941971558001;0.0264343610275821;-0.0440009924105612;0.0214035200419317;1.39191491856994;0.519821154457042;1.01257179133641;0.958018208669666]
% Par=[855.677254593824;2.86609565035937e-05;1.26699569065173e-06;0.0649440447259673;5.79031265233158e-06;1.06559927123140e-05;0.112193151915659;-0.0861690926380724;1.08132163100991;0.713206193418926;0.0299650141587846;-0.00122423672633126;0.0477321806283407;-0.0185806941293281;0.231549207331267;0.924510378580812;1.37788488747339;-0.385239786990464;1.59965852889060;0.950647473342846;-0.0663299192042424;-0.0299032434285010;-0.0429982451906487;0.0203036887441242;0.984511073666365;0.503250681033503;1.01755715848537;0.973074169376385]

%v2.5 with 2% control Grass and sheep and 1.5% control wolves
% Par=[1016.38741682611;2.89167391782216e-05;7.01025625592541e-07;0.427794193591778;9.85877053188545e-07;1.14824410397981e-05;0.0583830159230668;-0.0573859616447017;1.07252981513768;0.769316788541306;-0.273611272905820;-0.00358720407640660;0.0514770656446575;0.00634045787264099;0.154063176972102;0.915172740032156;1.37265278887422;-0.301604210204683;1.67166755333122;0.945867497668825;-0.0668164888948919;0.00836281968516820;-0.0481336301274764;0.0134460564477915;0.305393966556063;0.657587084900424;0.974219652485646;1.02923546201775]
% Par=[1026.31596480568;2.95723820711493e-05;7.85595256131513e-07;1.78703404207264;1.21311666299074e-06;1.18068610754779e-05;0.0582460465652736;-0.0537600491422890;1.07191139640985;0.777967321456126;-0.462056862223238;0.00550295414095926;0.0490471988648991;0.0298785826964045;0.154775050324967;0.915961227727342;1.38082637565161;0.803478734015288;1.64874430704353;0.938778628275464;-0.0702274833055161;-0.00170026589261001;-0.0431139402584057;0.0249928412239857;0.233292781031578;0.622585382228043;0.978676018157081;1.02928004264924]
Par=[1101.56066143358;2.95182065883021e-05;9.72853126647962e-07;1.91385037901842;7.84362600961635e-07;1.26570409541109e-05;0.0613654242587333;-0.0611477933738020;1.08667353832568;0.776810343311747;-0.486323213123282;0.0130619420299988;0.0629454415177019;0.0262145077893222;0.150276670828555;0.891745745356395;1.36118282775825;0.830813549752694;1.66332965702456;0.920665169816947;-0.0699308359026207;0.000705888820201751;-0.0418167654121998;0.0299791926048226;0.276134902900581;0.628362758945057;0.972602707658072;1.02455180137201]




% close all
lst_ODE_model(Par,1)
drawnow; pause(.1)

% options = optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000);
% new_par=fminsearch(@lst_ODE_model,Par,options,0)
% lst_ODE_model(new_par,1)
% z=reshape(new_par,[7 4]);
% % zp=reshape(Par,[7 4]);


% zp=reshape(Par,[7 4]);
% double(int16(zp(:,2:4)*100))/100

return

function z=lst_ODE_model(Par,g)
Par(1:7)=abs(Par(1:7));
% nzeros=sum(abs(Par(8:end))<=0.03);
Par(8:end)=Par(8:end).*(abs(Par(8:end))>0.03);
ok=0;lst1=0;lst2=0;
% Par([11 12 ])=0;

load('Dados01')
Tspan=50:1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:) ],options,Par);
if t(end)<Tspan(end)
    ok=ok+1e12
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
%     y(end,:)
% y
   
end
dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst1=sum(sum(lst));

load('Dados02')
Tspan=1051:2000;
     
[t,y] = ode45(@ODE_SWG,Tspan,[Ya(Tspan(1),:)],options,Par);
if t(end)<Tspan(end)
    ok=ok+1e12
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
% y

end
dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst2=sum(sum(lst));

if g==1
load('Dados03')
load('Dados03s')
Tspan=1000:2000;
     
[t,y] = ode45(@ODE_SWG_wC,Tspan,[ 2.4236    0.4001    0.1898]*1e4,options,Par);
if t(end)<Tspan(end)
    ok=ok+1e12
end


    figure(3)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(1000:2000),Ya(1000:2000,:),'--', 'LineWidth', 1.5)
    plot(Ta(1000:2000),Ya(1000:2000,:)+3*Yas(1000:2000,:),'--', 'LineWidth', 1)
    plot(Ta(1000:2000),Ya(1000:2000,:)-3*Yas(1000:2000,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([1000 2000])
%     y(end,:)
% y


dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst3=sum(sum(lst));
end

%%
load('DadosConGrass2')
load('DadosConGrass2s')
Tspan=500:2500;
     
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,[0.020 0 0]);
if t(end)<Tspan(end)
    ok=ok+1e12
end

if g==1
    figure(4)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:2500),Ya(500:2500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:2500),Ya(500:2500,:)+3*Yas(500:2500,:),'--', 'LineWidth', 1)
    plot(Ta(500:2500),Ya(500:2500,:)-3*Yas(500:2500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 2500])
    
end

dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst4=sum(sum(lst));

%%
load('DadosConSheep2')
load('DadosConSheep2s')
Tspan=500:2500;
     
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,[0 0.020 0]);
if t(end)<Tspan(end)
    ok=ok+1e12
end

if g==1
    figure(5)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:2500),Ya(500:2500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:2500),Ya(500:2500,:)+3*Yas(500:2500,:),'--', 'LineWidth', 1)
    plot(Ta(500:2500),Ya(500:2500,:)-3*Yas(500:2500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 2500])

end

dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst5=sum(sum(lst));

%%
load('DadosConWolves1.5.mat')
load('DadosConWolves1.5s.mat')
Tspan=500:2500;
     
[t,y] = ode45(@ODE_SWG_wC2,Tspan,Ya(Tspan(1),:),options,Par,[0 0 0.015]);
if t(end)<Tspan(end)
    ok=ok+1e12
end

if g==1
    figure(6)
    colororder(newcolors)
    plot(t,y, 'LineWidth', 1.5)
    hold on
    plot(Ta(500:2500),Ya(500:2500,:),'--', 'LineWidth', 1.5)
    plot(Ta(500:2500),Ya(500:2500,:)+3*Yas(500:2500,:),'--', 'LineWidth', 1)
    plot(Ta(500:2500),Ya(500:2500,:)-3*Yas(500:2500,:),'--', 'LineWidth', 1)
    
    hold off
    legend('grass','sheep','wolves')
    xlim([500 2500])
%         clc
% [t,y]
% %     [Ta,Ya]
%     ihihiiuhi
end

dadosS=interp1(t,y,Ta);
lst=(dadosS-Ya).^2;
lst(isnan(lst))=0;
lst6=sum(sum(lst));
%%



lstl1=(reshape(Par,[7 4]));
lstl1=((lstl1(:,2:4))*1000).^2; %-[NaN,0,0;1,1,0;0,1,0;-1,1,0;0,1,1;0,0,1;0,-1,NaN]
lstl1(isnan(lstl1))=0;
lstl1=sum(sum(lstl1))*10;

if g==1 disp([lst1 lst2 lstl1*0 lst4 lst5 lst6 lstl1]);end

z=lst1+lst2+lstl1+lst4+lst5+lst6+ok;
end
function dxdt=ODE_SWG(t,x,p) %This is the ODE function

K=255*255;


Pn=reshape(1:28,[7 4]);
Pr(:,1)=p(Pn(:,1));
Pr(:,2:4)=repmat([x(1) x(2) x(3)],7,1).^p(Pn(:,2:4));
Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)
      Pr(3)-Pr(4)-Pr(5)
      Pr(6)-Pr(7)];

end
function dxdt=ODE_SWG_wC(t,x,p) %This is the ODE function

K=255*255;

u=0;
if t>1500 u=.01; end

Pn=reshape(1:28,[7 4]);
Pr(:,1)=p(Pn(:,1));
Pr(:,2:4)=repmat([x(1) x(2) x(3)],7,1).^p(Pn(:,2:4));
Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)
      Pr(3)-Pr(4)-Pr(5)
      Pr(6)-Pr(7)-u*x(3)];


end
function dxdt=ODE_SWG_wC2(t,x,p,uu) %This is the ODE function

K=255*255;

u=[0 0 0];
if t>1000 u=uu; end

Pn=reshape(1:28,[7 4]);
Pr(:,1)=p(Pn(:,1));
Pr(:,2:4)=repmat([x(1) x(2) x(3)],7,1).^p(Pn(:,2:4));
Pr=prod(Pr,2);

dxdt=[Pr(1)-Pr(2)-u(1)*x(1)
      Pr(3)-Pr(4)-Pr(5)-u(2)*x(2)
      Pr(6)-Pr(7)-u(3)*x(3)];


end