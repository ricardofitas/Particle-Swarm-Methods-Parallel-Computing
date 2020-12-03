% application of Hybrid R. Particle Swarm Optimisation + Firefly Algorithm (FAPSO)
% Parallel Off. Computing - parallel swarms for each method
% Ricardo Fitas, 2020

%% initial code
try
    myCluster = parcluster('local');
    parpool(myCluster.NumWorkers)
end
main_PSOst = evalin('base','main_PSOst');
if main_PSOst.selectfunc == 1
    main_PSOst = Part1(main_PSOst);
end
if main_PSOst.selectfunc == 2
    main_PSOst = Part2(main_PSOst);
end
if main_PSOst.selectfunc == 3
    main_PSOst = Part3(main_PSOst);
end
assignin('base','main_PSOst',main_PSOst)

%% functions
function main_PSOst = Part1(main_PSOst)
%% handle of variables
X_min = main_PSOst.X_min;
X_max = main_PSOst.X_max;
P = main_PSOst.P;
D = main_PSOst.D;
main_PSOst.alfa = 0.2; main_PSOst.gam = 1;
main_PSOst.beta0 = 1; main_PSOst.beta_min =0.1;

%% code
main_PSOst.X= round(rand(D,P).*(X_max-X_min)+X_min);
main_PSOst.xpbest = main_PSOst.X;
end

function main_PSOst = Part2(main_PSOst)
%% handle of variables
K = main_PSOst.K;
X = main_PSOst.X;
P = main_PSOst.P;
XG = main_PSOst.XG;
YG = main_PSOst.YG;
ZG = main_PSOst.ZG;
f1 = main_PSOst.f1;
iter = main_PSOst.iter;
detaux = str2func(string(main_PSOst.KDSstr{1}));
user_push1 = evalin('base','user_push1');
%% code
parfor ii = 1:P
    assignin('base','user_push1','user_push1');
    [Scherwinkel,XT,YT,ZT]=detaux(K, X(:,ii), XG, YG, ZG);
    f1(ii) = max(max(abs(Scherwinkel)));
    Scherwinkelmat(:,:,ii) = Scherwinkel;
    XTmat(:,:,ii) =  XT;
    YTmat(:,:,ii) =  YT;
    ZTmat(:,:,ii) =  ZT;
    
end
%% return of handle
if iter == 1
    [main_PSOst.gbest,k] = min(f1);
    main_PSOst.xgbest = X(:,k);
    main_PSOst.XT = XTmat(:,:,k);
    main_PSOst.YT = YTmat(:,:,k);
    main_PSOst.ZT = ZTmat(:,:,k);
    main_PSOst.Scherwinkel = Scherwinkelmat(:,:,k);
end
main_PSOst.f1 = f1;
main_PSOst.XTmat = XTmat;
main_PSOst.YTmat = YTmat;
main_PSOst.ZTmat = ZTmat;
main_PSOst.Scherwinkelmat = Scherwinkelmat;
end

function main_PSOst = Part3(main_PSOst)
%% handle of variables
K = main_PSOst.K;
X = main_PSOst.X;
P = main_PSOst.P;
XG = main_PSOst.XG;
YG = main_PSOst.YG;
ZG = main_PSOst.ZG;
XT = main_PSOst.XT;
YT = main_PSOst.YT;
ZT = main_PSOst.ZT;
Scherwinkel = main_PSOst.Scherwinkel;
f1 = main_PSOst.f1;
X_min = main_PSOst.X_min;
X_max = main_PSOst.X_max;
P = main_PSOst.P;
D = main_PSOst.D;
iter = main_PSOst.iter;
v = main_PSOst.v;
main_PSOst.omega = main_PSOst.o_max - iter*...
    (main_PSOst.o_max-main_PSOst.o_min)/(main_PSOst.maxiter);
c1 = main_PSOst.c1;
c2 = main_PSOst.c2;
alfa = main_PSOst.alfa; gam = main_PSOst.gam;
beta0 = main_PSOst.beta0; beta_min = main_PSOst.beta_min;

%% code
K= 2/abs(2-c1-c2-sqrt((c1+c2)^2-4*(c1+c2)));
part_n = mod(iter,P)+1;
X(:,part_n) = round(rand(D,1)*(X_max-X_min)+X_min);
main_PSOst = Part2(main_PSOst);
f1 = main_PSOst.f1;
gbest = main_PSOst.gbest;
pbest = main_PSOst.pbest;
xpbest = main_PSOst.xpbest;
g = (abs(pbest)> abs(f1));
        pbest(g)=f1(g);
        xpbest(:,g) = X(:,g);
        gbest2 = min(abs(pbest));
        
            if gbest>gbest2
                gbest = gbest2;
                string1 = "Found a new gbest at iteration "+ string(iter)+": "+ string(gbest);
                message2 ="----- Notification from MATLAB ------" + newline + string1 + newline + "Kind regards,"+ newline + "Ricardo Fitas";
                message = "MATLAB: Notification of new gbest";
%                 run('sendemail.m')
                [~,k]=min(pbest);
                xgbest = X(:,k);
                Scherwinkel = main_PSOst.Scherwinkelmat(:,:,k);
                XT = main_PSOst.XTmat(:,:,k);
                YT = main_PSOst.YTmat(:,:,k);
                ZT = main_PSOst.ZTmat(:,:,k);
            end
        
        pbesta = pbest(1:round(size(f1,2)/2));
        xpbesta = xpbest(:,1:round(size(f1,2)/2));
        xpbestb = xpbest(:,(round(size(f1,2)/2)+1):end);
        Xa = X(:,1:round(size(f1,2)/2));
        Xb = X(:,(round(size(f1,2)/2)+1):end);
        vb = v(:,(round(size(f1,2)/2)+1):end);
        Ia = (1+pbesta).^(-1);
        for ik = 1:size(Ia)
            for jk = 1:(ik-1)
                if Ia(jk)> Ia(ik)
                    r = norm(Xa(:,jk)-Xa(:,ik));
                    Xa(:,ik) = Xa(:,ik) + ((beta0 - beta_min)*exp(-gam*r^2)...
                                +beta_min)*(Xa(:,jk)-Xa(:,ik))+...
                                alfa*(rand(1,1)-0.5);
                end
            end
        end
               
        r1 = rand(1,size(Xb,2));
        r2 = rand(1,size(Xb,2));
        vb = K*(vb+c1.*r1.*(xpbestb-Xb) + c2.*r2.*(xgbest-Xb));
        Xb = Xb + vb;
        X =[Xa,Xb];
        X = round(min(X_max, max(X_min,X)));
        rnd1 = rand(1,size(X,2));
        [~,coord] = sort(rnd1);
        X = X(:,coord);
        v = v(:,coord);pbest= pbest(:,coord);xpbest=xpbest(:,coord);
               
%% return of handle
main_PSOst.X = X;
main_PSOst.v = v;
main_PSOst.gbest =gbest;
main_PSOst.pbest =pbest;
main_PSOst.xpbest =xpbest;
main_PSOst.xgbest = xgbest;
main_PSOst.XT = XT;
main_PSOst.YT = YT;
main_PSOst.ZT = ZT;
main_PSOst.Scherwinkel = Scherwinkel;
end
