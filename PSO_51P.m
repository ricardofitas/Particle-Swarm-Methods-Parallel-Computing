% application of Darwinian PSO - Tillett et al. 2005 - Parallel Off.
% Computing
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
D = main_PSOst.D;
main_PSOst.Sw = 100;
main_PSOst.P=1000;
main_PSOst.v = ones(D,main_PSOst.P) + X_max;
main_PSOst.Nkill = [];
Sw = 100;main_PSOst.flag = 0;
main_PSOst.cont_sw = [1:Sw:P, P+1];
main_PSOst.X= round(rand(D,main_PSOst.P).*(X_max-X_min)+X_min);
main_PSOst.xpbest = main_PSOst.X;
main_PSOst.m_max = 400;main_PSOst.m_min = 25;
main_PSOst.perc_x = 0;main_PSOst.percfact=0.6;main_PSOst.pfact=0.2;
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
pbest = main_PSOst.pbest;
xpbest = main_PSOst.xpbest;
gbest = main_PSOst.gbest;
xgbest = main_PSOst.xgbest;
cont_sw = main_PSOst.cont_sw;
%% code
K_K= 2/abs(2-c1-c2-sqrt((c1+c2)^2-4*(c1+c2)));
main_PSOst = Part2(main_PSOst);
g = (abs(pbest)> abs(f1));
pbest(g)=f1(g);
xpbest(:,g) = X(:,g);

for nsw = 1:size(cont_sw,2)-1 % for each swarm in the collection
    if nsw < size(cont_sw,2)
        gbest2 = min(abs(pbest(cont_sw(nsw):(cont_sw(nsw+1)-1))));
        if iter>0
            if gbest(nsw)>gbest2
                gbest(nsw) = gbest2;
                string1 = "Found a new gbest at iteration "+ string(iter)+": "+ string(gbest);
                message2 ="----- Notification from MATLAB ------" + newline + string1 + newline + "Kind regards,"+ newline + "Ricardo Fitas";
                message = "MATLAB: Notification of new gbest";
%                 run('sendemail.m')
                [~,k]=min(pbest(cont_sw(nsw):(cont_sw(nsw+1)-1)));
                xgbest(:,nsw) = X(:,cont_sw(nsw) + k - 1);
                Nkill(nsw) = 0;
                if cont_sw(nsw+1)-cont_sw(nsw) <=m_max
                    % creation of a particle in the current swarm
                    X = [X(:,1:cont_sw(nsw)), round(rand(D,1)*(X_max-X_min)+X_min),...
                        X(:,(cont_sw(nsw)+1):end)];
                    pbest = [pbest(1:cont_sw(nsw)), worst, pbest((cont_sw(nsw)+1):end)];
                    xpbest = [xpbest(:,1:cont_sw(nsw)), X(:,cont_sw(nsw)+1),...
                        xpbest(:,(cont_sw(nsw)+1):end)];
                    v = [v(:,1:cont_sw(nsw)), X_max*ones(D,1), v(:,(cont_sw(nsw)+1):end)];
                    cont_sw((nsw+1):end)=cont_sw((nsw+1):end)+1;
                    curr_swarm = size(X(:,cont_sw(nsw):(cont_sw(nsw+1)-1)),2);                           
                    P = P +1;
                    if size(cont_sw,2) < rand(1,1)/rand(1,1)
                        % creation of a new swarm
                        curr_swarm = size(X(:,cont_sw(nsw):(cont_sw(nsw+1)-1)),2);
                        raios = rand(1,curr_swarm);
                        Xc = mean(X(1,cont_sw(nsw):(cont_sw(nsw+1)-1)));
                        rsort = sort(norm(Xc-X(:,cont_sw(nsw):(cont_sw(nsw+1)-1))));
                        r_max = prctile(rsort,percfact*100)*pfact;
                        xiN_sw = raios.*r_max.*cos(rand(1,curr_swarm)*2*pi)...
                            +X(1,cont_sw(nsw):(cont_sw(nsw+1)-1));
                        yiN_sw = raios.*r_max.*sin(rand(1,curr_swarm)*2*pi)...
                            +X(2,cont_sw(nsw):(cont_sw(nsw+1)-1));
                        X=[X, [xiN_sw;yiN_sw]];
                        v = [v, ones(D,curr_swarm) + X_max];
                        P  = P + curr_swarm;
                        cont_sw = [cont_sw(1:(end-1)), cont_sw((end-1)) + curr_swarm, P+1];
                        pbest = [pbest, zeros(1,curr_swarm) + worst];
                        gbest = [gbest, worst];
                        xgbest = [xgbest, X(:,cont_sw(nsw))];
                        xpbest = [xpbest, [xiN_sw;yiN_sw]];
                    end


                end

            else
                % deleting a particle from the current swarm
                [~,k] = max(pbest(cont_sw(nsw):(cont_sw(nsw+1)-1)));
                X(:,cont_sw(nsw) + k - 1) = [];
                v(:,cont_sw(nsw) + k - 1) = [];
                pbest(cont_sw(nsw) + k - 1)=[];
                xpbest(:,cont_sw(nsw) + k - 1) = [];
                cont_sw((nsw+1):end) = cont_sw((nsw+1):end) -1;
                Nkill(nsw) = Nkill(nsw) +1 ;
                curr_swarm = size(X(:,cont_sw(nsw):(cont_sw(nsw+1)-1)),2);
                P  = P -1;                      
                if cont_sw(nsw+1) -cont_sw(nsw) <m_min
                    % deleting the entire swarm
                    flag = 1;
                    X(:,cont_sw(nsw):(cont_sw(nsw+1)-1)) = [];
                    v(:,cont_sw(nsw):(cont_sw(nsw+1)-1)) = [];
                    xpbest(:,cont_sw(nsw):(cont_sw(nsw+1)-1)) = [];
                    P  = P - curr_swarm;
                    pbest(cont_sw(nsw):(cont_sw(nsw+1)-1))=[];
                    gbest(nsw) = [];
                    xgbest(:,nsw)=[];
                    cont_sw((nsw+1):end) = cont_sw((nsw+1):end) -curr_swarm;
                    cont_sw(nsw) =[];
                    nsw = nsw -1;
                end
            end
        else
            % creation of swarm gbest for iteration 0
            gbest = [gbest, gbest2];
            [gbest3,k]=min(pbest(cont_sw(nsw):(cont_sw(nsw+1)-1)));
            worst = max(pbest(cont_sw(nsw):(cont_sw(nsw+1)-1)));
            xgbest = [xgbest, X(:,cont_sw(nsw) + k - 1)];
            Nkill =[Nkill, 0];
            curr_swarm = size(X(:,cont_sw(nsw):(cont_sw(nsw+1)-1)),2);

        end
        if flag == 0 % if the swarm was not deleted
            r1 = rand(D,curr_swarm);
            r2 = rand(D,curr_swarm);
            v(:,cont_sw(nsw):(cont_sw(nsw+1)-1)) =...
                K*(v(:,cont_sw(nsw):(cont_sw(nsw+1)-1))...
                +c1.*r1.*(xpbest(:,cont_sw(nsw):(cont_sw(nsw+1)-1))...
                -X(:,cont_sw(nsw):(cont_sw(nsw+1)-1))) + c2.*r2.*...
                (xgbest(:,nsw)-X(:,cont_sw(nsw):(cont_sw(nsw+1)-1))));
            X(:,cont_sw(nsw):(cont_sw(nsw+1)-1)) = X(:,cont_sw(nsw):(cont_sw(nsw+1)-1))...
                + v(:,cont_sw(nsw):(cont_sw(nsw+1)-1));
            X = round(min(X_max, max(X_min,X)));
        else
            flag = 0;
        end
    end
end

% update global best (stated as gbest3) between all nbest (vectorized as gbest)
[gbest4,k3] = min(gbest);
if gbest3 > gbest4
    gbest3 = gbest4;
    xgbest3 = xgbest(:,k3);
    Scherwinkel = main_PSOst.Scherwinkelmat(:,:,k3);
    XT = main_PSOst.XTmat(:,:,k3);
    YT = main_PSOst.YTmat(:,:,k3);
    ZT = main_PSOst.ZTmat(:,:,k3);
end
if isempty(xgbest3)
    xgbest = [0;0];
    gbest = 0;
else
    xgbest = xgbest3;
end
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


       