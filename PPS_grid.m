%%初期化処理
% f_init = input('initialize?[y,n]','s');  %'y'の時,初期化処理実行
% if isempty(f_init)
%     f_init = 'n' ;
% end
% 
% 
% if f_init == 'y'
    clear all;
%     disp('初期化中...');
    %% 
    
%%グラフの定義###

% エッジ集合にて定義

% 最小グラフ
E = [ 1,2; 2,1];

num_x = max( max(E) );


% 隣接行列 N
    N = zeros([num_x,num_x]);
    for i=1:length(E);
        N( E(i,1) , E(i,2) ) = 1;
        N( E(i,2) , E(i,1) ) = 1;  %%無向グラフの場合
    end
    
    %グラフラプラシアン　L
    L_diag=zeros([1 num_x]);
    for i=1:num_x
        for j=1:length(E)
            if E(j,1)==i || E(j,2)==i   %無向グラフ
                L_diag(i) = L_diag(i)+1;
            end
        end
    end

    Lp = -N;
    for i=1:num_x
        Lp(i,i)=L_diag(i);
    end
    
    %% xの決定
    x = sym('x',[num_x 1]);
    
    %% Gの決定
    G_sym = x(2)-x(1);
    num_G = length(G_sym);
    
%     disp('Gの準備中');
    
    %% λの定義
    num_lambda = num_G;
    lambda = sym('lambda');
    
    % λGの決定
    lG_sym = lambda.'*G_sym;
    
    %% d(λG)/dx　（手打ち）
    dlGdx_sym = sym('dlGdx_sym',[num_x 1]);
    dlGdxi = cell(num_x);
    for n=1:num_x
        dlGdx_sym(n) = diff(lG_sym, x(n));
        dlGdxi{n} = matlabFunction(dlGdx_sym(n));
    end
    dlGdx = matlabFunction(dlGdx_sym);
%% 
% end
% disp('完了');


%% パラメータ設定
A =2;
B = .1;
gamma = .01;
c = .1./L_diag;

B_p = B/sum(1./c);  %スーパバイザ用のB

day = 24;
c_delay = 5;
stp_max = day*3+1;    %s(実行step数)の最大
eps_x = .001;   %x[k]の更新の打ち切り基準:dx[k]<eps_x
dx_max = 1000;    %x[k]の更新の計算中止dx


%% シミュレーション実行
f_run = input('run?[y,n]','s');  %yで実行
if isempty(f_run)
    f_run = 'n';
end
if f_run == 'y'
    %% x,λの推移を記憶
    
    
    X = ones(num_x, stp_max);
    X_min = ones(num_x,stp_max*60);
    LAMBDA = ones(num_lambda, stp_max);  % スーパバイザ方式

    
    %% 初期条件(step = 1)
    X(:,1) = rand([num_x 1]);
    for i=1:num_x
        for mi=1:60
            X_min(:,mi) = X(:,1);
        end
    end
    
    LAMBDA(:,1) = rand(1);

    
    
    %% ステップ実行(step >=2)
%     disp('実行中...');
    for step = 2:stp_max
        
        % xの更新
        % x[0]の準備
        x = X(:, step-1);
        
        %各ノードについてのループ
        for i=1:num_x
            kx=0;
            while kx < 60
                lambda = LAMBDA(1,step-1);
                df = 2*X_min(i,(step-1)*60 + kx);
                dg = dlGdxi{i}(lambda);
                
                x(i) = x(i) -A* ( gamma*df + dg);
                kx=kx+1;
                
                X_min(i,(step-1)*60+kx) = x(i);
            end
        end
        % xの更新
        X(:,step) = x;
        
        %λの更新
        LAMBDA(1,step) = LAMBDA(1,step-1) + B_p*(X(2,step-1) - X(1,step-1)) ;
    end
end
clear f_run;


%% 結果の表示
f_plot = input('plot?[y,n]', 's');  %'y'で実行
if isempty(f_plot)
    f_plot = 'n';
end
if f_plot == 'y'
    
    
    LAMBDA_min = zeros([num_lambda stp_max*60]);
    for t=1:num_lambda
        for m=1:60
            LAMBDA_min(:,(t-1)*60+m)=LAMBDA(:,t);
        end
    end
    
    time_min = 1:stp_max*60;
    time_h = time_min./60;
    
    figure(1);
    title('xiの推移');
    plot(time_h,X_min(:,:));
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
    axis([0,20,-50,50]);
    grid on;
    
    
    figure(2);
    title('λの推移');
    plot(time_h,LAMBDA_min(:,:));
    axis([0,1.5,0,1]);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
    grid on;
    
end
clear f_plot;












