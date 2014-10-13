%% 初期化処理
f_init = input('initialize?[y,n]','s');  %'y'の時,初期化処理実行
if isempty(f_init)
    f_init = 'n' ;
end


if f_init == 'y' || exist('define.mat','file')~=2
    clear all;
    disp('初期化中...');
      
    %% グラフの定義###
 
    % 最小グラフ
    E = [   
        1,4;
        2,5;
        3,6;
    ];

    % エージェントの数, 変数の数決定
    time_agt=3; 
    num_x = max( max(E) );
    agt_num = num_x/time_agt; 


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
    
    %エージェントタイプ
    %1:供給家
    %2:需要家
    
    agt_type=[
        2,2,2,...%需要家1
        1,1,1,...%供給家1
    ];
    
    agt_sym = sym('x',[agt_num time_agt]);
    
       
    %% Gの決定
    G_hat_sym = [
        x(4)-x(1);
        x(5)-x(2);
        x(6)-x(3);
    ];

    TIME_CONSTRAINT= 20; %各時間の制約
    G_sym = [G_hat_sym.' (-G_hat_sym).'].';
    for i=1:num_x
        G_sym(i) = G_sym(i) - TIME_CONSTRAINT;
    end
    
    num_G = length(G_sym);
    
    G = cell(num_G);
    for m=1:num_G
        G{m} = matlabFunction(G_sym(m),'vars',{x});
    end
    
  
    %% Hの決定
    % 奇数が需要家, 偶数が供給家
    H_sym = -sym('x',1);
    for i=2:num_x
        if agt_type(i)==2
            H_sym = H_sym - x(i);
        elseif agt_type(i)==1
            H_sym = H_sym + x(i);
        end
    end
    num_H = length(H_sym);
    
    H = cell(num_H);
    for m=1:num_H
        H{m} = matlabFunction(H_sym(m),'vars',{x});
    end
%     disp('Gの準備中');
    
   
    %% λの定義
    num_lambda = num_H + num_G;
    lambda = sym('lambda',[num_lambda 1]);
    lG_sym = sym('lG_sym',[6 1]);
    % λGの決定
    for i=1:num_G
        lG_sym(i) = lambda(i)*G_sym(i);
    end
    lG_sym
    
    % λHの決定
    lH_sym = lambda(num_lambda).'*H_sym
    
   
    %% d(λG)/dx の決定
    dlGdx_sym = sym('dlGdx_sym',[num_x num_x]);
    dlGdxi = cell(num_x,1);
    for n=1:num_x
        for m=1:num_x
        dlGdx_sym(n,m) = diff(lG_sym(m), x(n));
        dlGdxi{n,m} = matlabFunction(dlGdx_sym(n,m));
        end
    end
    dlGdx = matlabFunction(dlGdx_sym);
    
    % Region定義
    lambda_matrix = zeros([num_x num_lambda]);
    for i = 1:num_x
        for j=1:num_lambda
        end
    end
    %% d(λH)/dx　の決定
    dlHdx_sym = sym('dlHdx_sym',[num_x 1]);
    dlHdxi = cell(num_x,1);
    for n=1:num_x
        dlHdx_sym(n) = diff(lH_sym, x(n));
        dlHdxi{n} = matlabFunction(dlHdx_sym(n));
    end
    dlHdx = matlabFunction(dlHdx_sym);
        
    save('define','time_agt','num_x','agt_num','N','L_diag','Lp',...
        'agt_type','G_hat_sym','G_sym','G','H_sym','H','num_lambda','lambda','lG_sym','lH_sym',...
        'dlGdx_sym','dlGdxi','dlHdx_sym','dlHdxi');
    clear all;
    disp('初期化完了')
       
end
clear f_init;
clear all;

load('define');

%% パラメータ設定
A =2;
B = .1;
gamma = .01;
c = .1./L_diag;

B_p = B/sum(1./c);  %スーパバイザ用のB

stp_max = 50;    %s(実行step数)の最大
opt_max = 100;
eps_x = .001;   %x[k]の更新の打ち切り基準:dx[k]<eps_x
x_max = 1000;    %x[k]の更新の計算中止dx


%% シミュレーション実行
f_run = input('run?[y,n]','s');  %yで実行
if isempty(f_run)
    f_run = 'n';
end

if f_run == 'y'
    %% x,λの推移を記憶
    
    X = ones(num_x, stp_max);
    X_loop = ones(num_x,opt_max);
    LAMBDA = ones(num_lambda, stp_max);  % スーパバイザ方式

    
    %% 初期条件(step = 1)
    
    X(:,1) =100*rand([num_x 1]);
    X_loop(:,1,1) = X(:,1);
    first_x=X_loop(:,1)
    
    LAMBDA(:,1) = 10*rand(1); %λの初期値
    d = 100*rand([num_x 1]); % xiの所望量
    
%     for i=1:num_x
%         if agt_type(i)==2
%             d(i)=d(1);
%         elseif agt_type(i)==1;
%             d(i)=d(4);
%         end
%     end
    d
    
    
    %% ステップ実行(step >=2)
%     disp('実行中...');
    for step = 2:stp_max
        
        % x[0]の準備
        x = X(:, step-1);
        
        %各ノードについてのループ
        for i=1:num_x
            k=2;
            dg = 0;
            while true
                % xの更新
                df = 2*gamma*(X_loop(i,k-1)-d(i));
                lambda = LAMBDA(1,step-1);
                
                dh = dlHdxi{i}(lambda);
                x(i) = x(i) -A* ( df + dg + dh);

                X_loop(i,k)=x(i);
                dx = X_loop(i,k)-X_loop(i,k-1);
                if abs(dx) < eps_x
                    for n=k:opt_max
                        X_loop(i,n)=X_loop(i,k);
                    end
                    break;
                elseif abs(dx)>x_max
                    disp('発散');
                    break;
                end

                k=k+1;
                disp([i j k])
            end
        end
        X(:,step) = x;
        
        %λの更新:スーパバイザ方式
        for m = 1:num_lambda
            LAMBDA(m,step) =(max(0, LAMBDA(m,step-1) + B_p*H{m}(X(:,step))));
        end
        
    end
    last_X = X(:,step)
    last_LAMBDA = LAMBDA(:,step)
    
    save ('result','stp_max','opt_max','d','X','last_X','LAMBDA','last_LAMBDA','gamma');
    clear all;
end
clear f_run;

%% 結果の表示
f_plot = input('plot?[y,n]', 's');  %'y'で実行
if isempty(f_plot)
    f_plot = 'n';
end
if f_plot == 'y'
    
    load('define');
    load('result');
    
    FX = zeros([stp_max 1]);
    HX = zeros([num_lambda stp_max]);
    for step=1:stp_max
        for i=1:num_x
            FX(step)= gamma*(X(i,step)-d(i))^2;
        end
        for m=1:num_lambda
            HX(m,step) = H{m}(X(:,step));
        end
    end
    
    time=0:stp_max-1;
    time2=1:stp_max;
    opt_time=opt_max;
    
    % x(i)の最適解
    x_Time = [0; 10; 20]; 
    plot_X = zeros(time_agt, agt_num);
    i = 0;
    for j = 1:agt_num
        for i = 1:time_agt
            if j ==1;
                plot_X(i,j) = last_X(i);
            else
                plot_X(i,j) = last_X(i+time_agt);
            end
        end
    end
    
    figure(1);
    title('x(i)の最適化推移');
    plot(time,X(1,:),time+10,X(2,:),time+20,X(3,:),time,X(4,:),time+10,X(5,:),time+20,X(6,:),'LineWidth',1.5);
%     plot(time,X(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     xlim([0 20]);
    xlabel('step');
    ylabel('x(i)');
    grid on;
   
    figure(5);
    title('x(i)の最適解')
    stem(x_Time, plot_X,'fill');
%     ylim([0 100]);
    xlabel('時間')
    ylabel('変数量')
    grid on;
    
    
%     figure(1);
%     title('需要家1のx');
%     plot(time,Dp(1,:),time,d(1),'LineWidth',1.5);
%     set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);

    
    figure(2);
    title('Hの推移');
    plot(time,HX(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     axis([0 40 -700 100]);
    xlabel('step');
    ylabel('G(x)');
    grid on;
    
    figure(3);
    title('λの推移');
    plot(time,LAMBDA(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     xlim([0 40]);
    xlabel('step');
    ylabel('lambda');
    grid on;
    
    figure(4);
    title('F(x)の推移');
    plot(time,FX(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     xlim([0 40]);
    xlabel('step');
    ylabel('F(x)');
    grid on;
      
    
    
end
clear f_plot;












