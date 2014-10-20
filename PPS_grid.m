%% ����������
f_init = input('initialize?[y,n]','s');  %'y'�̎�,�������������s
if isempty(f_init)
    f_init = 'n' ;
end


if f_init == 'y' || exist('define.mat','file')~=2
    clear all;
    disp('��������...');
      
    %% �O���t�̒�`###
 
    % �ŏ��O���t
    E = [   
        1,4;
        2,5;
        3,6;
    ];

    % �G�[�W�F���g�̐�, �ϐ��̐�����
    time_agt=3; 
    num_x = max( max(E) );
    agt_num = num_x/time_agt; 


    % �אڍs�� N
    N = zeros([num_x,num_x]);
    for i=1:length(E);
        N( E(i,1) , E(i,2) ) = 1;
        N( E(i,2) , E(i,1) ) = 1;  %%�����O���t�̏ꍇ
    end
    
    %�O���t���v���V�A���@L
    L_diag=zeros([1 num_x]);
    for i=1:num_x
        for j=1:length(E)
            if E(j,1)==i || E(j,2)==i   %�����O���t
                L_diag(i) = L_diag(i)+1;
            end
        end
    end

    Lp = -N;
    for i=1:num_x
        Lp(i,i)=L_diag(i);
    end
    
        
    %% x�̌���
    x = sym('x',[num_x 1]);
    
    %�G�[�W�F���g�^�C�v
    %1:������
    %2:���v��
    
    agt_type=[
        2,2,2,...%���v��1
        1,1,1,...%������1
    ];
    
    agt_sym = sym('x',[agt_num time_agt]);
    
       
    %% G�̌���
    G_hat_sym = [
        x(4)-x(1);
        x(5)-x(2);
        x(6)-x(3);
    ];

    TC = 10  %�e���Ԃ̐���
    G_sym = [G_hat_sym.' (-G_hat_sym).'].';
    for i=1:num_x
        G_sym(i) = G_sym(i) + TC;
    end
    
    num_G = length(G_sym);
    
    G = cell(num_G);
    for m=1:num_G
        G{m} = matlabFunction(G_sym(m),'vars',{x});
    end
    
  
    %% H�̌���
    % ������v��, ������������
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
%     disp('G�̏�����');
    
   
    %% �ɂ̒�`
    G_lambda = sym('G_lambda',[num_G 1]);
    lG_sym = 0;
    % ��G�̌���
    for i=1:num_G
        lG_sym = lG_sym + G_lambda(i)*G_sym(i);
    end
    lG_sym
    
    % ��H�̌���
    H_lambda = sym('H_lambda');
    lH_sym = H_lambda.'*H_sym
    
   
    %% d(��G)/dx �̌���
    dlGdx_sym = sym('dlGdx_sym',[num_x 1]);
    dlGdxi = cell(num_x,1);
    for n=1:num_x
        dlGdx_sym(n) = diff(lG_sym, x(n));
        dlGdxi{n} = matlabFunction(dlGdx_sym(n),'vars',{G_lambda});
    end
    dlGdx = matlabFunction(dlGdx_sym,'vars',{G_lambda});

    
    %% d(��H)/dx�@�̌���
    dlHdx_sym = sym('dlHdx_sym',[num_x 1]);
    dlHdxi = cell(num_x,1);
    for n=1:num_x
        dlHdx_sym(n) = diff(lH_sym, x(n));
        dlHdxi{n} = matlabFunction(dlHdx_sym(n),'vars',{H_lambda});
    end
    dlHdx = matlabFunction(dlHdx_sym,'vars',{H_lambda});
        
    
    
    save('define','time_agt','num_x','agt_num','N','L_diag','Lp',...
        'agt_type','G_hat_sym','G_sym','G','num_G','H_sym','H','num_H','G_lambda','H_lambda','lG_sym','lH_sym',...
        'dlGdx_sym','dlGdxi','dlHdx_sym','dlHdxi');
    clear all;
    disp('����������')
       
end
clear f_init;
clear all;

load('define');

%% �p�����[�^�ݒ�
A =2;
B = .1;
gamma = .01;
c = .1./L_diag;

B_h = B/sum(1./c);  %�X�[�p�o�C�U�p��B
B_g = .001;
stp_max = 70;    %s(���sstep��)�̍ő�
opt_max = 100;
eps_x = .001;   %x[k]�̍X�V�̑ł��؂�:dx[k]<eps_x
x_max = 1000;    %x[k]�̍X�V�̌v�Z���~dx


%% �V�~�����[�V�������s
f_run = input('run?[y,n]','s');  %y�Ŏ��s
if isempty(f_run)
    f_run = 'n';
end

if f_run == 'y'
    %% x,�ɂ̐��ڂ��L��
    
    X = ones(num_x, stp_max);
    X_loop = ones(num_x,opt_max);
    LAMBDA_G = zeros(num_G, stp_max);  % �X�[�p�o�C�U����
    LAMBDA_H = zeros(num_H,stp_max);
    
    %% ��������(step = 1)
    rng(100); % �����Œ�p
    X(:,1) = rand([num_x 1]);
    X_loop(:,1,1) = X(:,1);
    first_x=X_loop(:,1)
    
    for i=1:num_G
        LAMBDA_G(i,1) = rand(1);
    end
    LAMBDA_H(:,1) = rand(1); %�ɂ̏����l
    
    d = rand([num_x 1]); % xi�̏��]��
    
%     for i=1:num_x
%         if agt_type(i)==2
%             d(i)=d(1);
%         elseif agt_type(i)==1;
%             d(i)=d(4);
%         end
%     end
    d
    
    disp('x,�ɍX�V��')
    %% �X�e�b�v���s(step >=2)
%     disp('���s��...');
    for step = 2:stp_max
        
        % x[0]�̏���
        x = X(:, step-1);
        
        %�e�m�[�h�ɂ��Ẵ��[�v
        for i=1:num_x
            k=2;
            
            dg = dlGdxi{i}(LAMBDA_G(:,step-1));
            dh = dlHdxi{i}(LAMBDA_H(:,step-1));
            while true
                df = 2*gamma*(X_loop(i,k-1)-d(i));
                % x�̍X�V                                     
                x(i) = x(i) -A* ( df + dg + dh);
                
                X_loop(i,k)=x(i);
                dx = X_loop(i,k)-X_loop(i,k-1);
                if abs(dx) < eps_x
                    for n=k:opt_max
                        X_loop(i,n)=X_loop(i,k);
                    end
                    break;
                elseif abs(dx)>x_max
                    disp('���U');
                    break;
                end

                k=k+1;
                
            end
        end
        X(:,step) = x;
        %�ɂ̍X�V:�X�[�p�o�C�U����
        for m=1:num_G
            LAMBDA_G(m,step) =max(0, LAMBDA_G(m,step-1) + B_g*G{m}(X(:,step)));
        end
        
        for m=1:num_H
%             LAMBDA_H(m,step) =(max(0, LAMBDA_H(m,step-1) + B_p*H{m}(X(:,step))));
            LAMBDA_H(m,step) =LAMBDA_H(m,step-1) + B_h*H{m}(X(:,step));
        end
        
    end
    last_X = X(:,step)
    last_LAMBDA_G = LAMBDA_G(1:3,step) -LAMBDA_G(4:6,step)
    last_LAMBDA_H = LAMBDA_H(:,step)
    
    save ('result','stp_max','opt_max','d','X','last_X','LAMBDA_G',...
        'LAMBDA_H','last_LAMBDA_G','last_LAMBDA_H','gamma');
    clear all;
end
clear f_run;

%% ���ʂ̕\��
f_plot = input('plot?[y,n]', 's');  %'y'�Ŏ��s
if isempty(f_plot)
    f_plot = 'n';
end
if f_plot == 'y'
    
    load('define');
    load('result');
    
    FX = zeros([1 stp_max]);
    HX = zeros([num_H stp_max]);
    GX = zeros([num_G stp_max]);
    for step=1:stp_max
        for i=1:num_x
            FX(:,step)= FX(:,step) + gamma*(X(i,step)-d(i))^2;
        end
        for m=1:num_H
            HX(m,step) = H{m}(X(:,step));
        end
        for n=1:num_G
            GX(n,step) = G{n}(X(:,step));
        end
    end
    last_GX = GX(:,step)
    last_HX = HX(:,step)
    last_FX = FX(:,step)
    
    time=0:stp_max-1;
    time2=1:stp_max;
    opt_time=opt_max;
    
    % x(i)�̍œK��
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
    title('x(i)�̍œK������');
    plot(time,X(1,:),time+10,X(2,:),time+20,X(3,:),time,X(4,:),time+10,X(5,:),time+20,X(6,:),'LineWidth',1.5);
%     plot(time,X(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     xlim([0 20]);
    xlabel('step');
    ylabel('x(i)');
    grid on;
   
%     figure(2);
%     title('x(i)�̍œK��')
%     stem(x_Time, plot_X,'fill');
% %     ylim([0 100]);
%     xlabel('����')
%     ylabel('�ϐ���')
%     grid on;
    
    figure(3);
    title('H�̐���');
    plot(time,HX(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     axis([0 40 -700 100]);
    xlabel('step');
    ylabel('H(x)');
    grid on;
    
    figure(4);
    list = [1 4];
    title('G�̐���');
%     plot(time,GX(1,:),time+10,GX(2,:),time+20,GX(3,:),time,GX(4,:),time+10,GX(5,:),time+20,GX(6,:),'LineWidth',1.5);
    plot(time,GX(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     axis([0 80 -100 100]);
    xlabel('step');
    ylabel('G(x)');
    grid on;
    
    LAMBDA_plot = LAMBDA_G(1:3,:) - LAMBDA_G(4:6,:);
    figure(5)
    title('G�̃ɂ̐���');
    plot(time,LAMBDA_plot(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     xlim([0 40]);
    xlabel('step');
    ylabel('lambda_G');
    grid on;
    
    figure(6);
    title('H�̃ɂ̐���');
    plot(time,LAMBDA_H(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     xlim([0 40]);
    xlabel('step');
    ylabel('lambda_H');
    grid on;
    
    figure(7);
    title('F(x)�̐���');
    plot(time,FX(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
%     xlim([0 40]);
    xlabel('step');
    ylabel('F(x)');
    grid on;
      
    
    
end
clear f_plot;












