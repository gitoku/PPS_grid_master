%% ����������
f_init = input('initialize?[y,n]','s');  %'y'�̎�,�������������s
if isempty(f_init)
    f_init = 'n' ;
end


if f_init == 'y'
    clear all;
    disp('��������...');
    
%% �O���t�̒�`###

% �G�b�W�W���ɂĒ�`

% �ŏ��O���t
E = [   
    1,2;2,1;
];

num_x = max( max(E) );
num_agt = num_x/3;

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
        2,1,...%Region1
    ];
    
    %% G�̌���
    % ������v��, ������������
    G_sym = -sym('x',1);
    for i=2:num_x
        if agt_type(i)==2
            G_sym = G_sym - x(i);
        elseif agt_type(i)==1
            G_sym = G_sym + x(i);
        end
    end
    num_G = length(G_sym);
    G_sym
    
    G = cell(num_G);
    for m=1:num_G
        G{m} = matlabFunction(G_sym(m),'vars',{x});
    end
%     disp('G�̏�����');
    
    %% �ɂ̒�`
    num_lambda = num_G;
    lambda = sym('lambda');
    
    % ��G�̌���
    lG_sym = lambda.'*G_sym
    
    %% d(��G)/dx�@�i��ł��j
    dlGdx_sym = sym('dlGdx_sym',[num_x 1]);
    dlGdxi = cell(num_x,1);
    for n=1:num_x
        dlGdx_sym(n) = diff(lG_sym, x(n));
        dlGdxi{n} = matlabFunction(dlGdx_sym(n));
    end
    dlGdx = matlabFunction(dlGdx_sym);
    
end

%% �p�����[�^�ݒ�
A =2;
B = .1;
gamma = .01;
c = .1./L_diag;

B_p = B/sum(1./c);  %�X�[�p�o�C�U�p��B

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
    LAMBDA = ones(num_lambda, stp_max);  % �X�[�p�o�C�U����

    
    %% ��������(step = 1)
    
    X(:,1) =rand([num_x 1]);
    X_loop(:,1,1) = X(:,1);
    first_x=X_loop(:,1)
    
    LAMBDA(:,1) = rand(1); %�ɂ̏����l
    d = rand([num_x 1]) % xi�̏��]��
    
    %% �X�e�b�v���s(step >=2)
%     disp('���s��...');
    for step = 2:stp_max
        
        % x[0]�̏���
        x = X(:, step-1);
        
        %�e�m�[�h�ɂ��Ẵ��[�v
        for i=1:num_x
            k=2;
            while true
                % x�̍X�V
                df = 2*gamma*(X_loop(i,k-1)-d(i));
                lambda = LAMBDA(1,step-1);
                dg = dlGdxi{i}(lambda);
                x(i) = x(i) -A* ( df + dg);
                
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
        
        %�ɂ̍X�V
        for m = 1:num_lambda
            LAMBDA(m,step) =(max(0, LAMBDA(m,step-1) + B_p*G{m}(X(:,step))));
        end
        
    end
    last_X = X(:,step)
    
    save result;
end
clear f_run;


%% ���ʂ̕\��
f_plot = input('plot?[y,n]', 's');  %'y'�Ŏ��s
if isempty(f_plot)
    f_plot = 'n';
end
if f_plot == 'y'
    
    FX = zeros([stp_max 1]);
    GX = zeros([num_lambda stp_max]);
    for step=1:stp_max
        for i=1:num_x
            FX(step)= gamma*(X(i,step)-d(i))^2;
        end
        for m=1:num_lambda
            GX(m,step) = G{m}(X(:,step));
        end
    end
    
    time=0:stp_max-1;
    opt_time=opt_max;
   
    
    
    figure(1);
    title('x�̐���');
    plot(time,X(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
    xlim([1 30]);
    xlabel('step');
    ylabel('x(i)');
    grid on;
    
    figure(2);
    title('G�̐���');
    plot(time,GX(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
    xlim([0 70]);
    xlabel('step');
    ylabel('G(x)');
    grid on;
    
    figure(3);
    title('�ɂ̐���');
    plot(time,LAMBDA(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
    xlim([0 70]);
    xlabel('step');
    ylabel('lambda');
    grid on;
    
    figure(4);
    title('F(x)�̐���');
    plot(time,FX(:,:),'LineWidth',1.5);
    set(gca,'FontName','Times','Fontsize',18,'LineWidth',1.5);
    xlim([0 30]);
    xlabel('step');
    ylabel('F(x)');
    grid on;
   
end
clear f_plot;












