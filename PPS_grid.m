%%����������
% f_init = input('initialize?[y,n]','s');  %'y'�̎�,�������������s
% if isempty(f_init)
%     f_init = 'n' ;
% end
% 
% 
% if f_init == 'y'
    clear all;
%     disp('��������...');
    %% 
    
%%�O���t�̒�`###

% �G�b�W�W���ɂĒ�`

% �ŏ��O���t
E = [ 1,2; 2,1];

num_x = max( max(E) );


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
    
    %% G�̌���
    G_sym = [
            x(1)-x(2)
        ];
    num_G = length(G_sym);
    
%     disp('G�̏�����');
    
    %% �ɂ̒�`
    num_lambda = num_G;
    lambda = sym('lambda');
    
    % ��G�̌���
    lG_sym = lambda.'*G_sym;
    
    %% d/dx(��G)�@�i��ł��j
    dlGdx_sym = sym('dlGdx_sym',[num_x 1]);
    
    dlGdx_sym(1) = 2*x(1)+lambda;
    dlGdx_sym(2) = 2*x(2)-lambda;
    
%% 
% end
% disp('����');


%% �p�����[�^�ݒ�
A =2;
B = .1;
gamma = .01;
c = .1./L_diag;

B_p = B/sum(1./c);  %�X�[�p�o�C�U�p��B

day = 24;
c_delay = 5;
stp_max = day*3+1;    %s(���sstep��)�̍ő�
eps_x = .001;   %x[k]�̍X�V�̑ł��؂�:dx[k]<eps_x
dx_max = 1000;    %x[k]�̍X�V�̌v�Z���~dx


%% �V�~�����[�V�������s
f_run = input('run?[y,n]','s');  %y�Ŏ��s
if isempty(f_run)
    f_run = 'n';
end
if f_run == 'y'
    %% x,�ɂ̐��ڂ��L��
    
    
    X = ones(num_x, stp_max);
    X_min = ones(num_x,stp_max*60);
    LAMBDA = ones(num_lambda, stp_max);  % �X�[�p�o�C�U����

    
    %% ��������(step = 1)
    X(:,1) = rand([num_x 1]);
    for i=1:num_x
        for mi=1:60
            X_min(:,mi) = X(:,1);
        end
    end
    
    LAMBDA(:,1) = rand(1);

    
    
    %% �X�e�b�v���s(step >=2)
%     disp('���s��...');
    for step = 2:stp_max
        
        % x�̍X�V
        % x[0]�̏���
        x = X(:, step-1);
        
        %�e�m�[�h�ɂ��Ẵ��[�v
        for i=1:num_x
            kx=0;
            while kx < 60
                
%                 x(i) = x(i) - A* (
                
                kx=kx+1;
                
                X_min(i,(step-1)*60+kx) = x(i);
            end
        end
    end
end

















