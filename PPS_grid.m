%%����������
f_init = input('initialize?[y,n]','s');  %'y'�̎�,�������������s
if isempty(f_init)
    f_init = 'n' ;
end

if f_init == 'y'
    clear all;
    disp('��������...');
    

    
    
%%�O���t�̒�`###

% �G�b�W�W���ɂĒ�`

% �ŏ��O���t
E = [ 1,2; 2,1];

num_x = max( max(E) )


% �אڍs�� N
    for i=1:length(E);
        N( E(i,1) , E(i,2) ) = 1;
        N( E(i,2) , E(i,1) ) = 1;  %%�����O���t�̏ꍇ
    end