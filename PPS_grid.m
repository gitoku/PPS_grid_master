%%初期化処理
f_init = input('initialize?[y,n]','s');  %'y'の時,初期化処理実行
if isempty(f_init)
    f_init = 'n' ;
end

if f_init == 'y'
    clear all;
    disp('初期化中...');
    

    
    
%%グラフの定義###

% エッジ集合にて定義

% 最小グラフ
E = [ 1,2; 2,1];

num_x = max( max(E) )


% 隣接行列 N
    for i=1:length(E);
        N( E(i,1) , E(i,2) ) = 1;
        N( E(i,2) , E(i,1) ) = 1;  %%無向グラフの場合
    end