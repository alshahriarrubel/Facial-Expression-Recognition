function [ result,ldnCode ] = LDN( img )
img=rgb2gray(img);
img = imresize( img, [400 400 ]);
[ r c ] = size( img );
% img=double(img);

result = [];
%
kmask = zeros( 3,3,8 );
kmask( :,:,1 ) = [ -3 -3 5 ; -3 0 5; -3 -3 5 ];
kmask( :,:,2 ) = [ -3 5 5 ; -3 0 5; -3 -3 -3 ];
kmask( :,:,3 ) = [ 5 5 5 ; -3 0 -3; -3 -3 -3 ];
kmask( :,:,4 ) = [ 5 5 -3 ; 5 0 -3; -3 -3 -3 ];
kmask( :,:,5 ) = [ 5 -3 -3 ; 5 0 -3; 5 -3 -3 ];
kmask( :,:,6 ) = [ -3 -3 -3 ; 5 0 -3; 5 5 -3 ];
kmask( :,:,7 ) = [ -3 -3 -3 ; -3 0 -3; 5 5 5 ];
kmask( :,:,8 ) = [ -3 -3 -3 ; -3 0 5; -3 5 5 ];


% var2 = rgb2gray( img ) ;
var2 = double(img);

% fid = fopen('out1.txt' , 'w' );
dx = [ 0, -1, -1, -1, 0, 1, 1, 1 ];
dy = [ -1, -1, 0, 1, 1, 1, 0, -1 ];

%Matrix for storing the result of LDP
ldnCode = zeros( r, c );

%For storing the 8 sorted masked values
stemp = zeros( 8,1);

%For storing masked value for each mask
temp = zeros( 8,1 );

%Mask Position
mr = 2;
mc = 2;
K = 3;
localWindowr =10;
localWindowc=10;
ratio_r = r / localWindowr;
ratio_c = c / localWindowc;

%Loop for each window
for l = 1:localWindowr
    %Local Window row start and end
    sr = ( ratio_r * ( l - 1 ) ) + 1;
    er = ratio_r * l ;
    for col = 1:localWindowc
        %Local Window column start and end
        sc = ( ratio_c * ( col - 1 ) ) + 1;
        ec = ratio_c * col ;
        his = zeros( 1, 256 );
        cnt = 0;
        %Loops for a single local window
        for i = sr:er
            for j = sc:ec
                %Loop for selecting Kirsch Masks
                for k = 1:8
                    sum = 0;
                    %Loop for direction in each Mask
                    for z = 1:8
                        nr = i + dx( z );
                        nc = j + dy( z ) ;
                        if( nr >= 1 && nr <= r && nc >= 1 && nc <= c ) 
                            % nr and nc for image and mr and mc for mask
                            sum = sum + ( var2( nr,nc ) * kmask( mr + dx( z ) , mc + dy( z ) , k ) ) ;
                        end

                    end
                    temp( k,1 ) = sum;

                end
                %% Find Max_Edge_Response_index
                [Max_value,Max_value_index]=max(temp(:));
                [Max_value_index_row, Max_value_index_col] = ind2sub(size(temp),Max_value_index);
                
                %%%% Find Min_Edge_Response_index
                [Min_value,Min_value_index]=min(temp(:));
                [Min_value_index_row, Min_value_index_col] = ind2sub(size(temp),Min_value_index);
                
                % Generate LDN code
                ldnCode( i,j ) = 8*(Max_value_index_row-1)+(Min_value_index_row-1) ;

                his( 1, ldnCode( i,j ) + 1 ) = his( 1,ldnCode( i,j ) + 1 ) + 1;

            end
        end

        result = [result his];
    end
end



