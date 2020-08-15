function [ result] = LDP( img,region_x,region_y )

[r,c]=size(img);
% img=double(img);
%img=rgb2gray(img);
result = [];
%[ r c k ] = size( img );
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
ldpCode = zeros( r, c );

%For storing the 8 sorted masked values
stemp = zeros( 8,1);

%For storing masked value for each mask
temp = zeros( 8,1 );
temp2 = zeros( 8,1 );
%Mask Position
mr = 2;
mc = 2;
K = 3;

ratio_r = r / region_x;
ratio_c = c / region_y;

%Loop for each window
for l = 1:region_x
    %Local Window row start and end
    sr = ( ratio_r * ( l - 1 ) ) + 1;
    er = ratio_r * l ;
    for col = 1:region_y
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

                s = 0;
                temp2=abs(temp);
                [max1,max1_idx]=max(temp2);
                temp2(max1_idx)=-inf;
                [max2,max2_idx]=max(temp2);
                temp2(max2_idx)=-inf;
                [max3,max3_idx]=max(temp2);
                temp2(max3_idx)=-inf;
                for z = 1:8
                    
                    if( abs(temp(z,1))==max1 || abs(temp(z,1))==max2 || abs(temp(z,1))==max3) 
                        temp(z,1) = 1;
                    else
                        temp( z,1) = 0;
                    end
                    %Convert to decimal
                    s = s + ( ( 2 .^ ( z - 1 ) )  * temp( z, 1 ) ) ;

                end  
                
                ldpCode( i,j ) = s ;
                his( 1, ldpCode( i,j ) + 1 ) = his( 1, ldpCode( i,j ) + 1 ) + 1;

            end
        end

        result = [result his];
    end
end




