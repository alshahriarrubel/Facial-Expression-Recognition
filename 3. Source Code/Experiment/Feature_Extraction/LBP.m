function [ result ] = LBP( img,localWindowr,localWindowc )

%img = imresize( img, [400 400 ]);
[r,c,ch]=size(img);
% img=double(img);
if (ch>1)
	img=rgb2gray(img);
end
result = [];

img = double(img);

% fid = fopen('out1.txt' , 'w' );
dx = [ 0,-1,-1, -1, 0, 1, 1, 1 ];
dy = [ 1, 1, 0, -1,-1,-1, 0, 1 ];

%Matrix for storing the result of LDP
lbpCode = zeros( r, c );

%For storing the 8 sorted masked values
stemp = zeros( 8,1);

%For storing masked value for each mask
temp = zeros( 8,1 );

%Mask Position
mr = 2;
mc = 2;
K = 3;
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
        
        %Loops for a single local window
        for i = sr:er
            for j = sc:ec

                %Loop for Neighbours
                temp = zeros( 8,1 );
                for z = 1:8
                    nr = i + dx( z );
                    nc = j + dy( z ) ;
                    if( nr >= 1 && nr <= r && nc >= 1 && nc <= c ) 
                        % nr and nc for image and mr and mc for mask
                        temp( z,1 )=img(nr,nc)-img(i,j);
                    end

                end

                sum = 0;
                for z = 1:8            
                    if( temp(z,1) >= 0 ) 
                        temp(z,1) = 1;
                    else
                        temp( z,1) = 0;
                    end
                    %Need to check again
                    sum = sum + ( ( 2 .^ ( z - 1 ) )  * temp( z, 1 ) ) ;

                end
                lbpCode( i,j ) = sum ;
                his( 1, lbpCode( i,j ) + 1 ) = his( 1, lbpCode( i,j ) + 1 ) + 1;

            end
        end

        result = [result his];
    end
end




