function [ result] = LBC( img ,Region_r, Region_c)

[~,~,ch]=size(img);
if(ch>1)
    img=rgb2gray(img);
end
[r,c]=size(img);
% img=double(img);
result = [];

% fid = fopen('out1.txt' , 'w' );
dx = [ 0,-1,-1, -1, 0, 1, 1, 1 ];
dy = [ 1, 1, 0, -1,-1,-1, 0, 1 ];

%Matrix for storing the result of LDP
LBCCode = zeros( r, c );

%For storing masked value for each mask
temp = zeros( 8,1 );

ratio_r = r / Region_r;
ratio_c = c / Region_c;

%Loop for each window
for l = 1:Region_r
    %Local Window row start and end
    sr = ( ratio_r * ( l - 1 ) ) + 1;
    er = ratio_r * l ;
    for col = 1:Region_c
        %Local Window column start and end
        sc = ( ratio_c * ( col - 1 ) ) + 1;
        ec = ratio_c * col ;
        his = zeros( 1, 9 );
        
        %Loops for a single local window
        for i = sr:er
            for j = sc:ec

                %Loop for Neighbours
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
                    sum = sum +  temp( z, 1 )  ;

                end
                LBCCode( i,j ) = sum ;
                his( 1, LBCCode( i,j ) + 1 ) = his( 1, LBCCode( i,j ) + 1 ) + 1;

            end
        end

        result = [result his];
    end
end




