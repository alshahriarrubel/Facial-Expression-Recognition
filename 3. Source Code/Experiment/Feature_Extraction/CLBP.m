function [ result] = CLBP( img,Region_r,Region_c )
[~,~,ch]=size(img);
if(ch>1)
    img=rgb2gray(img);
end
[r,c]=size(img);
img=double(img);
result = [];

% fid = fopen('out1.txt' , 'w' );
dx = [ 0,-1,-1, -1, 0, 1, 1, 1 ];
dy = [ 1, 1, 0, -1,-1,-1, 0, 1 ];

%Matrix for storing the result of CLBP
CLBP_Sign = zeros( r, c );
CLBP_Magnitude = zeros( r, c );
CLBP_Center = zeros( r, c );
%For storing value for each direction
Local_differnce = zeros( 8,1 );
Sign_P=Local_differnce;
Magnitude_P=Local_differnce;
Magnitude_B=Local_differnce;
ratio_r = r / Region_r;
ratio_c = c / Region_c;
CLBPCode=zeros(1,256*10000+256*10+2);

Magnitude_sum=0;
Gray_level_sum=0;
for i=1:r
   for j=1:c
      Gray_level_sum=Gray_level_sum+img(i,j); 
      if(i>1 && i<r && j>1 && j<c)
          for k=1:8
              nrr = i + dx( k );
              ncc = j + dy( k ) ;
              Magnitude_sum=Magnitude_sum+abs(img(nrr,ncc)-img(i,j));                   
          end
      end
   end
end

AVG_Gray_level=Gray_level_sum/(r*c);
AVGM=Magnitude_sum/((r-2)*(c-2)*8);


%%%%%%%%%%%%%%%%%

Hist_Sign=zeros(1,256);
Hist_Magnitude=zeros(1,256);
Hist_Center=zeros(1,2);
FeatureCode=[];
%%%%%%%%%%%%%%%


%Loop for each window
for l = 1:Region_r
    %Local Window row start and end
    sr = ( ratio_r * ( l - 1 ) ) + 1;
    er = ratio_r * l ;
    for col = 1:Region_c
        %Local Window column start and end
        sc = ( ratio_c * ( col - 1 ) ) + 1;
        ec = ratio_c * col ;
        his = zeros( 256, 256, 2 );
        
        %Loops for a single local window
        Hist_Sign=zeros(1,256);
        Hist_Magnitude=zeros(1,256);
        Hist_Center=zeros(1,2);
        
        for i = sr:er
            for j = sc:ec

                %Loop for Neighbours
                for z = 1:8
                    nr = i + dx( z );
                    nc = j + dy( z ) ;
                    if( nr >= 1 && nr <= r && nc >= 1 && nc <= c ) 
                        % nr and nc for image and mr and mc for mask
                        Local_differnce( z,1 )=img(nr,nc)-img(i,j);
                        if(Local_differnce(z,1)>=0)
                            Sign_P(z,1)=1;
                        else 
                            Sign_P(z,1)=0;
                        end
                        Magnitude_P(z,1)=abs(Local_differnce( z,1 ));
                        if (Magnitude_P(z,1)>=AVGM)
                            Magnitude_B(z,1)=1;
                        else
                            Magnitude_B(z,1)=0;
                        end
                        
                    end

                end

                sum = 0;
                sum2=0;
                for z = 1:8                               
                    sum = sum + ( ( 2 .^ ( z - 1 ) )  * Sign_P( z, 1 ) ) ;
                    sum2 = sum2 + ( ( 2 .^ ( z - 1 ) )  * Magnitude_B( z, 1 ) ) ;
                end
                CLBP_Sign( i,j ) = sum ;
                CLBP_Magnitude(i,j)=sum2;
                if (img(i,j)>=AVG_Gray_level)
                   CLBP_Center(i,j)=1;
                else
                   CLBP_Center(i,j)=0;
                end
                %%%%%%%%%%%%%
%                CLBPCode(i,j)=[CLBP_Sign( i,j ) CLBP_Magnitude(i,j) CLBP_Center(i,j)];
                
                Hist_Sign(CLBP_Sign( i,j )+1)=Hist_Sign(CLBP_Sign( i,j )+1)+1;
                Hist_Magnitude(CLBP_Magnitude(i,j)+1)=Hist_Magnitude(CLBP_Magnitude(i,j)+1)+1;
                Hist_Center(CLBP_Center(i,j)+1)=Hist_Center(CLBP_Center(i,j)+1)+1;
                %%%%%%%%%%%%%%%%%
                
                %his(CLBP_Sign( i,j ) + 1,CLBP_Magnitude(i,j)+1,CLBP_Center(i,j)+1 ) = his(CLBP_Sign( i,j ) + 1,CLBP_Magnitude(i,j)+1,CLBP_Center(i,j)+1 ) + 1;
                %his(1,CLBPCode(i,j)+1)=his(1,CLBPCode(i,j)+1)+1;
                
            end
        end
        FeatureCode=[FeatureCode Hist_Sign Hist_Magnitude Hist_Center];

    end
end


result=FeatureCode;
end

