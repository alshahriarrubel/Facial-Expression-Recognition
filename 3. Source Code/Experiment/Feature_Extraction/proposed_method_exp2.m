function [result] = proposed_method_exp2(img,left_eyebrow,right_eyebrow,lower_eye,upper_lip,Region_r,Region_c)
[row,col,channel]=size(img);
if(channel>2)
    img=rgb2gray(img);
end
img=double(img);
%result = [];

kmask = zeros( 3,3,8 );
kmask( :,:,1 ) = [ -3 -3 5 ; -3 0 5; -3 -3 5 ];
kmask( :,:,2 ) = [ -3 5 5 ; -3 0 5; -3 -3 -3 ];
kmask( :,:,3 ) = [ 5 5 5 ; -3 0 -3; -3 -3 -3 ];
kmask( :,:,4 ) = [ 5 5 -3 ; 5 0 -3; -3 -3 -3 ];
kmask( :,:,5 ) = [ 5 -3 -3 ; 5 0 -3; 5 -3 -3 ];
kmask( :,:,6 ) = [ -3 -3 -3 ; 5 0 -3; 5 5 -3 ];
kmask( :,:,7 ) = [ -3 -3 -3 ; -3 0 -3; 5 5 5 ];
kmask( :,:,8 ) = [ -3 -3 -3 ; -3 0 5; -3 5 5 ];

dx = [ 0, -1, -1, -1, 0, 1, 1, 1 ];
dy = [ -1, -1, 0, 1, 1, 1, 0, -1 ];

Edge_Responses = zeros( row,col,8);
Original_Edge_Responses = zeros( row,col,8);
mr = 2;
mc = 2;

Primary_Edge_Response=zeros(row,col);
Secondary_Edge_Responses=zeros(row,col);
Primary_direction=zeros(row,col);
Secondary_direction=zeros(row,col);

hist=zeros(1,3825);




for i = 1:row
    for j = 1:col
        %Loop for selecting Kirsch Masks
        for k = 1:8
            sum = 0;
            %Loop for direction in each Mask
            for z = 1:8
                nr = i + dx( z );
                nc = j + dy( z ) ;
                if( nr >= 1 && nr <= row && nc >= 1 && nc <= col ) 
                    % nr and nc for image and mr and mc for mask
                    sum = sum + ( img( nr,nc ) * kmask( mr + dx( z ) , mc + dy( z ) , k ) ) ;
                end
            end
            Edge_Responses(i,j,k) = sum; 
            Original_Edge_Responses(i,j,k)=sum;
        end
        [Primary_Edge_Response(i,j),Primary_direction(i,j)]=max(abs(Edge_Responses(i,j,:)));
        Edge_Responses(i,j,Primary_direction(i,j))=0;
        [Secondary_Edge_Responses(i,j),Secondary_direction(i,j)]=max(abs(Edge_Responses(i,j,:)));    
        
        front_neighbour_of_Primary=Primary_direction(i,j)+1;
        back_neighbour_of_Primary=Primary_direction(i,j)-1;
        if(front_neighbour_of_Primary>8)
            front_neighbour_of_Primary=1;
        end
        if(back_neighbour_of_Primary<1)
            back_neighbour_of_Primary=8;
        end

        if(Secondary_direction(i,j)==front_neighbour_of_Primary || Secondary_direction(i,j)==back_neighbour_of_Primary)
            Edge_Responses(i,j,front_neighbour_of_Primary)=0;
            Edge_Responses(i,j,back_neighbour_of_Primary)=0;
            [Secondary_Edge_Responses(i,j),Secondary_direction(i,j)]=max(Edge_Responses(i,j,:));
        end
        
    end
end


%% get the threshold sigma & average pixel value
sum_pixel=0;
sum_edge_response=0;
min_edge_res_cheek=9999999;
max_edge_res_cheek=0;
for ii=left_eyebrow:right_eyebrow
    for jj=lower_eye:upper_lip
        primary_edge_res=Primary_Edge_Response(ii,jj);
        hist(1,primary_edge_res+1)=hist(1,primary_edge_res+1)+1;
        sum_edge_response=sum_edge_response+primary_edge_res; 
        if(primary_edge_res>max_edge_res_cheek)
            max_edge_res_cheek=primary_edge_res;
        end
        sum_pixel=sum_pixel+img(ii,jj);
        
    end
end
%[~,Sigma]=max(hist);
%Sigma=round(sum_edge_response/((right_eyebrow-left_eyebrow+1)*(upper_lip-lower_eye+1)));
Avg_edge_local=round(sum_edge_response/((right_eyebrow-left_eyebrow+1)*(upper_lip-lower_eye+1)));
Avg_edge_global=mean(Primary_Edge_Response);
%%%%%%%%%%%%%%%%%%%%%%%
%% Date: 04/10/18
Min_Primary_Edge_Response=min(reshape(Primary_Edge_Response(left_eyebrow:right_eyebrow,lower_eye:upper_lip),1,[]));
Median_Primary_Edge_Response=median(reshape(Primary_Edge_Response(left_eyebrow:right_eyebrow,lower_eye:upper_lip),1,[]));
%fprintf('Median: %f\n',Median_Primary_Edge_Response);
%%%%%%%%%%%%%%%%%%%%%%%
%Sigma=max_edge_res_cheek;
Sigma=Median_Primary_Edge_Response;
Average_pixel=sum_pixel/((right_eyebrow-left_eyebrow+1)*(upper_lip-lower_eye+1));
%fprintf('Average_pixel: %f\n',Average_pixel);


%% Code for Edge Response
Sign_Code=zeros(row,col);
Edge_code=zeros(row,col);
Pixel_code=zeros(row,col);
Direction_code=zeros(row,col);
FeatureCode=[];

ratio_r=row/Region_r;
ratio_c=col/Region_c;

for rr=1:Region_r
    sr=( ratio_r * ( rr - 1 ) ) + 1;
    er = ratio_r * rr ;
    for cc=1:Region_c
        sc = ( ratio_c * ( cc - 1 ) ) + 1;
        ec = ratio_c * cc ;
        Hist_Sign_code=zeros(1,256);
        Hist_Edge_code=zeros(1,256);
        Hist_Pixel_code=zeros(1,2);
        Hist_Direction_code=zeros(1,88);
        for i = sr:er
            for j = sc:ec

                %% Sign Code
                Sign_Code(i,j)=0;
                for k=1:8
                  if(Original_Edge_Responses(i,j,k)>=0)
                      bit=1;
                  else
                      bit=0;             
                  end
                  Sign_Code(i,j)=Sign_Code(i,j)+bit*2^(k-1);           
               end       
               Hist_Sign_code(1,Sign_Code(i,j)+1)=Hist_Sign_code(1,Sign_Code(i,j)+1)+1;
                
                
                
                
               %% Edge Code
               Edge_code(i,j)=0;
               for k=1:8
                  if(abs(Original_Edge_Responses(i,j,k)>=Sigma))
                      bit=1;
                  else
                      bit=0;             
                  end
                  Edge_code(i,j)=Edge_code(i,j)+bit*2^(k-1);           
               end       
               Hist_Edge_code(1,Edge_code(i,j)+1)=Hist_Edge_code(1,Edge_code(i,j)+1)+1;

               %% Pixel Code
               if(img(i,j)>=Average_pixel)
                   Pixel_code(i,j)=1;
               else 
                   Pixel_code(i,j)=0;
               end
               Hist_Pixel_code(1,Pixel_code(i,j)+1)=Hist_Pixel_code(1,Pixel_code(i,j)+1)+1;

               %% Direction Code
               Direction_code(i,j)=10*Primary_direction(i,j)+Secondary_direction(i,j);
               Hist_Direction_code(1,Direction_code(i,j)+1)=Hist_Direction_code(1,Direction_code(i,j)+1)+1;
           end
        end
        FeatureCode=[FeatureCode Hist_Sign_code Hist_Edge_code Hist_Direction_code Hist_Pixel_code];
       
    end
    
end

%% Histogram for each of the 10 blocks 
% size of each histogram is 256 as 2^8=256
% r=row;
% c=col;

%% Histogram for Edge Code

% [hist1,edges] = histcounts(Edge_code,256); 
% [hist2 edges] = histcounts(Edge_code(1:round(r/2),1:round(c/2)),256);
% [hist3,edges] = histcounts(Edge_code(1:round(r/2),round(c/2+1):c),256);
% [hist4,edges] = histcounts(Edge_code(round(r/2+1):r,1:round(c/2)),256);
% [hist5,edges] = histcounts(Edge_code(round(r/2+1):r,round(c/2+1):c),256);
% [hist6,edges] = histcounts(Edge_code(round(r/4+1):round(3*r/4),round(c/4+1):round(3*c/4)),256);
% [hist7,edges] = histcounts(Edge_code(round(r/4+1):round(r/2),round(c/4+1):round(c/2)),256);
% [hist8,edges] = histcounts(Edge_code(round(r/4+1):round(r/2),round(c/2):round(3*c/4)),256);
% [hist9,edges] = histcounts(Edge_code(round(r/2+1):round(3*r/4),round(c/4+1):round(c/2)),256);
% [hist10,edges] = histcounts(Edge_code(round(r/2+1):round(3*r/4),round(c/2+1):round(3*c/4)),256);
% 
% % Concatenate 10 histograms
% 
% Hist_Edge_code=[hist1 hist2 hist3 hist4 hist5 hist6 hist7 hist8 hist9 hist10];
% 
% %% Histogram for Direction Code
% 
% [hist1,edges] = histcounts(Direction_code,16); 
% [hist2 edges] = histcounts(Direction_code(1:round(r/2),1:round(c/2)),16);
% [hist3,edges] = histcounts(Direction_code(1:round(r/2),round(c/2+1):c),16);
% [hist4,edges] = histcounts(Direction_code(round(r/2+1):r,1:round(c/2)),16);
% [hist5,edges] = histcounts(Direction_code(round(r/2+1):r,round(c/2+1):c),16);
% [hist6,edges] = histcounts(Direction_code(round(r/4+1):round(3*r/4),round(c/4+1):round(3*c/4)),16);
% [hist7,edges] = histcounts(Direction_code(round(r/4+1):round(r/2),round(c/4+1):round(c/2)),16);
% [hist8,edges] = histcounts(Direction_code(round(r/4+1):round(r/2),round(c/2):round(3*c/4)),16);
% [hist9,edges] = histcounts(Direction_code(round(r/2+1):round(3*r/4),round(c/4+1):round(c/2)),16);
% [hist10,edges] = histcounts(Direction_code(round(r/2+1):round(3*r/4),round(c/2+1):round(3*c/4)),16);
% 
% % Concatenate 10 histograms
% 
% Hist_Direction_code=[hist1 hist2 hist3 hist4 hist5 hist6 hist7 hist8 hist9 hist10];
% 
% %% Histogram for Pixel Code
% 
% [hist1,edges] = histcounts(Pixel_code,2); 
% [hist2 edges] = histcounts(Pixel_code(1:round(r/2),1:round(c/2)),2);
% [hist3,edges] = histcounts(Pixel_code(1:round(r/2),round(c/2+1):c),2);
% [hist4,edges] = histcounts(Pixel_code(round(r/2+1):r,1:round(c/2)),2);
% [hist5,edges] = histcounts(Pixel_code(round(r/2+1):r,round(c/2+1):c),2);
% [hist6,edges] = histcounts(Pixel_code(round(r/4+1):round(3*r/4),round(c/4+1):round(3*c/4)),2);
% [hist7,edges] = histcounts(Pixel_code(round(r/4+1):round(r/2),round(c/4+1):round(c/2)),2);
% [hist8,edges] = histcounts(Pixel_code(round(r/4+1):round(r/2),round(c/2):round(3*c/4)),2);
% [hist9,edges] = histcounts(Pixel_code(round(r/2+1):round(3*r/4),round(c/4+1):round(c/2)),2);
% [hist10,edges] = histcounts(Pixel_code(round(r/2+1):round(3*r/4),round(c/2+1):round(3*c/4)),2);
% 
% % Concatenate 10 histograms 
% Hist_Pixel_code=[hist1 hist2 hist3 hist4 hist5 hist6 hist7 hist8 hist9 hist10];


% FeatureCode=[Hist_Edge_code Hist_Direction_code Hist_Pixel_code];

result=FeatureCode;


end










