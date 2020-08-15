function [result] = ARBP(img, region_x, region_y)
%%%%%%% Size of Feature Vctor=512
patch_size_p_x=3;
patch_size_p_y=3;
k_max=2;
alpha=1;
[~, ~ ,channel]=size(img);
if(channel>1)
    img=rgb2gray(img);
end
img=double(img);
[r,c]=size(img);
result = [];


AMBP_code=zeros(r,c);

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
        his = zeros( 1, 2^(patch_size_p_x*patch_size_p_y));
        
        %Loops for a single local window
        for i = sr:er
            for j = sc:ec

                for k=1:k_max

                   if(i-k<1 && j-k<1)
                       S=img(1:i+k,1:j+k);
                   elseif(i+k>r && j+k>c)
                       S=img(i-k:r,j-k:c);
                   elseif(i-k<1 && j+k>c)
                       S=img(1:i+k,j-k:c);
                   elseif(i+k>r && j-k<1)
                       S=img(i-k:r,1:j+k);
                   elseif(i-k<1)
                       S=img(1:i+k,j-k:j+k); 
                   elseif( j-k<1)
                       S=img(i-k:i+k,1:j+k);                   
                   elseif(i+k>r)
                       S=img(i-k:r,j-k:j+k);
                   elseif(j+k>c)   
                       S=img(i-k:i+k,j-k:c);
                   else 
                       S=img(i-k:i+k,j-k:j+k); 
                   end
                   
                   Z_med=median(S);
                   Z_min=min(S);
                   Z_max=max(S);
                   Z_avg=mean(S);
                   [rr, cc]=size(S);
                   d=sum(img(i,j)-S(:,:))/(rr*cc-1);
                   if(abs(Z_med-img(i,j))<abs(Z_avg-img(i,j)))
                       Z_para=Z_avg;
                   else
                       Z_para=Z_med;
                   end
                   if((Z_min+alpha*d<Z_para) & (Z_max-alpha*d>Z_para))
                       break;
                   end
                  
                end
                if((Z_min+alpha*d<img(i,j)) & (Z_max-alpha*d>img(i,j)))
                    Threshold=img(i,j);
                else 
                    Threshold=Z_para;
                end
                
                p_x=floor(patch_size_p_x/2);
                p_y=floor(patch_size_p_y/2);
                
                if(i-p_x<1 && j-p_y<1)
                       img_patch=img(1:i+p_x,1:j+p_y);
                elseif(i+p_x>r && j+p_y>c)
                       img_patch=img(i-p_x:r,j-p_y:c);
                elseif(i-p_x<1 && j+p_y>c)
                       img_patch=img(1:i+p_x,j-p_y:c);
                elseif(i+p_x>r && j-p_y<1)
                       img_patch=img(i-p_x:r,1:j+p_y);
                elseif(i-p_x<1)
                       img_patch=img(1:i+p_x,j-p_y:j+p_y); 
                elseif( j-p_y<1)
                       img_patch=img(i-p_x:i+p_x,1:j+p_y);                   
                elseif(i+p_x>r)
                       img_patch=img(i-p_x:r,j-p_y:j+p_y);
                elseif(j+p_y>c)   
                       img_patch=img(i-p_x:i+p_x,j-p_y:c);
                else 
                       img_patch=img(i-p_x:i+p_x,j-p_y:j+p_y); 
                end
                               
                sum_value=0;
                for z = 1:size(img_patch)
                    if(img_patch(z)-Threshold>=0)
                        ss=1;
                    else
                        ss=0;
                    end
                    sum_value=sum_value+ss*2^(z-1);

                end
                

                AMBP_code( i,j ) = sum_value ;

                his( 1, AMBP_code( i,j ) + 1 ) = his( 1, AMBP_code( i,j ) + 1 ) + 1;

            end
        end

        result = [result his];
    end
end






