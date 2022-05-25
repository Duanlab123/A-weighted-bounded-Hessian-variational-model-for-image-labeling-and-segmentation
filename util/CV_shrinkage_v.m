
function v = CV_shrinkage_v(temp_v)

[l1,l2] = size(temp_v);

v = temp_v;

for i = 1 : l1
    for j = 1 : l2        
        if( temp_v(i,j)>1 )
            v(i,j) = 1;
        else
            if(temp_v(i,j)<0)
                v(i,j) = 0;
            end
        end
    end
end



