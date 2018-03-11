function [ p ] = diag_sum( M )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
    r  = size(M, 1);
    c = size(M, 2);
    num = r + c - 1;
    index_r = 1;
    index_c = 1;
    for k = 1: num
        sum = M(index_r, index_c);
        i_r = index_r;
        i_c = index_c;
        while (1)
            i_r = i_r - 1;
            i_c = i_c + 1;
            if i_r < 1
                break;
            end
            if i_c > c
                break;
            end
            sum = sum + M(i_r, i_c); 
        end
        p(k) = sum;
        
        if index_r < r
            index_r = index_r +1;
        else
            index_c = index_c + 1;
        end
        if index_c > c
            break;
        end     
    end
end

