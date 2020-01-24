function precalculation_DY_mult(maxORD)

    %maxORD = 5;
    stopping = 3;
    %%% indexes on multiplications DY terms
    DY_ij=cell(2,maxORD+stopping);
    DY_ijk=cell(3,maxORD+stopping);

    for k=2:maxORD+stopping-1
        mm1=k:-1:1;
        mm2=1:k;
        DY_ij{1,k-1}=mm1;
        DY_ij{2,k-1}=mm2;

        for jj=1:k
            mm1=k+1-jj:-1:1;
            mm2=(1:k-jj+1);
            DY_ijk{1,k-1}=[DY_ijk{1,k-1},mm1];
            DY_ijk{2,k-1}=[DY_ijk{2,k-1},mm2];
            DY_ijk{3,k-1}=[DY_ijk{3,k-1},jj*ones(1,length(mm1))];
        end
   end % k (ORD)

    filename=['DY_indexes_maxORD_',int2str(maxORD)];
    save(filename,'DY_ij','DY_ijk');
    
    DY_ijk = DY_reverse_order(DY_ijk);
    
    filename=['DY_indexes_maxORD_GN_',int2str(maxORD)];
    save(filename,'DY_ij','DY_ijk');
   
end

function var_flipped = DY_reverse_order(var)
    %% ordering coefficients (pro vytykani)
    % change order
    DY_edited=var(end:-1:1,:);
    rows = size(DY_edited,1);
    cols = size(DY_edited,2);
    var_flipped = cell(rows,cols);
    
    % flip individual cells 
    for i=1:rows
        for j=1:cols
            var_flipped{i,j} = flip(DY_edited{i,j});     
        end
    end
end

% reodering of columns - descending order, based on the last row, for DY4,
% DY5
function var_reordered = reorder_columns(var, numTerms, maxORD)

    var_reordered = var;
    for k=1:maxORD
        if k == 1
           t_new_start = 2; 
        else
            t_new_start = size(var{1,k-1},2)+1;
        end
        DY_mtx=cell2mat(var(:,k));
        DY_mtx_toSort = DY_mtx(:,t_new_start:end);
        DY_mtx_t = DY_mtx_toSort';


        DY_ordered=sortrows(DY_mtx_t,-numTerms);  % minus - 3,2,1, kladne - 1,2,3
        DY_ordered=DY_ordered';
        
        for j=1:numTerms
            var_reordered{j,k}(t_new_start:end)=DY_ordered(j,:);
        end
    %var_reordered{:,k}

    end
    
end
