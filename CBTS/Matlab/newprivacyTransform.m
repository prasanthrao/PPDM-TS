
function [result,corr_diff] = newprivacyTransform(inputData,privateColumns,type,comptype)

%privacyTransform
% Normalized transformation for continous numeric data accounting for privacy.
% Normalization range [0.1 0.9]
% [correlationMatrix, normalizedData, COEFF, latent, tsquare] = privacyTransform(inputData,privateColumns)
% correlationMatrix = correlation between columns
% normalizedData = normalized inputData for computation
% COEFF = PCA Co-effecients
    [N M] = size(inputData);
    inputData(isnan(inputData)) = 0;
    min_d = min(inputData(:));
    max_d = max(inputData(:));
    ra = 1;
    rb = 0;
    normalizedData = (((ra-rb) * (inputData - min_d)) / (max_d - min_d)) + rb;
    result = inputData;
    %Correlation Analysis of Data
    %correlationMatrix = corrcoef(normalizedData);
        
    correlationMatrix = corrcoef(normalizedData);
    
    for i=privateColumns
        this_corr = correlationMatrix(i,:);
        sort(this_corr);
        threshold = this_corr(ceil(M/4));
        indices = correlationMatrix(i,:)>=threshold;
        %inputData(1:15,indices)
        %if sum(indices) == 1
         %   [m n] = size(indices);
         %  n= round(n/10)+ 1 % %Percentage of Correlated Columns in case of threshold miss
         %   sortedArray=sort(abs(correlationMatrix(i,:)),'descend');
         %    sortedArray(isnan(sortedArray))=0;
         %   sortedArray = sort(sortedArray,'descend');
         %   indices = abs(correlationMatrix(i,:))>=sortedArray(1,n)
        %end  
        
        %transformedColumn = normalizedData(:,i);
        %val = abs(correlationMatrix(i,:))==1;
        switch comptype
            case 'pca'
                [COEFF,SCORE,latent,tsquare] = princomp(normalizedData(:,indices));
                normalizedData(:,i)=(SCORE(:,1));
                %converting mixed values to positive
                normalizedData(:,i) = normalizedData(:,i) + abs(min(normalizedData(:,i)));
                %transformedColumn = SCORE(:,1);
            case 'svd'
                x = 1;
                while N-x >= 99
                    st = x;
                    ed = x+99;
                    [u,s,v] = svd(normalizedData(st:ed,indices));
                    res = normalizedData(st:ed,indices)*v';
                    [n m] = size(res);
                    sel_col = normalizedData(st:ed,i);
                    max_cor = 0;
                    for j = 1:m
                        corr_mat = corrcoef(normalizedData(st:ed,i),res(:,j));
                        if(abs(corr_mat(2:2))>max_cor)
                            max_cor = abs(corr_mat(2:2));
                            sel_col = res(:,j);
                        end
                    end
                    normalizedData(st:ed,i) = abs(sel_col);
                    x = x+100;
                end
                
                st = x;
                [u,s,v] = svd(normalizedData(st:end,indices));
                res = normalizedData(st:end,indices)*v';
                [n m] = size(res);
                sel_col = normalizedData(st:end,i);
                max_cor = 0;
                for j = 1:m
                    corr_mat = corrcoef(normalizedData(st:end,i),res(:,j));
                    if(abs(corr_mat(2:2))>max_cor)
                        max_cor = abs(corr_mat(2:2));
                        sel_col = res(:,j);
                    end
                end
                normalizedData(st:end,i) = abs(sel_col);
            case 'nnmf' 
                %normalizedData=abs(normalizedData);
                [r,c] = size(normalizedData(:,indices));
                max_corr = 0;
                sel_col = normalizedData(:,i);
                for j = 1:(min(r,c))
                    [a,b] = nnmf(normalizedData(:,indices),j);
                    res=a*b;
                    corr_mat = corrcoef(res(:,1),normalizedData(:,i));
                    if(abs(corr_mat(2:2))>max_corr && abs(corr_mat(2:2))~=1)
                        max_corr = abs(corr_mat(2:2));
                        sel_col = res(:,1);
                    end 
                end
                normalizedData(:,i) = sel_col;
        end
        ra = max(inputData(:,i));
        rb = min(inputData(:,i));
        max_d = max(normalizedData(:,i));
        min_d = min(normalizedData(:,i));
        if strcmp(type,'int')
            result(:,i) = round((((ra - rb) * (normalizedData(:,i) - min_d)) / (max_d - min_d)) + rb);
        elseif strcmp(type,'real')
            result(:,i) = ((((ra - rb) * (normalizedData(:,i) - min_d)) / (max_d - min_d)) + rb);
        end
    end
    
    % denormalize

    corr_ds = corrcoef(inputData);
    corr_res = corrcoef(result);
    corr_ds(isnan(corr_ds)) = 0;
    corr_res(isnan(corr_res)) = 0;
    % average change in correlation matrix
    corr_diff = sum(sum(abs(corr_ds-corr_res)))/(M*M);
    e = entropy(result);
    
    % Results as on 14-sept-2014
    % 
    % Data : 'C:\Users\hkbharath\Documents\Data sets\bcancer.csv'
    % Private Col : [2 5 7 9]
    %
    %
    % Threshold : 0.75
    % Method : 'nnmf' Result(avg correlation difference) : 0.0442-0.0454 (changes)
    % Method : 'svd(with negetive to positive conversion with abs)' 
    %           Result(avg correlation difference) : 0.0450
    % Method : 'pca(with negetive to positive conversion with addition)' 
    %           Result(avg correlation difference) : 0.0576
    % Problems : 'pca' method produces numbers like [8.67361737988404e-15]
    %             in [9] which is undesirable and unstable.
    %             values of [9] does not change in nnmf and svd
    %
    % Threshold : 0.9
    % Method : 'nnmf' Result(avg correlation difference) : 0.0034 (not changes)
    % Method : 'svd(with negetive to positive conversion with abs)' 
    %           Result(avg correlation difference) : 0.0034
    % Method : 'pca(with negetive to positive conversion with addition)'
    %           Result(avg correlation difference) : 0.0034
    % Problems: 'pca' method produces numbers like [8.67361737988404e-15]
    %           in [9] and [1] which is undesirable and unstable.
    %           'nnmf' and 'svd' methods didnt changed values of 
    %           private coloums [5 7 9]
    % Conclusion : less threshold provides more possiblity of change in all
    %           private coloums
    %
    %
    % Threshold : 0.6
    % Method : 'nnmf' Result(avg correlation difference) : 0.0707-0.0727 (changes)
    % Method : 'svd(with negetive to positive conversion with abs)' 
    %           Result(avg correlation difference) : 0.0642
    % Method : 'pca(with negetive to positive conversion with addition)'
    %           Result(avg correlation difference) : 0.0763
    % Problems: 'pca' method produces numbers like [8.67361737988404e-15]
    %           in [9] which is undesirable and unstable.
    %           values of [9] does not change in nnmf and svd
    %
    %
    % Threshold : 0.7
    % Method : 'nnmf' Result(avg correlation difference) : 0.0507-0.0523 (changes)
    % Method : 'svd(with negetive to positive conversion with abs)' 
    %           Result(avg correlation difference) : 0.0369
    % Method : 'pca(with negetive to positive conversion with addition)'
    %           Result(avg correlation difference) : 0.0675
    % Problems: 'pca' method produces numbers like [8.67361737988404e-15]
    %           in [9] which is undesirable and unstable.
    %           values of [9] does not change in nnmf and svd