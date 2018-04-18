

% List of Variables we shocked
list=zeros(numel(paths_list),1);
for ii=1:numel(paths_list)
    eval(['aux=paths.' paths_list{ii} ';']);
    if sum(sum(aux~=aux(:,1)*ones(1,length(aux))))>1;
       list(ii)=1;
    end
end

% Experiment along the transition
if sum(list)==0
    disp('No Shocks');
else
    disp('There are Shocks to the following variables');
    for ii=1:numel(paths_list)
        if list(ii)==1
            disp(['Shock to ' paths_list_tit{ii}]);
            eval(['aux=paths.' paths_list{ii} ';']);
            if plotit==1
                figure(cc);
                plot(1:T,aux(1,1:T));
                title(['Path for ' paths_list_tit{ii}]); ylim([0 25]);
                axis tight;
                cc=cc+1;
            end
        end
    end
end
