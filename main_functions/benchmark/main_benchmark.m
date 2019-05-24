
for i = 1:4
    for n = 1:25
        X = data{i,n};
        
        freq = computeSpectra(data,[],3,0);

        
    end
end

save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\main_functions\benchmark\data\data','data')
save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\main_functions\benchmark\data\CMat','CMat')