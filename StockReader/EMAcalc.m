function EMAvec1 = EMAcalc(data,values)

s = 1;
i = 1;
    for j = values
       for j1 = 1:length(data) 
           if j1 == j
                EMAvec(s) = mean(data(1:j));
                s = s + 1;
           elseif j1 > j
                   EMAvec(s) = data(j1)*2/(j+1)+ EMAvec(s-1)*(1-2/(j+1));
                   s = s + 1;
           end 
       end
       

       EMAvec1{i} = EMAvec;  
       i = i+1;
       
       clear EMAvec
       s = 1;
       
    end