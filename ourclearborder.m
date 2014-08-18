% This function eliminates cells on the edge of the image
function Lout = ourclearborder(Lin)

Lout = zeros(size(Lin));
ind = 1;
for i = 1:max(max(Lin)) 
    L1 = Lin;
    L1(find(L1~=i))=0;
    L1 = imclearborder(L1);
    Lout = Lout + L1; 
end