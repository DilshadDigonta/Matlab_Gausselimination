%Gauss elemination
function x = gausselimintation(A,b)
A = input('Please Enter the Co efficient Mattrix:' );
b = input('Please Enter the constant:');

Aug = [A b];                                             %Augmented Matrix
[m, n]= size(Aug);                                      %Size of Aug
for i = 1: m                                            %Loop for Row Operation
    for j = i+1 : m                                    %Loop for the multiplier for Pivoting 
        mul = Aug(j,i) ./ Aug(i,i);                     %Multiplier for Pivoting
        Aug(j,:) = Aug(j,:) - mul.* Aug(i,:);           %Finding the upper triangular Matrix 
    end
end

x = zeros(1,size(A,2));
for i = size(Aug,1):-1:1
    d = sum(Aug(i,i+1:end-1).*x(i+1:end));
    x(i) = (Aug(i,end)-d)./ Aug(i,i);
end
fprintf('solve for x in gauss elimination : %d \n',x)