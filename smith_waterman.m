%144133E - M.G.K.C.Piyarathna
%Assignment 2 - BioInfomatics
%Implementing matlab code for Smith Waterman Algorithm
function [] =smith_waterman(sequence_1,sequence_2)
%function is a group of statements that together perform a task.The name of the file and of the function should be the same.
%I have input two parameters.Those are Sequences we have to provide
%Take in inputs (sequence_1, sequence_2)
%assign sequence_1 to variabe A
A = sequence_1;
%assign sequence_2 to variabe A
B = sequence_2;

%Here I have used gap penalty = -2,match = 2, mismatch = -3
d = -2;%gap penalty
m = 2;%match
s = -1;%mismatch

col = length(A);%return the length of A and assign it to the column
row = length(B);%return the length of B and assign it to the row
% rows|columns

%preallocating matrix F all zeros
F = zeros(row+1,col+1); %create an array of Zeros/preallocating matrix F all zeros

%Filling in the matrix
for i=2:row+1
    for j=2:col+1        
        %j is going through the columns which is A
        %if the two positions match,scoring is done accordingly
        
        if (A(j-1) == B(i-1))
            scores(A(j-1),B(i-1)) = m;
        else
            scores(A(j-1),B(i-1)) = s;
        end        
        %Filling-in partial alignments here  
        Diagonal     = F(i-1,j-1) + scores(A(j-1),B(i-1));
        Left = F(i, j-1) + d;
        Right = F(i-1, j) + d;
        %computing the final score of the alignment and assigning it to F
        Temp = [Diagonal Left Right 0];
        F(i,j) = max(Temp);% return the maximum of Temp
          % end
    end
end
%traceback part
%find max element
[rr,cc]=(find(F==max(F(:))));

Alignment_A = '';
Alignment_B = '';
i = rr; %traceback starts with the maximum element of F
j = cc; %max column
TotalScore = 0;
while (F(i,j) ~= 0 )
   Score = F(i,j);
   DIAG = F(i-1,j-1);
   %LEFT = F(i-1,j);
   UP   = F(i,j-1);
   
   %if scores are equal to the diagonal, there is no gap.
   if (Score == DIAG + scores(A(j-1),B(i-1)))
       Alignment_A = strcat(Alignment_A, A(j-1));
       Alignment_B = strcat(Alignment_B, B(i-1));
       %computes score, checks to see if the alignment are the same
       %characters
       if (A(j-1) == B(i-1))
           TotalScore = TotalScore + 1;
       else
           %mismatch
           TotalScore = TotalScore + s;
       end
       i = i-1;
       j = j-1;
   
   %gap in sequence B
   elseif (Score == UP + d)
       Alignment_A = strcat(Alignment_A, A(j-1));%String concatenation
       Alignment_B = strcat(Alignment_B, '-');%String concatenation
       j = j-1;
       
   %gap in sequence A
   else
       Alignment_A = strcat(Alignment_A, '-');%String concatenation
       Alignment_B = strcat(Alignment_B, B(i-1));%String concatenation
       i = i-1;
   end
end


%displays the alignment score and alignments
disp('Score Matrix');%print output
i = 1;
while(i<=row+1)
    disp(F(i,1:end));%print output
    i=i+1;
end
matches = 0;
for i=1:length(Alignment_A) 
    if( Alignment_A(i) == Alignment_B(i))
        matches = matches + 1;
    end
end   
matches = 0;
for i=1:length(Alignment_A)
    if( Alignment_A(i) == Alignment_B(i))
        matches = matches + 1;
    end
end    
disp (strcat('Number of Matches :',num2str(matches)));%String concatenation --print output %num2str- converts the array A into a string

disp('Final alignment for given two sequences');%print output

disp(Alignment_A);%print output
disp(Alignment_B);%print output