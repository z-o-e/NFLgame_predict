% The goals of this problem are to use historical football scores to estimate values 
% of ¦Ìi and ¦Ái for all 32 teams in the National Football League, then use the resulting
 % models to attempt to predict the winners of all 16 football games that will be
 % played between 11/22/2012-11/26/2012

% The Matlab file games.mat contains the following data: 

 % The matrix teams contains the names of the 32 NFL teams. The index of a team 
 % is given by the team¡¯s position in this matrix.

 % The matrix scores contains all scores for the 160 games that have been played
 % in the 2012 NFL season. Each row contains the matchup and final score of a 
 % game. Column 1 contains the index of the visiting team. Column 2 contains the
 % index of the home team. Column 3 contains the number of points scored by the 
 % visiting team. Column 4 contains the number of points scored by the home team.
 
 % The matrix games contains all matchups for the 16 games to be played between
 % 11/22/2012-11/26/2012. Each row contains a pair of teams that is playing a 
 % game. Column 1 contains the index of the visiting team. Column 2 contains the
 % index of the home team.
%
% author huiyingzhang

load ('games.mat')
% set initial point

n1 = size(32,1);
n2 = size(32,1);
x1 = ones(32,1);
x2 = ones(32,1);
I1 = zeros(32);
I2 = zeros(32);
I3 = zeros(32);
I4 = zeros(32) ;


%home team indicator
for i = 1:160
    for j = 1:32
        if (scores(i,2) == j)
            I1(i,j) = 1;
            
        end
    end
end

  %away team indicator
  for i = 1:160
    for j = 1:32
        if (scores(i,1) == j)
            I2(i,j) = 1;
            
        end
    end
  end

  %home team score
  for i = 1:160
    for j = 1:32
        if (scores(i,2) == j)
            I3(i,j) = -scores(i,3);
            
        end
    end
  end

  %away team score
  for i = 1:160
    for j = 1:32
        if (scores(i,1) == j)
            I4(i,j) = -scores(i,4);
            
        end
    end
  end

  
c1 = [I1 I4];
c2 = [I2 I3];
x = [x1;x2];
m1 = zeros(32,1);
m2 = ones(32,1);
m = [m1;m2];
N=zeros(1,64);
for i=1:32
    N(i+32)=sum(I1(:,i))+sum(I2(:,i));
end


Lnalph = log(x);
sigsq=0.1;
cons=1/(2*sigsq);
y=cons*x'*(c1'*c1)*x + cons*x'*(c2'*c2)*x+N*Lnalph;
g=2*cons*(c1'*c1)*x+2*cons*(c2'*c2)*x-diag(1./x)*N';
H=2*cons*(c1'*c1)+2*cons*(c2'*c2)+diag(N)*diag((1./x).^2);

 

while(g'*g>0.01)   %check norm square of gradient 
    
    u=-H\g;   %choose search direction
    
    tl=0;
    tu=1;   %pick initial interval for bisection
    
     gt=2*cons*(c1'*c1)*x+2*cons*(c2'*c2)*x-diag(1./x)*N';   %evaluate the gradient of the current direction at tu
    
    if(max((x+tu*u)'*m)>0 || gt'*u>0)%if derivative is positive or new point is infeasible(note that in S A*x<b) at tu, then start bisection
        while((tu-tl)>0.01)   %check if the interval of bisection is within tolerance
            t=(tl+tu)/2;    %find the midpoint
            
            gt=2*cons*(c1'*c1)*x+2*cons*(c2'*c2)*x-diag(1./x)*N';  %calculte the gradient at midpoint
            
            if(max((x+t*u)'*m)>0 || gt'*u>0)    %check bisection criteria
                tu=t;
            else
                tl=t;
            end
        end
    else
        t=tu;   %if derivative is negtive 
    end
            
    x=x+t*u;    %update x
   
   
g=2*cons*(c1'*c1)*x+2*cons*(c2'*c2)*x-diag(1./x)*N';
H=2*cons*(c1'*c1)+2*cons*(c2'*c2)+diag(N)*diag((1./x).^2); %compute gradient for next iteration     
end

xstar=x;

mu=xstar(1:32)
alpha=xstar(33:64)
disp('Winners:');
for k=1:16
    ik=games(k,2);
    jk=games(k,1);
    if(mu(ik)/alpha(jk)>mu(jk)/alpha(ik))
        disp(teams(ik,:));
    else
        disp(teams(jk,:));
    end;
end;

