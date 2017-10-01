function PlotIt2MFs(varargin)
if nargin<=3
    n=varargin{1};
    J_max=varargin{2};
    Net=varargin{3};
    %xIntval=Net(1).temp;
elseif nargin<=4
    n=varargin{1};
    J_max=varargin{2};
    Net=varargin{3};
    xIntval=varargin{4};
end 
        
x=linspace(xIntval(1),xIntval(2),500)';
line_color=lines(n);
for i=1:n
    subplot(n,1,i)
    for j=1:J_max
        line_color_j=lines(J_max);
        UMFVal=exp(-( (x-Net(2).m{1}(i,j)) ./ Net(2).sigma{1}(i,j) ).^2 ) ;
        LMFVal=exp(-( (x-Net(2).m{2}(i,j)) ./ Net(2).sigma{2}(i,j) ).^2 ) ;
        plot(x,UMFVal,'-', 'Color',line_color_j(j,:))
        hold on
        plot(x,LMFVal,'-.','LineWidth',0.8, 'Color',line_color_j(j,:))
    end
    
end
xlabel('x(t-3),x(t-2),x(t-1),x(t)')
end