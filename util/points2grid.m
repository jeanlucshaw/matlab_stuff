function [varargout] = points2grid(x,y,m_map,varargin)
% Synthax :       points2grid(x,y,m_map,varargin)
%
% Takes as input a grid of coordinates contained in 'x' and 'y' with the
% horizontal coordinate being the columns of rows of 'x' and the vertical
% coordinate being the columns of 'y'. Values with a real coordinate to be
% the center of a grid cell and nan coordinate values to correspond to a
% dead grid cell. Draws a regularly spaced grid around the live grid
% coordinates on the current active figure. Puts hold to on.
%
% The 'm_map' input variable indicated if the function is called for a
% figure drawn in an m_map projection. It should be the string 'y' if 
% this is the case, and 'n' (or anything else) if it is not.
%
% The subsequent arguments are taken as the property/value pair of the
% plotted grid itself.
%

% Options 
% default
% 

% input size
[M,N] = size(x) ;

% index grid
[indexx,indexy] = meshgrid(1:2*M+1,1:2*N+1) ;

% Interpolation
[I,J]   = find(isfinite(x) & isfinite(y)) ;
L       = find(isfinite(x) & isfinite(y)) ;
xvec    = x(L) ; 
yvec    = y(L) ; 
INTERPx = scatteredInterpolant(2*I,2*J,xvec(:), ...
            'linear','linear') ;
INTERPy = scatteredInterpolant(2*I,2*J,yvec(:), ...
            'linear','linear') ;

% Extrapolation
X = INTERPx(indexx',indexy') ;
Y = INTERPy(indexx',indexy') ;

% Save full grid for output
FX = X ;
FY = Y ;

% Replace initial nan values
[I,J]   = find(isnan(x) | isnan(y)) ;
for ii = 1:numel(I)
    X(2*I(ii),2*J(ii)) = nan ;
    Y(2*I(ii),2*J(ii)) = nan ;
end

[Inanx,Inany] = find(isnan(X) | isnan(Y)) ;

% nans in the grid lines
for ii = 1:numel(Inanx)
    % ========================================
    if Inany(ii) == 2 % left border node
        X(Inanx(ii),Inany(ii)-1) = nan ;
        if Inanx(ii) == 2    % corner bot node
            X(Inany(ii)-1,Inanx(ii)  ) = nan ;
            X(Inany(ii)-1,Inanx(ii)-1) = nan ;
        
        elseif Inanx(ii) == 2*M % corner top node
            X(Inanx(ii)+1,Inany(ii)  ) = nan ;
            X(Inanx(ii)+1,Inany(ii)-1) = nan ;
        else
            if isnan(X(Inanx(ii)-2,Inany(ii)))% adjacent up is also nan
                X(Inanx(ii)-1,Inany(ii)  ) = nan ;
                X(Inanx(ii)-1,Inany(ii)-1) = nan ;
            end
            if isnan(X(Inanx(ii)+2,Inany(ii)))% adjacent dw is also nan
                X(Inanx(ii)+1,Inany(ii)  ) = nan ;
                X(Inanx(ii)+1,Inany(ii)-1) = nan ;
            end
            if isnan(X(Inanx(ii),Inany(ii)+2))% adjacent rg is also nan
                X(Inanx(ii),Inany(ii)+1  ) = nan ;
            end
        end
    

    % ========================================
    elseif Inany(ii) == 2*N % right border node
        X(Inanx(ii),Inany(ii)+1) = nan ;
        if Inanx(ii) == 2    % corner bot node
            X(Inanx(ii)-1,Inany(ii)  ) = nan ;
            X(Inanx(ii)-1,Inany(ii)+1) = nan ;
        
        elseif Inanx(ii) == 2*M % corner top node
            X(Inanx(ii)+1,Inany(ii)  ) = nan ;
            X(Inanx(ii)+1,Inany(ii)+1) = nan ;
        else
            if isnan(X(Inanx(ii)-2,Inany(ii)))% adjacent up is also nan
                X(Inanx(ii)-1,Inany(ii)  ) = nan ;
                X(Inanx(ii)-1,Inany(ii)+1) = nan ;
            end
            if isnan(X(Inanx(ii)+2,Inany(ii)))% adjacent dw is also nan
                X(Inanx(ii)+1,Inany(ii)  ) = nan ;
                X(Inanx(ii)+1,Inany(ii)+1) = nan ;
            end    
            if isnan(X(Inanx(ii),Inany(ii)-2))% adjacent lf is also nan
                X(Inanx(ii),Inany(ii)-1  ) = nan ;
            end
        end
    
    
    
    % ========================================
    elseif Inanx(ii) == 2 & Inany(ii)~=2 & Inany(ii)~=2*N % bottom node
            X(Inanx(ii)-1,Inany(ii)  ) = nan ;
        if isnan(X(Inanx(ii),Inany(ii)+2))
            X(Inanx(ii)  ,Inany(ii)+1) = nan ;
            X(Inanx(ii)-1,Inany(ii)+1) = nan ;
        end
        if isnan(X(Inanx(ii),Inany(ii)-2))
            X(Inanx(ii)  ,Inany(ii)-1) = nan ;
            X(Inanx(ii)-1,Inany(ii)-1) = nan ;
        end   
        if isnan(X(Inanx(ii)+2,Inany(ii)))% adjacent up is also nan
            X(Inanx(ii)+1,Inany(ii)  ) = nan ;
        end
        if isnan(X(Inanx(ii),Inany(ii)+2))% adjacent rg is also nan
            X(Inanx(ii),Inany(ii)+1  ) = nan ;
            X(Inanx(ii)-1,Inany(ii)+1  ) = nan ;
        end
        if isnan(X(Inanx(ii),Inany(ii)-2))% adjacent lf is also nan
            X(Inanx(ii),Inany(ii)-1  ) = nan ;
            X(Inanx(ii)-1,Inany(ii)-1  ) = nan ;
        end
    
    
    % ========================================
    elseif Inanx(ii) == 2*M & Inany(ii)~=2 & Inany(ii)~=2*N % top node
        X(Inanx(ii)+1,Inany(ii)  ) = nan ;
        if isnan(X(Inanx(ii),Inany(ii)+2))
            X(Inanx(ii)  ,Inany(ii)+1) = nan ;
            X(Inanx(ii)+1,Inany(ii)+1) = nan ;
        end
        if isnan(X(Inanx(ii),Inany(ii)-2))
            X(Inanx(ii)  ,Inany(ii)-1) = nan ;
            X(Inanx(ii)+1,Inany(ii)-1) = nan ;
        end        
        if isnan(X(Inanx(ii)-2,Inany(ii)))% adjacent dw is also nan
            X(Inanx(ii)-1,Inany(ii)  ) = nan ;
        end
        if isnan(X(Inanx(ii),Inany(ii)+2))% adjacent rg is also nan
            X(Inanx(ii),Inany(ii)+1  ) = nan ;
            X(Inanx(ii)+1,Inany(ii)+1  ) = nan ;
        end
        if isnan(X(Inanx(ii),Inany(ii)-2))% adjacent lf is also nan
            X(Inanx(ii),Inany(ii)-1  ) = nan ;
            X(Inanx(ii)+1,Inany(ii)-1  ) = nan ;
        end
    else % INSIDE THE GRID
            if isnan(X(Inanx(ii),Inany(ii)+2)) % right neighboor
                X(Inanx(ii)  ,Inany(ii)+1) = nan ;
            end
            if isnan(X(Inanx(ii),Inany(ii)-2)) % left neighboor
                X(Inanx(ii)  ,Inany(ii)-1) = nan ;
            end
            if isnan(X(Inanx(ii)-2,Inany(ii))) % up   neighboor
                X(Inanx(ii)-1,Inany(ii)) = nan ;
            end
            if isnan(X(Inanx(ii)+2,Inany(ii))) % down neighboor
                X(Inanx(ii)+1,Inany(ii)) = nan ;
            end
    end
    
end

% same nan values in Y
I    = find(isnan(X)) ;
Y(I) = nan ;

if startsWith(m_map,'y')
for ii = 1:2:2*M+1
    m_plot(X(ii,:),Y(ii,:),varargin{:})
end

for jj = 1:2:2*N+1
    m_plot(X(:,jj),Y(:,jj),varargin{:})
end

else

for ii = 1:2:2*M+1
    plot(X(ii,:),Y(ii,:),varargin{:})
end

for jj = 1:2:2*N+1
    plot(X(:,jj),Y(:,jj),varargin{:})
end
end

% Manage output
if 	nargout==2
	varargout{1} = X ; 
	varargout{2} = Y ; 
elseif	nargout==4
	varargout{1} = X ; 
	varargout{2} = Y ; 
	varargout{3} = FX ; 
	varargout{4} = FY ; 
elseif	nargout==0
else
	varargout = [] ;
end
end
