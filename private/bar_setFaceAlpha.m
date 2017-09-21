function handles = bar_setFaceAlpha(ax,alpha)
%BAR_SETFACEALPHA      set the FaceAlpha property of BAR and BARH objects
%
%    This is a simple workaround for the fact that Matlab's BAR and BARH
%    plotting functions do not pass the 'FaceAlpha' property to the actual
%    plotting functions.
%    Usefully, e.g.
%        * in case you want to have the grid lines shine through the bars
%        * in general in case you want to have some graphics behind your
%          bars
%        * or, in case you simply prefer the faced out colors.
%
%    BAR_SETFACEALPHA ... set FaceAlpha of any barplot child of the
%    current axes to 0.5
%    BAR_SETFACEALPHA(AX) ... set FaceAlpha of any barplot child of
%    axeshandle AX to 0.5 where AX may also be an array of axes handles.
%    BAR_SETFACEALPHA(BAR_HANDLES) ... set FaceAlpha of any barplot to 0.5
%    BAR_SETFACEALPHA(AXES,ALPHA) ... sets FaceAlpha to ALPHA
%    BAR_SETFACEALPHA(BAR_HANDLES,ALPHA) ... sets FaceAlpha to ALPHA
%    HANDLES = BAR_SETFACEALPHA(...)  ... returns the handles to the
%    objects changed
%
%    Hopefully, this function gets obsolete in some future Matlab release.
%
%    Examples:
%
%        bar(magic(6));
%        bar_setFaceAlpha
%
%        for i=1:4
%          ax(i) = subplot(2,2,i) ;
%          bar(magic(i+2));
%          set(gca,'YGrid','on');
%        end
%        set(gcf,'Position',round(get(0,'ScreenSize')/2))
%        bar_setFaceAlpha(ax)
%
%        h1 = bar(magic(4),'stacked');
%        bar_setFaceAlpha(h1,0.2);
%        hold on;
%        h2 = bar(magic(4));
%        bar_setFaceAlpha(h2,0.8);
%        hold off;
%
%
%    NOTE:  For BAR3 and BAR3H this is not necessary (and doesn't have an effect
%    either).
%    Simple do something like this:
%        h1 = bar3(magic(4),'detached');
%        hold on
%        h2 = bar3(magic(4).','detached');
%        h3 = bar3h(magic(4));
%        hold off
%        set(h1,'FaceColor','b','FaceAlpha',0.4)
%        set(h2,'FaceColor','r','FaceAlpha',0.6)
%        set(h3,'FaceColor','k','FaceAlpha',0.3)
%
%    See also BAR BARH FILL PATCH.

% Copyright (c) georg ogris ::: spantec.at ::: fall 2011
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the author nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

if ~exist('ax','var'), ax = gca; end
if ~exist('alpha','var'), alpha = 0.5; end
if alpha > 1
    alpha = 1 ;
    warning('Alpha should not be greater 1 nor smaller 0.') ;
end
if alpha < 0
    alpha = 0 ;
    warning('Alpha should not be greater 1 nor smaller 0.') ;
end

bh = ax ;
bh(arrayfun(@(x)~isprop(x,'BarLayout'),bh)) = [] ;
if isempty(bh)
    bh = get(ax,'Children') ;
    if iscell(bh)
        bh = cell2mat(bh) ;
    end
    bh(arrayfun(@(x)~isprop(x,'BarLayout'),bh)) = [] ;
end
bh = get(bh,'Children') ;
if iscell(bh)
    bh = cell2mat(bh) ;% FIXME test this!
end

set(bh,'FaceAlpha',alpha) ;

if nargout
    handles = bh ;
end
