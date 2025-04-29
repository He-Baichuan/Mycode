%当一二阶导数解析解位于工作区后，可用于计算数值解
%
%


nfx = zeros(size(fx));
nfx(:) = eval(fx(:));

nfxp = zeros(size(fxp));
nfxp(:)= eval(fxp(:));

nfy = zeros(size(fy));
nfy(:) = eval(fy(:));

nfyp = zeros(size(fyp));
nfyp(:)= eval(fyp(:));

nf = zeros(size(f));
nf(:)=eval(f(:));

if approx==1%若近似为1阶，二阶导数全为0
   

nfypyp=0; nfypy=0; nfypxp=0; nfypx=0; nfyyp=0; nfyy=0; nfyxp=0; nfyx=0; nfxpyp=0; nfxpy=0; nfxpxp=0; nfxpx=0; nfxyp=0; nfxy=0; nfxxp=0; nfxx=0;
   
   else

nfypyp=zeros(size(fypyp));
nfypyp(:)=eval(fypyp(:));

nfypy=zeros(size(fypy));
nfypy(:)=eval(fypy(:));

nfypxp=zeros(size(fypxp));
nfypxp(:)=eval(fypxp(:));

nfypx=zeros(size(fypx));
nfypx(:)=eval(fypx(:));

nfyyp=zeros(size(fyyp));
nfyyp(:)=eval(fyyp(:));

nfyy=zeros(size(fyy));
nfyy(:)=eval(fyy(:));

nfyxp=zeros(size(fyxp));
nfyxp(:)=eval(fyxp(:));

nfyx=zeros(size(fyx));
nfyx(:)=eval(fyx(:));

nfxpyp=zeros(size(fxpyp));
nfxpyp(:)=eval(fxpyp(:));

nfxpy=zeros(size(fxpy));
nfxpy(:)=eval(fxpy(:));

nfxpxp=zeros(size(fxpxp));
nfxpxp(:)=eval(fxpxp(:));

nfxpx=zeros(size(fxpx));
nfxpx(:)=eval(fxpx(:));

nfxyp=zeros(size(fxyp));
nfxyp(:)=eval(fxyp(:));

nfxy=zeros(size(fxy));
nfxy(:)=eval(fxy(:));

nfxxp=zeros(size(fxxp));
nfxxp(:)=eval(fxxp(:));

nfxx=zeros(size(fxx));
nfxx(:)=eval(fxx(:));


end 