

function HH_decision(wseq, rseq, Penseq, trseq, τ_w_seq, τ_b_seq, sp1, age_effciency, param)
    @unpack periods_working, periods_retire, nages, nA, n_ϵ, n_η,
            vGridAsset, τ_k, a_max, a_min, nages, γ, λ, σ, β,
            vGrids_ϵ, vGrids_η, mTrans_η = param 
    V = zeros(nages, n_ϵ, n_η, nA)
    A_policy = zeros(nages, n_ϵ, n_η, nA)
    C_policy = zeros(nages, n_ϵ, n_η, nA)
    L_policy = zeros(nages, n_ϵ, n_η, nA)
    psi1 = 0.0001
    #retirement's dicision
    for ieps in 1:n_ϵ, ieta in 1:n_η, ia in 1:nA
        c = Penseq[end,ieps] + trseq[end] + (1.0 + (1.0 - τ_k)*rseq[end])*vGridAsset[ia] 
        C_policy[end,ieps,ieta,ia] = c 
        #u = utility(c,0.0,param)
        if c > 0.0
            V[end,ieps,ieta,ia] = utility(c,0.0,param) 
        else
            V[end,ieps,ieta,ia] = -1e+12
        end
    end
    #向前迭代，由s = J, 到s = R 
    Gpow = (1.0 + γ)^(λ*(1.0 - σ)) 
    for s in (periods_retire-1):-1:1
        j = periods_working + s
        jp = j + 1
        for ieta in 1:n_η
            for ieps in 1:n_ϵ
                # 下一期退休者值函数 V_{j+1}(a')，长度 na
                Vnext = V[jp, ieps, ieta, :]
                
                    # （可选）检查 V_{j+1} 随资产单调非减
                    
                    @inbounds for k in 1:nA-1
                        if Vnext[k+1] < Vnext[k]
                            error("value function not monotone at (j=$jp, ieta=$ieta, ieps=$ieps)")
                        end
                    end
                    
                 # ======= 单调性热启动索引：必须放在 ia 循环外 =======
                m0 = 1
                # 固定 (j,ieps) 的回报与非资产收入
                R = 1.0 + (1.0 - τ_k) * rseq[j]
                y = Penseq[j, ieps] + trseq[j]
                cont_pref = β * sp1[j] * Gpow   # beta1*sp1(t+it)*(1+g)^(gam*(1-sigma))
                for ia in 1:nA
                    a = vGridAsset[ia]
                    cmax = y + R*a 
                    if cmax <= 0.0
                        A_policy[j, ieps, ieta, ia] = 0.0
                        C_policy[j, ieps, ieta, ia] = 0.0
                        V[j, ieps, ieta, ia]  = -1e+12
                        continue
                    end
                    # 可行上界：c>=0 ⇒ a' <= cmax/(1+g)，再截断到 assetmax
                    a_upper  = min(a_max, cmax / (1.0 + γ))
                    id_max  = searchsortedlast(vGridAsset, a_upper)     # 最大可行网格索引
                    if id_max < 1
                        A_policy[j, ieps, ieta, ia] = a_min
                        C_policy[j, ieps, ieta, ia] = cmax
                        V[j, ieps, ieta, ia]  = -1e+12
                        continue
                    end 
                        # rvalue(a'): 用你给的线性插值 lip（并在边界处钳制）
                        rvalue = ap -> begin
                            if ap <= a_min
                                return Vnext[1]
                            elseif ap >= a_max
                                return Vnext[end]
                            else
                                return lip(vGridAsset, Vnext, ap; on_extrap=:error)
                            end
                        end
                        # value1(a'): 直接复刻你 MATLAB 的 Bellman RHS
                        value1 = ap -> begin
                            c = cmax - (1.0 + γ)*ap
                            if c <= 0.0
                                return -1e+12
                            else
                                return utility(c, 0.0, param) + cont_pref * rvalue(ap)
                            end
                        end
                        # 用二分法搜索最优资产 a'，使得 value1(a') 最大# ======= 利用“最优索引单调右移”的爬坡括区间 =======
                        # 从 m0-2 附近开始向右爬，直到 value1 下降
                        m = clamp(m0 - 2, 1, max(1, id_max-1))
                        v = value1(vGridAsset[m])
                        while m < id_max
                            vnext = value1(vGridAsset[m+1])
                            if vnext >= v
                                m += 1
                                v  = vnext
                            else
                                break
                            end
                        end
                        #m0 = m  # 更新热启动索引（下一个 ia 会从这里附近开始）
                        # 峰附近小区间 [xl,xr]
                        l = max(1,  m-1)
                        r = min(id_max,  m+1)
                        xl = vGridAsset[l]
                        xr = vGridAsset[r]
                        # ======= 角点处理（与 MATLAB 的 psi1 测试同思路） =======
                        a1 = 0.0
                        if m == 1
                            # 左端角点 a'=0 ?
                            if value1(xl + psi1) >= value1(xl)
                                a1 = GoldenSearch(value1, xl, xr)
                            else
                                a1 = xl
                            end
                        elseif m == id_max
                            # 右端角点 a'=U(≈agrid[μ]) ?
                            xr_eps = max(a_min, xr - psi1)
                            if value1(xr) > value1(xr_eps)
                                a1 = xr
                            else
                                a1 = GoldenSearch(value1, xl, xr)
                            end
                        else
                            a1 = GoldenSearch(value1, xl, xr)
                        end
                        c = cmax - (1.0 + γ)*a1

                        A_policy[j, ieps, ieta, ia] = a1
                        C_policy[j, ieps, ieta, ia] = c
                        V[j, ieps, ieta, ia]  = value1(a1)
                        # （可选）检查本期 V_j(a) 随 a 单调
                        
                        if ia > 1 && V[j, ieps, ieta, ia] < V[j, ieps, ieta, ia-1]
                            error("no monotonicity in V at (j=$j, ieta=$ieta, ieps=$ieps)")
                        end
                        
                    
                end


            end    
        end
    end
    
            @views for s in periods_working:-1:1
                for ieps in 1:n_ϵ
                    
                    Vnext = V[s+1, ieps, :, :] 
                    #（可选）检查每个 η' 下 V_{t+1}(·) 对资产单调
                    
                    @inbounds for ieta in 1:n_η
                        for ia in 1:nA-1
                            if Vnext[ieta, ia+1] < Vnext[ieta, ia]
                                error("Vnext not monotone: s=$(s+1), ieps=$ieps, ieta=$ieta")
                            end
                        end
                    end
                    
                    for ieta in 1:n_η
                        m0 = 1
                        w0 = wseq[s] * age_effciency[s] * vGrids_ϵ[ieps] * vGrids_η[ieta]
                        coef = β * sp1[s] * Gpow
                        for ia in 1:nA
                            a0 = vGridAsset[ia]
                            # 目标函数 wvalue(ap)（闭包），内部用 π_eta 做 η' 的期望
                            wvalue = ap -> begin
                                    # labor + consumption（与你 MATLAB 同式）
                                
                                    c = λ * ( (1.0 + (1.0-τ_k)*rseq[s]) * a0 +
                                                trseq[s] +
                                                (1.0 - τ_w_seq[s] - τ_b_seq[s]) * w0 -
                                                (1.0 + γ) * ap )
                                    if c <= 0.0
                                        return -1e+12
                                    end
                                    l = 1.0 - c / ((1.0 - τ_w_seq[s] - τ_b_seq[s]) * w0) * (1.0-λ)/λ
                                    if l < 0.0
                                        l = 0.0
                                        c = (1.0 + (1.0-τ_k)*rseq[s]) * a0 + trseq[s] - (1.0+γ)*ap
                                        if c <= 0.0
                                            return -1e+12
                                        end
                                    end
                                    # EV = Σ_{η'} π(η→η') V_{t+1}(ap,η',ε)
                                    EV = 0.0
                                    @inbounds for ηp in 1:n_η
                                        Vηp = lip(vGridAsset, Vnext[ηp, :], ap; on_extrap=:error)
                                        EV += mTrans_η[ieta, ηp] * Vηp
                                    end

                                return utility(c, l, param) + coef * EV
                
                            end
                            # 爬坡括峰（用标量 m0 热启动）
                            m = clamp(m0 - 2, 1, nA-1)
                            v = wvalue(vGridAsset[m])
                            while m < nA
                                vnext = wvalue(vGridAsset[m+1])
                                if vnext >= v
                                    m += 1
                                    v  = vnext
                                else
                                    break
                                end
                            end
                            #m0 = m
                            # 小区间 [xl,xr]
                            lidx = max(1, m-1)
                            ridx = min(nA, m+1)
                            xl = vGridAsset[lidx]
                            xr = vGridAsset[ridx]
                            # ======= 角点处理 =======
                            a1 = 0.0
                            if m == 1
                                v0 = wvalue(xl)
                                v1 = wvalue(xl + psi1)
                                if v1 >= v0
                                    a1 = GoldenSearch(wvalue, xl, xr)
                                else
                                    a1 = xl
                                end
                            elseif m == nA
                                v0 = wvalue(xr)
                                v1 = wvalue(xr - psi1)
                                if v0 >= v1
                                    a1 = xr
                                else
                                    a1 = GoldenSearch(wvalue, xl, xr)
                                end
                            else
                                a1 = GoldenSearch(wvalue, xl, xr)
                            end
                            #=
                            if l < 0.0
                                l = 0.0
                                c = (1.0 + (1.0-τ_k)*rseq[s]) * a0 + trseq[s] - (1.0+γ)*a1
                            end
                            =#
                            A_policy[s, ieps, ieta, ia] = a1
                            V[s, ieps, ieta, ia] = wvalue(a1)
                            Y = (1.0 + (1.0-τ_k)*rseq[s]) * a0 +
                                trseq[s] +
                                (1.0 - τ_w_seq[s] - τ_b_seq[s]) * w0 -
                                (1.0 + γ) * a1

                            c = λ * Y
                            if c <= 0.0
                                c = 0.0
                                l = 0.0
                            else
                                l = 1.0 - c / ((1.0 - τ_w_seq[s] - τ_b_seq[s]) * w0) * (1.0-λ)/λ
                                if l < 0.0
                                    l = 0.0
                                    c = (1.0 + (1.0-τ_k)*rseq[s]) * a0 + trseq[s] - (1.0 + γ) * a1
                                    if c <= 0.0
                                        c = 0.0
                                    end
                                end
                            end
                            C_policy[s, ieps, ieta, ia] = c
                            L_policy[s, ieps, ieta, ia] = l
                            # (可选) 检查 V_s(·) 单调
                            
                            if ia > 1 && V[s, ieps, ieta, ia] < V[s, ieps, ieta, ia-1]
                                error("V not monotone: s=$s, ieps=$ieps, ieta=$ieta")
                            end
                            
                        end






                    end







                end
           end


    return A_policy, C_policy, L_policy, V     
end

function wage(kbar, Lbar, param)
    @unpack α = param
    w = (1.0 - α)*kbar^α * Lbar^(-α)
    return w 
end

function interest(kbar, Lbar, param)
    @unpack α, δ = param 
    r = α * kbar^(α - 1.0 ) * Lbar^(-α) - δ
    return r 
end
function utility(c,l,param)
    @unpack σ, λ = param 
        # 安全检查：消费必须为正，休闲必须在 [0,1] 之间
    if c < 0.0 || l < 0.0 || l > 1.0
        return -1e12  # 返回一个非常大的负数
    end
    if σ == 1.0
       u = λ*log(c) + (1.0 - λ)*log(1.0 - l) 
    else
        u = (c^λ * (1.0 - l)^(1.0 - λ))^(1.0 - σ) / (1.0 - σ)
    end
    return u
end

function ginid(x::AbstractVector, g::AbstractVector)
    ng = length(x)
    ng == length(g) || throw(DimensionMismatch("x and g must have the same length"))

    # ensure 1D vectors (Julia vectors already are 1D; this is mostly for safety)
    xvec = collect(x)
    gvec = collect(g)

    # MATLAB: x = max([x zeros(ng,1)],[],2)  => clamp negatives to 0
    xvec = max.(xvec, 0)

    xmean = sum(xvec .* gvec)
    xmean > 0 || throw(DomainError(xmean, "mean of x under weights g is nonpositive; Gini undefined"))

    # MATLAB: y=[x g]; y=sort(y,1);  => sort by x (and permute g accordingly)
    p = sortperm(xvec)           # indices that sort x ascending
    x_sorted = xvec[p]
    g_sorted = gvec[p]

    f = zeros(Float64, ng)       # accumulated "income share" (as in the MATLAB code)
    f[1] = g_sorted[1] * x_sorted[1] / xmean

    gini = 1 - f[1] * g_sorted[1]
    @inbounds for i in 2:ng
        f[i] = f[i-1] + g_sorted[i] * x_sorted[i] / xmean
        gini -= (f[i] + f[i-1]) * g_sorted[i]
    end

    return gini
end

function getssvalues(kbar,Lbar,average_hours,trbar0,τ_b,τ_w0,sp1_init,
                    mass_init, g_n_init, age_effciency, param)
    # Calculate steady-state values based on inputs
    @unpack θ, τ_w, τ_k, n_ϵ, n_η, nA, nages, n_income,
        vGridAsset, vGridIncome, income_max, income_min, τ_k,
        vGrids_ϵ, vGrids_η, n, mTrans_η, periods_working, a_max,
        a_min, G_Y  = param
    
    wbar = wage(kbar, Lbar, param)
    rbar = interest(kbar, Lbar, param)
    Penbar = θ * (1.0 - τ_w0 - τ_b)*wbar * average_hours* ones(n_ϵ)
    wseq0 = ones(nages)*wbar
    rseq0 = ones(nages)*rbar
    Penseq0 = ones(nages)*Penbar'#矩阵
    trseq0 = trbar0 * ones(nages)
    τ_w_seq = τ_w0*ones(nages)
    τ_b_seq = τ_b *ones(nages)
    sp1 = sp1_init
    mass = mass_init
    A_policy, C_policy, L_policy, V = HH_decision(wseq0, rseq0, Penseq0, trseq0, τ_w_seq, τ_b_seq, sp1, age_effciency, param)
    ga = zeros(nages, n_ϵ, n_η, nA)
    measure_tol = 0.0 
    # age-profiles
	agen=zeros(nages,2);		# assets workers
	cgen=zeros(nages,2);		# consumption epsilon=1,2
	lgen=zeros(nages);		# working hours workers
	lgeneps=zeros(nages,2);		# working hours, for each epsilon
	hoursnew=0;
	bigl=0;					# aggregate labor 
	biga=0;					# aggregate assets
	bigtax=0;				# aggregate taxes
	bigcontrib=0;			# aggregate contributions to pension system
	bigbequests=0;			# aggregate accidental bequests
	bigpensions=0;			# total pensions
	xbar0=0;				# average contributions to pension system
	bigc=0;					# aggregate consumption

    # distribution functions for Gini coefficients
	fa=zeros(nA);         # distribution of wealth a
	fnetincome=zeros(n_income);    # net income
	fgrossincome=zeros(n_income);	# gross income
	fwage=zeros(n_income);			# wage income
	fpension=zeros(n_income);		#  pension income
	fconsumption=zeros(n_income);	# consumption

    #initialization at age 1
	# equal distribution of abilities
	ga[1,1,1,1]=1/2*1/2*mass[1];
	ga[1,2,1,1]=1/2*1/2*mass[1];
	ga[1,1,2,1]=1/2*1/2*mass[1];
	ga[1,2,2,1]=1/2*1/2*mass[1];
	#trbar0 = trbar
	g_n = g_n_init
        for it=1:1:nages-1
            for ia=1:1:nA   # asset holding in period t 
                asset0=vGridAsset[ia];
                for ieps=1:1:n_ϵ
                    for ieta=1:1:n_η
                        measure=ga[it,ieps,ieta,ia];    # measure of household it,ieta,ieps,ia
                        
                                
                        c=C_policy[it,ieps,ieta,ia];
                        bigc=bigc+c*measure;
                        
                        if c<=0
                            fconsumption[1]=fconsumption[1]+measure;
                        elseif c>=income_max
                            fconsumption[n_income]=fconsumption[n_income]+measure;
                        else
                            in1=sum(c.>vGridIncome)
                            lambda=(vGridIncome[in1+1]-c) / (vGridIncome[in1+1]-vGridIncome[in1] );
                            fconsumption[in1]=fconsumption[in1]+lambda*measure;
                            fconsumption[in1+1]=fconsumption[in1+1]+(1-lambda)*measure;
                        end
                        
                        
                        fa[ia]=fa[ia]+measure;	# wealth distribution
                        if it<=periods_working
                            l0=L_policy[it,ieps,ieta,ia];	# working hours at age it	
                        end
                        
                        a1=A_policy[it,ieps,ieta,ia]		# optimal next-period assets
                        
                        cgen[it,ieps]=cgen[it,ieps]+measure*c;
                        agen[it,ieps]=agen[it,ieps]+measure*asset0;
                        bigtax=bigtax+τ_k*rbar*asset0*measure;		# interest rate tax worker
                        measure_tol=measure_tol+measure;
                            
                            
                        if it<=periods_working
                            lgen[it]=lgen[it]+measure*l0;
                            lgeneps[it,ieps]=lgeneps[it,ieps]+l0*measure;
                            bigl=bigl+l0*vGrids_η[ieta]*vGrids_ϵ[ieps]*age_effciency[it]*measure;		# effective labor supply				
                            bigtax=bigtax+τ_w*l0*age_effciency[it]*vGrids_ϵ[ieps]*vGrids_η[ieta]*wbar*measure;		# wage income tax
                            bigcontrib=bigcontrib+τ_b*l0*age_effciency[it]*vGrids_ϵ[ieps]*vGrids_η[ieta]*wbar*measure;	# pension contributions			
                            wageincome=l0*age_effciency[it]*vGrids_ϵ[ieps]*vGrids_η[ieta]*wbar;
                            netincome=(1-τ_w-τ_b)*wageincome+trbar0+(1-τ_k)*rbar*asset0;
                            
                            # wage income distribution: only employed workers
                            if wageincome <=0
                                fwage[1]=fwage[1]+measure;
                            elseif wageincome>=income_max
                                fwage[n_income]=fwage[n_income]+measure;
                            else
                                in1=sum(wageincome.>vGridIncome);
                                lambda=(vGridIncome[in1+1]-wageincome) / (vGridIncome[in1+1]-vGridIncome[in1] );
                                fwage[in1]=fwage[in1]+lambda*measure;
                                fwage[in1+1]=fwage[in1+1]+(1-lambda)*measure;
                            end
                                
                            grossincome=wageincome+trbar0+rbar*asset0;
                        
                        else
                            bigpensions=bigpensions+Penbar[ieps]*measure;
                            pensionincome=Penbar[ieps];
                            grossincome=Penbar[ieps]+trbar0+rbar*asset0;
                            netincome=Penbar[ieps]+trbar0+(1-τ_k)*rbar*asset0;	
                        end
                            
                        
                        biga=biga+vGridAsset[ia]*measure;

                        # income distributions
                        
                        if grossincome<=0
                            fgrossincome[1]=fgrossincome[1]+measure;
                        elseif grossincome>=income_max
                            fgrossincome[n_income]=fgrossincome[n_income]+measure;
                        else
                            in1=sum(grossincome.>vGridIncome);
                            lambda=(vGridIncome[in1+1]-grossincome) / (vGridIncome[in1+1]-vGridIncome[in1] );
                            fgrossincome[in1]=fgrossincome[in1]+lambda*measure;
                            fgrossincome[in1+1]=fgrossincome[in1+1]+(1-lambda)*measure;
                        end

                        if netincome<=0
                            fnetincome[1]=fnetincome[1]+measure;
                        elseif netincome>=income_max
                            fnetincome[n_income]=fnetincome[n_income]+measure;
                        else
                            in1=sum(netincome.>vGridIncome);
                            lambda=(vGridIncome[in1+1]-netincome) / (vGridIncome[in1+1]-vGridIncome[in1] );
                            fnetincome[in1]=fnetincome[in1]+lambda*measure;
                            fnetincome[in1+1]=fnetincome[in1+1]+(1-lambda)*measure;
                        end

                        if it>periods_working
                            if pensionincome<=0
                                fpension[1]=fpension[1]+measure;
                            elseif pensionincome>=income_max
                                fpension[n_income]=fpension[n_income]+measure;
                            else
                                in1=sum(pensionincome .>vGridIncome);
                                lambda=(vGridIncome[in1+1]-pensionincome) / (vGridIncome[in1+1]-vGridIncome[in1] );
                                fpension[in1]=fpension[in1]+lambda*measure;
                                fpension[in1+1]=fpension[in1+1]+(1-lambda)*measure;
                            end
                        end
                                

                        # computation of next-periods distribution: linear interpolation 
                        #  for a(ia1)<= a1 <= a(ia1+1) 
                        
                        if a1<=a_min
                            ga[it+1,ieps,1,1]=ga[it+1,ieps,1,1]+mTrans_η[ieta,1]*sp1[it]/(1+g_n)*measure;
                            ga[it+1,ieps,2,1]=ga[it+1,ieps,2,1]+mTrans_η[ieta,2]*sp1[it]/(1+g_n)*measure;
                        elseif a1>=a_max
                            ga[it+1,ieps,1,nA]=ga[it+1,ieps,1,nA]+mTrans_η[ieta,1]*sp1[it]/(1+g_n)*measure;
                            ga[it+1,ieps,2,nA]=ga[it+1,ieps,2,nA]+mTrans_η[ieta,2]*sp1[it]/(1+g_n)*measure;	
                        else
                            ia1=sum(vGridAsset .< a1);
                            lambda1=(vGridAsset[ia1+1]-a1) / (vGridAsset[ia1+1]-vGridAsset[ia1] );
                            ga[it+1,ieps,1,ia1]=ga[it+1,ieps,1,ia1]	+lambda1 * mTrans_η[ieta,1]*sp1[it]/(1+g_n)*measure;
                            ga[it+1,ieps,2,ia1]=ga[it+1,ieps,2,ia1]	+lambda1 * mTrans_η[ieta,2]*sp1[it]/(1+g_n)*measure;
                            ga[it+1,ieps,1,ia1+1]=ga[it+1,ieps,1,ia1+1]	+ (1-lambda1) * mTrans_η[ieta,1]*sp1[it]/(1+g_n)*measure;
                            ga[it+1,ieps,2,ia1+1]=ga[it+1,ieps,2,ia1+1]+ (1-lambda1) * mTrans_η[ieta,2]*sp1[it]/(1+g_n)*measure;
                        end
                    end # ieta
                end # ieps
            end	# ia
        end # it
    # aggregation for the last period
        it=nages;
        for ia=1:1:nA
            asset0=vGridAsset[ia];
            for ieps= 1:1:n_ϵ
                for ieta=1:1:n_η
                    measure=ga[it,ieps,ieta,ia];	# measure of household it,ieta,ieps,ia,ix
                    measure_tol=measure_tol+measure;		
                    c=C_policy[it,ieps,ieta,ia];	
                    bigc=bigc+c*measure;
                        
                    if c<=0
                        fconsumption[1]=fconsumption[1]+measure;
                    elseif c>=income_max
                        fconsumption[n_income]=fconsumption[n_income]+measure;
                    else
                        in1=sum(c .>vGridIncome);
                        lambda=(vGridIncome[in1+1]-c) / (vGridIncome[in1+1]-vGridIncome[in1] );
                        fconsumption[in1]=fconsumption[in1]+lambda*measure;
                        fconsumption[in1+1]=fconsumption[in1+1]+(1-lambda)*measure;
                    end

                    fa[ia]=fa[ia]+measure;	# wealth distribution

                    grossincome=Penbar[ieps]+trbar0+rbar*asset0;
                    netincome=Penbar[ieps]+trbar0+(1-τ_k)*rbar*asset0;	
                    pensionincome=Penbar[ieps];
                    if grossincome<=0
                        fgrossincome[1]=fgrossincome[1]+measure;
                    elseif grossincome>=income_max
                        fgrossincome[n_income]=fgrossincome[n_income]+measure;
                    else
                        in1=sum(grossincome .>vGridIncome);
                        lambda=(vGridIncome[in1+1]-grossincome) / (vGridIncome[in1+1]-vGridIncome[in1] );
                        fgrossincome[in1]=fgrossincome[in1]+lambda*measure;
                        fgrossincome[in1+1]=fgrossincome[in1+1]+(1-lambda)*measure;
                    end

                    if netincome<=0
                        fnetincome[1]=fnetincome[1]+measure;
                    elseif netincome>=income_max
                        fnetincome[n_income]=fnetincome[n_income]+measure;
                    else
                        in1=sum(netincome .>vGridIncome);
                        lambda=(vGridIncome[in1+1]-netincome) / (vGridIncome[in1+1]-vGridIncome[in1] );
                        fnetincome[in1]=fnetincome[in1]+lambda*measure;
                        fnetincome[in1+1]=fnetincome[in1+1]+(1-lambda)*measure;
                    end
                                
                                
                    if pensionincome<=0
                        fpension[1]=fpension[1]+measure;
                    elseif pensionincome>=income_max
                        fpension[n_income]=fpension[n_income]+measure;
                    else
                        in1=sum(pensionincome .>vGridIncome);
                        lambda=(vGridIncome[in1+1]-pensionincome) / (vGridIncome[in1+1]-vGridIncome[in1] );
                        fpension[in1]=fpension[in1]+lambda*measure;
                        fpension[in1+1]=fpension[in1+1]+(1-lambda)*measure;
                    end
                                
                    cgen[it,ieps]=cgen[it,ieps]+measure*c;
                    agen[it,ieps]=agen[it,ieps]+measure*asset0;
                    bigtax=bigtax+τ_k*rbar*asset0*measure;		# interest rate tax 75-year old
                        
                    bigpensions=bigpensions+Penbar[ieps]*measure;
                    biga=biga+asset0*measure;
                        
                end	# ieta
            end	# ieps
        end # ia
        #println("aggregation done")
        # computation of the Ginis
	fwage=fwage/sum(fwage)     # normalization of mass of workers to one
	fpension=fpension/sum(fpension)
	giniwage=ginid(vGridIncome,fwage)
	ginigrossincome=ginid(vGridIncome,fgrossincome)
	gininetincome=ginid(vGridIncome,fnetincome)
	giniwealth=ginid(vGridAsset,fa)
	ginipension=ginid(vGridIncome,fpension)
	giniconsumption=ginid(vGridIncome,fconsumption)
	progressivityindex=1-ginipension/giniwage
	# total bequests: check if equilibrium condition is fine
	biga=sum(sum(agen))
	bigk=biga
	bigbequests=(1.0 .- sp1)'*agen
	bigbequests=sum(bigbequests')
	bigbequests=(1+(1-τ_k)*rbar)*bigbequests
	bigy=kbar^(α)*Lbar^(1-α);
	hoursnew=sum(lgen)/sum(mass[1:periods_working]);

	bigtax0=τ_k*rbar*biga+τ_w*wbar*bigl;
	welfare=1/4*(V[1,1,1,1]+V[1,1,2,1]+V[1,2,1,1]+V[1,2,2,1])
    taubnew = bigpensions/(wbar*bigl)


	#@unpack G_Y = param
	gbar = G_Y * bigy 
	trnew = bigtax0 + bigbequests - gbar
	tauwnew = τ_w
	y1 = bigk
    y2 = bigl
    y3 = hoursnew
    y4 = trnew
    y5 = taubnew
    y6 = tauwnew
    #y = [y1, y2, y3, y4, y5, y6, y7]
    return y1, y2, y3, y4, y5, y6
end

function obj_getss(x)
    kbar = x[1]
    Lbar = x[2]
    average_hours = x[3]
    trbar = x[4]
    τ_b = x[5]
    τ_w = x[6]
    y1, y2, y3, y4, y5, y6 = getssvalues(kbar,Lbar,average_hours,trbar,τ_b,τ_w,sp1_init,
                    mass_init, g_n_init, age_effciency, param)
    err1 = x[1] - y1
    err2 = x[2] - y2
    err3 = x[3] - y3
    err4 = x[4] - y4
    err5 = x[5] - y5
    err6 = x[6] - y6
    err = [err1, err2, err3, err4, err5, err6]
    return err 
endau
