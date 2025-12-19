#=
using Parameters, Plots
using XLSX, Statistics, NLsolve,LinearAlgebra
include("module FunctionApprox.jl")
include("./Numerical_method.jl")
include("Numerical Optim.jl")
include("./aux/auxf.jl")
using .FunctionApprox
using .MyNumerical_method
using .Numerical_Optim
=#

function SettingPar(;
                    n = 0.008,
                    β = 1.011,
                    σ = 2.0,
                    λ = 0.31,
                    α = 0.36,
                    δ = 0.08,
                    γ = 0.02,
                    ϵ1 = 0.57,
                    ϵ2 = 1.43,
                    G_Y = 0.195,
                    τ_w = 0.248,
                    τ_k = 0.429,
                    θ = 0.5,
                    year = 2015,
                    nA = 150,
                    a_min = 0.0,
                    a_max = 10.0,
                    l_max = 0.9,
                    l_min = 0.1,
                    n_income = 1000,
                    income_max = 4.0,
                    income_min = 0.0)

    vGridAsset = collect(range(a_min,a_max,nA))
    vGrids_η = [0.727, 1.273]
    mTrans_η = [0.98 0.02; 0.02 0.98]
    vGrids_ϵ = [ϵ1 , ϵ2]
    inv_dist = [1/4, 1/4, 1/4, 1/4]
    periods_retire = 25
    periods_working = 45
    nages = periods_retire + periods_working
    year_final_policy = year + 45
    vGridIncome = collect(range(income_min, income_max, n_income))
    return (n = n, β = β, σ = σ, λ = λ, α = α, δ = δ, γ = γ, vGrids_ϵ = vGrids_ϵ, vGrids_η = vGrids_η,
            mTrans_η = mTrans_η, G_Y = G_Y, τ_w = τ_w, τ_k = τ_k, θ = θ, inv_dist = inv_dist,
            periods_retire = periods_retire, periods_working = periods_working, nages = nages, year = year,
            year_final_policy = year_final_policy, vGridAsset = vGridAsset, l_min = l_min, l_max = l_max,
            n_η = 2, n_ϵ = 2, nA = nA, a_max = a_max, a_min = a_min,vGridIncome = vGridIncome, n_income = n_income,
            income_max = income_max, income_min = income_min)
end

# 设置参数
param = SettingPar()

data = XLSX.readdata("survival_probs_US.xlsx", "Data", "F22:AI37")
survivalprobs = reverse(Matrix{Float64}(data), dims=1)

popgrowth=XLSX.readdata("survival_probs_US.xlsx",2,"F41:AI43")
popgrowth = Matrix{Float64}(popgrowth)
timespan = collect(Int64, range(1950, 2095, 30))
nrate = [timespan popgrowth']

# 计算移动平均 
periods_moving = 4
@unpack year = param
year0 = Int((year - 1950)/5 + 1)
nage1 = 15

case_UNscen = 1
if (year - 1950)%5 != 0
    error("年份year必须是5的倍数")
else
    popgrowth = nrate[year0 - periods_moving + 1:year0, case_UNscen + 1]
    popgrowth = mean(popgrowth)/100
    sp = survivalprobs[1:nage1, year0-periods_moving+1:year0]
    sp = vec(mean(sp,dims = 2))
end
println(year)
println(popgrowth)


age = collect(range(20, 20+5*(nage1-1), nage1))
nage0 = 5*(nage1 - 1)+1
age1 = collect(range(20, 20+5*(nage1-1), nage0))
y2 = cspline(age,sp,1)
sp1 = splint(age, sp, y2, 1, age1)
sp1 = sp1.^(1/5)
p1 = plot(age, sp; label = false)
title!("5-annual survival probability")
p2 = plot(age1, sp1;label = false)
title!("Annual survival probability")
#plot(p1, p2, layout=(1,2), size = (900,500), linewidth = 1.5, xlims = (20, 90))


@unpack nages = param
sp1 = sp1[1:nages]
sp1_init = sp1
g_n_init = popgrowth
age_effciency = XLSX.readdata("efficiency_profile.xlsx",1,"A1:A45")
age_effciency = age_effciency./mean(age_effciency)
age2 = collect(range(20,64,45))
p3 = plot(age2, age_effciency; label = false)
xlabel!("age")
title!("Productivity")


mass = ones(nages)
for t in 2:nages
    mass[t] = mass[t-1]*sp1[t-1] /(1 + popgrowth)# sp1 按年度计算的存活概率
end
mass = mass./sum(mass)
mass_init = mass
p4 = plot(age1[1:nages], mass; label = false, xlims = (20, 90), linewidth = 1.5, title = "mass of saving Prob")


year1 = 2015		
year0 = Int((year1-1950)/5+1)
sptotal = sp1_init[1:nages]
gntotal = copy(g_n_init)

while year0<size(nrate, 1)
    global year0, popgrowth,sp, sp1, sptotal, gntotal, y2
	year0=year0+1
	if year0==round(year0)
		popgrowth=nrate[year0-periods_moving+1:year0,case_UNscen+1]
		popgrowth=mean(popgrowth,dims = 1)
		popgrowth=popgrowth/100	
		sp = survivalprobs[1:nage1, year0-periods_moving+1:year0]
		sp = vec(mean(sp,dims = 2))
    else
		error("year1 must be a multiple of 5")
		
    end
    y2 = cspline(age,sp,1)
    sp1 = splint(age, sp, y2, 1, age1)
    sp1 = sp1.^(1/5)         
	sptotal = hcat(sptotal, sp1[1:nages])
	gntotal = vcat(gntotal, popgrowth)
end

nsp = size(sptotal,2)
year_t = collect(range(year1, year1 + 5*(nsp - 1), nsp))
rows_sp_all = (nsp-1)*5+1
year_t_1 = collect(range(year1, year1+5*(nsp - 1), rows_sp_all))
g_n_all = zeros(rows_sp_all)
sp_all = zeros(nages, rows_sp_all)

sp_all[1:nages,rows_sp_all] = sptotal[1:nages,nsp]
g_n_all[rows_sp_all] = gntotal[nsp]

for i=1:(nsp-1)    
	g_n_all[(i-1)*5+1] = gntotal[i]
	g_n_all[(i-1)*5+2] = 4/5 *gntotal[i] + 1/5 * gntotal[i+1]
	g_n_all[(i-1)*5+3] = 3/5 * gntotal[i] + 2/5 * gntotal[i+1]
	g_n_all[(i-1)*5+4] = 2/5 * gntotal[i] + 3/5 *gntotal[i+1]
	g_n_all[(i-1)*5+5] = 1/5 * gntotal[i]+ 4/5 *gntotal[i+1]
	
	sp_all[1:nages,(i-1)*5+1] = sptotal[1:nages,i]
	sp_all[1:nages,(i-1)*5+2] = 4/5 * sptotal[1:nages,i] + 1/5 * sptotal[1:nages,i+1]
	sp_all[1:nages,(i-1)*5+3] = 3/5 * sptotal[1:nages,i] + 2/5 * sptotal[1:nages,i+1]
	sp_all[1:nages,(i-1)*5+4] = 2/5 * sptotal[1:nages,i] + 3/5 * sptotal[1:nages,i+1]
	sp_all[1:nages,(i-1)*5+5] = 1/5 * sptotal[1:nages,i] + 4/5 * sptotal[1:nages,i+1]
	
end
ntrans = 250
massvec = zeros(nages,ntrans)			# mass of the generation j in period tp
massvec0 = copy(massvec)					# just to check if the mass is equal to one during each transition period
massvec1 = zeros(ntrans)
mass0 = copy(mass_init)
massvec[1:nages,1] = mass0
mass1 = zeros(nages)
function gnproc(periodx)
	if periodx>rows_sp_all
		periodx = rows_sp_all		# all survival probs are equal to those in the year 2095
	elseif periodx<1
		periodx=1
    end
	y = g_n_all[periodx]
    return y
end	
function surviveprob(agex,periodx)
# computes the survival probability at age agex in period periodx;

	if periodx>rows_sp_all
		periodx=rows_sp_all		# all survival probs are equal to those in the year 2095
	elseif periodx<1
		periodx=1;
    end
	y=sp_all[agex,periodx]
    return y
end

for i=1:ntrans-1	
	mass1[1]=mass0[1]*(1+gnproc(i+1))
    for j = 2:nages
		mass1[j] = mass0[j-1]*surviveprob(j-1,i)
    end
	mass0 .= mass1
	massvec[1:nages,i+1] .= mass1
end

massyear = sum(massvec,dims=2)
dependency_ratio = zeros(ntrans)
@unpack periods_working = param
for i=1:ntrans
	dependency_ratio[i]= sum(massvec[periods_working+1:nages,i])/sum(massvec[1:periods_working,i])
end
periods_trans = collect(range(2015,2015+ntrans-1,ntrans))
p5 = plot(periods_trans, dependency_ratio; label = false)
title!("dependency_ratio")
xlims!(periods_trans[1], periods_trans[end])

plot(p1, p2,p3,p4, p5, layout=(3,2), size = (1000,1200), linewidth = 1.5)



#稳态计算的初始化
@unpack  α, β, δ, G_Y, τ_w = param
average_hours = 0.298
Lbar = 0.27463
kbar = 1.10991 
#Lbar = average_hours * sum(mass_init[1:periods_working])
#kbar = (α /(1/β -1+δ) )^(1/(1-α))*Lbar#使用最优增长模型的稳态作为初值
ybar = kbar^α * Lbar^(1-α)
Gbar = G_Y * ybar
trbar = 0.02293
τ_b = 0.08532
#y1, y2, y3, y4, y5, y6 = getssvalues(kbar,Lbar,average_hours,trbar,τ_b,τ_w,sp1_init,
 #                   mass_init, g_n_init, age_effciency, param)
#println(y)
x_init = [kbar, Lbar, average_hours, trbar, τ_b, τ_w]
sol = nlsolve(obj_getss,x_init)




