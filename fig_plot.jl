using DelimitedFiles, StatsBase, Plots
using HDF5, SparseArrays, XLSX, Dates

### Theorectical curves
f(x)=0.5 * (1.0 .-x .+ sqrt.(8.0*l .+ (x .-1).^2*NC)/sqrt(NC) )
ff=collect(0.0:0.02:10)

### Smartseq3
data = readdlm("Smartseq3.HEK.cleanup.UMIcounts.txt",'\t')
bc = data[1,1:end-1]
counts = data[2:end,2:end]
info = readdlm("Smartseq3.HEK.cleanup.sample_annotation.txt",'\t')
bc2 = info[2:end,4]
cond = zeros(length(bc2))
for i = 1 : length(bc2)
    if info[i+1,5] == "no"
        if info[i+1,6] == 0
            cond[i] = 1
        else
            cond[i] = 2
        end
    else
        if info[i+1,6] == 0
            cond[i] = 3
        else
            cond[i] = 4
        end
    end
end
ind = [findfirst(x->x==bc2[i],bc) for i = 1 : length(bc2)]
count_cur = counts[:,ind]
M = zeros(size(count_cur,1),4)
FF = zeros(size(count_cur,1),4)
nc = zeros(4)
for i = 1 : 4
    flg = cond .== i
    nc[i] = sum(flg)
    count_filter = count_cur[:,flg]
    for j = 1 : size(M,1)
        M[j,i] = mean(count_filter[j,:])
        FF[j,i] = var(count_filter[j,:],corrected=false)/mean(count_filter[j,:])
    end
end
i=4
p1 = scatter(FF[:,i],M[:,i],xlims=(0.7,2),ylims=(0,0.4),markersize=1.5,markerstrokewidth=0,xlabel="FF",ylabel="Mean", title="Smartseq-3",legend=false)
l = 0
NC = nc[i]
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### VASA seq
data = readdlm("GSM5369497_E6.5-1_i1_total.UFICounts.tsv",'\t','\n')[2:end,2:end]
M = zeros(size(data,1))
FF = zeros(size(data,1))
for i = 1 : size(data,1)
    M[i] = mean(data[i,:])
    FF[i] = var(data[i,:],corrected=false)/mean(data[i,:])
end
p2 = scatter(FF,M,xlims=(0.66,2),ylims=(0,0.5),markersize=1.5,markerstrokewidth=0,xlabel="FF",ylabel="Mean", title="VASA-seq",legend=false)
l = 0
NC = size(data,2)
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### 10x Genomics
fd = h5open("pbmc_1k_v3_filtered_feature_bc_matrix.h5","r")
barcodes = read(fd["matrix/barcodes"])
data = read(fd["matrix/data"])
indices = read(fd["matrix/indices"]).+1
indptr = read(fd["matrix/indptr"]).+1
num_rows = maximum(indices)
num_cols = length(indptr) - 1 
sparse_matrix = SparseMatrixCSC(num_rows, num_cols, indptr, indices, data)
matrix = Matrix(sparse_matrix)
annotation = readdlm("cluster.csv",',')[2:end,:]
hg = annotation[:,2] .== 6
hg_matrix = matrix[:,hg]
M_hg = zeros(size(hg_matrix,1))
FF_hg = zeros(size(hg_matrix,1))
for  i = 1 : size(hg_matrix,1)
    M_hg[i] = mean(hg_matrix[i,:])
    FF_hg[i] = var(hg_matrix[i,:],corrected=false)/mean(hg_matrix[i,:])
end
p3 = scatter(FF_hg,M_hg,xlims=(0.65,2),ylims=(0,0.4),markersize=1.5,markerstrokewidth=0,xlabel="FF",ylabel="Mean", title="10x Genomics",legend=false)
NC = sum(hg)
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### SCRB seq
xf = XLSX.readxlsx("41564_2018_330_MOESM3_ESM.xlsx")
sh = xf["Table_S4"]
data = sh["A1:DDA7051"][4:end,2:end]
barcode = sh["A1:DDA7051"][3,2:end]
sh2 = xf["Table_S3"]
cond = sh2["A1:M2824"][4:end-13,10]
u_cond = unique(cond)
#bc2 = sh2["A1:M2824"][4:end-13,1]
#ind = [findfirst(x->x==bc2[i],barcode) for i = 1 : length(bc2)]
#iszero.(sum(ind.-collect(1:length(bc2))))
M = zeros(size(data,1),length(u_cond))
FF = zeros(size(data,1),length(u_cond))
nc = zeros(length(u_cond))
for i = 1 : length(u_cond)
    flg = cond .==  u_cond[i]
    ctemp = data[:,flg]
    nc[i] = sum(flg)
    for j = 1 : size(data,1)
        M[j,i] = mean(ctemp[j,:])
        FF[j,i] = var(ctemp[j,:],corrected=false)/mean(ctemp[j,:])
    end
    print(i)
    print("\n")
end
i = 7
NC = nc[i]
u_cond[i]
p4 = scatter(FF[:,i],M[:,i],xlims=(0.64,2),ylims=(0,0.5),markersize=1.5,markerstrokewidth=0,xlabel="FF",ylabel="Mean", title="SCRB-seq",legend=false)
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### MERFISH
data = readdlm("pnas.1912459116.sd12.csv",',')
bc = data[1,2:end]
count = data[2:end,2:end]
flg = occursin.("B3",bc)
nc = sum(flg)
count_cur = count[:,flg]
M = zeros(size(count_cur,1))
FF = zeros(size(count_cur,1))
for i = 1 : size(count_cur,1)
    M[i] = mean(count_cur[i,:])
    FF[i] = var(count_cur[i,:],corrected=false)/mean(count_cur[i,:])
end
NC = nc
p5 = scatter(FF,M,xlims=(0.74,1.4),ylims=(0,0.3),markersize=1.5,markerstrokewidth=0,xlabel="FF",ylabel="Mean", title="MERFISH",legend=false)
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)


### FLASH-seq
data = readdlm("HEK_FS_lowAmplification_250K.csv",',')[2:end,2:end]
M = zeros(size(data,1))
FF = zeros(size(data,1))
for i = 1 : size(data,1)
    M[i] = mean(data[i,:])
    FF[i] = var(data[i,:],corrected=false)/mean(data[i,:])
end
NC = size(data,2)
p6 = scatter(FF,M,xlims=(0.75,2),ylims=(0,0.3),markersize=1.5,markerstrokewidth=0,xlabel="FF",ylabel="Mean", title="FLASH-seq",legend=false)
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### Integrate the first 6 plots
plot(p1,p2,p3,p4,p5,p6,layout=(2,3))

### Birds in North America
data = readdlm("PFW_all_2021_2023_June2023_Public.csv",',')
flg = data[:,10].==2022
datac = data[flg,:]
loc = unique(datac[:,1])
name = unique(datac[:,12])
function build_index(vector)
    index_map = Dict{Any, Int}()
    for (index, value) in enumerate(vector)
        index_map[value] = index
    end
    return index_map
end
function search_index(index_map, target)
    return get(index_map, target, -1)  # returns -1 if the target is not found
end
index_map = build_index(loc)
counts = zeros(length(loc),length(name))
for i = 1 : size(datac,1)
    ind1 = search_index(index_map, datac[i,1])
    ind2 = findfirst(x->x==datac[i,12],name)
    counts[ind1,ind2] = counts[ind1,ind2]+1
end
M = zeros(size(counts,1))
FF = zeros(size(counts,1))
for i = 1 : size(counts,1)
    M[i] = mean(counts[i,:])
    FF[i] = var(counts[i,:],corrected=false)/M[i]
end
p7 = scatter(FF,M,xlims=(0.85,2),ylims=(0,0.1),markersize=3,markerstrokewidth=0,xlabel="FF",ylabel="Mean", title="Birds in North America",legend=false)
NC = size(counts,2)
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### High-energy physics paper citation
### https://snap.stanford.edu/data/cit-HepTh.html
data = readdlm("Cit-HepTh.txt",'\t')[5:end,:]
data2 = readdlm("Cit-HepTh-dates.txt",'\t')[2:end,:]
date = data2[:,2]
mon = [month(Date(date[i],"yyyy-mm-dd")) for i = 1 : length(date)]
yr = [year(Date(date[i],"yyyy-mm-dd")) for i = 1 : length(date)]
yr_int = 1995
paper = data2[yr.==yr_int,1]
ng = length(paper)
nc = (2002-yr_int)*12
counts = zeros(ng,nc)
for i = 1 : ng
    paperid = paper[i]
    flg = data[:,2] .== paperid
    citation = data[flg,1]
    for j = 1 : length(citation)
        citation_ind = findfirst(x->x==citation[j],data2[:,1])
        if !isnothing(citation_ind)
            if yr[citation_ind] > yr_int 
                yr_temp = yr[citation_ind] - yr_int-1
                mon_temp = mon[citation_ind]
                cellid = yr_temp*12+mon_temp
                counts[i,cellid] += 1
            end
        end
    end
    if mod(i,50) == 0 
        print(i)
        print("\n")
    end
end
M = zeros(ng)
FF = zeros(ng)
for  i = 1 : ng
    M[i] = mean(counts[i,:])
    FF[i] = var(counts[i,:],corrected=false)/mean(counts[i,:])
end
p8 = scatter(FF,M,xlabel="FF",ylabel="Mean",title="HEP paper citation",xlims=(0.7,2),ylims=(0,0.5),markersize=3,markerstrokewidth=0,legend=false)
NC = nc
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### Colorado Butterfly
data_1 = readdlm("Boulder_Abundance.txt",'\t')
data = data_1[2:end,3:end]
M = zeros(size(data,1))
FF = zeros(size(data,1))
for i = 1 : size(data,1)
    M[i] = mean(data[i,:])
    FF[i] = var(data[i,:],corrected=false)/mean(data[i,:])
end
unique_matrix=hcat(FF,M)
p9 = scatter(FF,M,xlims=(0.6,4),ylims=(0,0.4),markersize=3,markerstrokewidth=0,xlabel="FF",ylabel="Mean", title="Colorado Butterfly",legend=false)
NC=size(data,2)
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### E-commerce in UK
#= #Pre-processing
data = readdlm("/media/dell/data/pattern/data.csv",',')
stock = data[2:end,3]
stock_u = unique(stock)
customer = data[2:end,5]
quant = data[2:end,4]
flg = quant .>0
country = data[2:end,8]
flg2 = country .== "United Kingdom"
sum(flg2)

stock_c = stock[flg.&flg2]
customer_c = customer[flg.&flg2]
country_c = country[flg.&flg2]
quant_c = quant[flg.&flg2]
customer_u = unique(customer_c)
nc = length(customer_u)
ng = length(stock_u)

counts = zeros(ng,nc)
#i=1
for i = 1 : length(quant_c)
    ind_c = findfirst(x->x==customer_c[i],customer_u)
    ind_s = findfirst(x->x==stock_c[i],stock_u)
    counts[ind_s,ind_c] = quant_c[i]
    if mod(i,1000) == 0
        print(i)
        print("\n")
    end
end
#writedlm("e-comm.csv",counts,',')
=#
counts = readdlm("e-comm.csv",',')
ng = size(counts,1)
M = zeros(ng)
FF = zeros(ng)
for  i = 1 : ng
    M[i] = mean(counts[i,:])
    FF[i] = var(counts[i,:],corrected=false)/mean(counts[i,:])
end
NC = size(counts,2)
p10 = scatter(FF,M,xlims=(0.9,2),ylims=(0,0.006),xlabel="FF",ylabel="Mean",title="E-Commerce data in UK", markersize=3,markerstrokewidth=0,legend=false)
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### Disease outbreak in US
data=readdlm("cdc.csv",',')[2:end,2:end]
M = zeros(size(data,1))
FF = zeros(size(data,1))
for i = 1 : size(data,1)
    M[i] = mean(data[i,:])
    FF[i] = var(data[i,:],corrected=false)/mean(data[i,:])
end
p11=scatter(FF,M,xlims=(0.6,2.5),ylims=(0,0.6),xlabel="FF",ylabel="Mean",title="Disease outbreak in US",markersize=3,markerstrokewidth=0,legend=false)
NC=56
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### Job posting in US
data = readdlm("postings.csv",',')
job = data[2:end,3]
u_job = unique(job)
region = data[2:end,7]
states = readdlm("us-states-territories.csv",',')[2:57,2:3]
counts = zeros(length(u_job),size(states,1))
for i = 1 : length(job)
    inds = findfirst(x->occursin(x,region[i]),states[:,1])
    inds2 = findfirst(x->occursin(x,region[i]),states[:,2])
    if !isnothing(inds) | !isnothing(inds2)
        if !isnothing(inds)
            indr = inds
        else
            indr = inds2
        end
        ind = findfirst(x->x==job[i],u_job)
        counts[ind, indr] = counts[ind,indr] + 1
    end
end
M = zeros(size(counts,1))
FF = zeros(size(counts,1))
for i = 1 : size(counts,1)
    M[i] = mean(counts[i,:])
    FF[i] = var(counts[i,:],corrected=false)/M[i]
end
p12=scatter(FF,M,xlims=(0.6,5),ylims=(0,0.2),title="LinkedIn jobs in US",markersize=3,markerstrokewidth=0,legend=false)
NC = size(counts,2)
l = 0
plot!(ff,f.(ff),linewidth=2)
l = 1
plot!(ff,f.(ff),linewidth=2)
l = 2
plot!(ff,f.(ff),linewidth=2)
l = 3
plot!(ff,f.(ff),linewidth=2)
l = 4
plot!(ff,f.(ff),linewidth=2)
l = 5
plot!(ff,f.(ff),linewidth=2)

### Integrate the second 6 plots
plot(p7,p8,p9,p10,p11,p12,layout=(2,3))
plot!(size=(2000,1200))