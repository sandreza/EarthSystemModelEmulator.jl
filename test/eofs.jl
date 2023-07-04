using EarthSystemModelEmulator
using NCDatasets
using GLMakie
using Statistics
using Printf, LinearAlgebra

data_directory = "/Users/andresouza/Desktop/CMIP6/CESM/SurfaceTemperatureHistorical/"
files = readdir(data_directory)

file_path = data_directory * files[1]
ds = Dataset(file_path)
t1_old = ds["ts"][:,:,:]
t1_r = extrema(ds["ts"][:,:,:]) 
lon = ds["lon"][:]
lat = ds["lat"][:]
times_old = ds["time"][:]
close(ds)
metric = reshape(cos.(lat / 180 * π), (1, length(lat)))
metric /= mean(metric)

t1 = t1_old .- mean(t1_old, dims=3)
##
s1, s2, s3 = size(t1)
rt = reshape(t1 , (s1 * s2, s3))
U, S, V = svd(rt) 
##
fig = Figure() 
scatter(log10.(S))
S[1:20] ./ sum(S) .* 100
shift = 0
for i in 1:6
    v1 = reshape(U[:, i+shift], (s1, s2)) 
    ax = Axis(fig[i,1])
    heatmap!(ax, v1, colormap = :balance)
    # ax2 = Axis(fig[i,2])
    # lines!(ax2, log10.(S))
    ax3 = Axis(fig[i,2])
    lines!(ax3, V[:,i+shift])
end
display(fig)

##
file_path = data_directory * files[end-1]
ds = Dataset(file_path)
t2_old = ds["ts"][:,:,:]
t2_r = extrema(ds["ts"][:,:,:]) 
lon = ds["lon"][:]
lat = ds["lat"][:]
times_old = ds["time"][:]
close(ds)
metric = reshape(cos.(lat / 180 * π), (1, length(lat)))
metric /= mean(metric)

t2 = t2_old .- mean(t2_old, dims=3)
##
s1, s2, s3 = size(t2)
rt = reshape(t2, (s1 * s2, s3))
U2, S2, V2 = svd(rt) 
##
fig = Figure() 
scatter(log10.(S))
S[1:20] ./ sum(S) .* 100
shift = 0
for i in 1:6
    v1 = reshape(U2[:, i+shift], (s1, s2))
    ax = Axis(fig[i,1])
    heatmap!(ax, v1, colormap = :balance)
    # ax2 = Axis(fig[i,2])
    # lines!(ax2, log10.(S))
    ax3 = Axis(fig[i,2])
    lines!(ax3, V2[:,i+shift])
end
display(fig)
##
fig = Figure() 
scatter(log10.(S))
S[1:20] ./ sum(S) .* 100
shift = 0
for i in 1:6
    v1 = reshape(U2[:, i+shift], (s1, s2))
    ax = Axis(fig[i,1])
    heatmap!(ax, v1, colormap = :balance)
    v1 = reshape(U[:, i+shift], (s1, s2))
    ax = Axis(fig[i,2])
    heatmap!(ax, v1, colormap = :balance)
end
display(fig)