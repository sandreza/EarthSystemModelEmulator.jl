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

t1 = t1_old .- mean(t1_old .* metric, dims=3)
##
s1, s2, s3 = size(t1)
rt = reshape(t1 .* sqrt.(metric), (s1 * s2, s3))
U, S, V = svd(rt) 
##
fig = Figure() 
scatter(log10.(S))
S[1:20] ./ sum(S) .* 100
shift = 1
for i in 1:6
    v1 = reshape(U[:, i+shift], (s1, s2)) ./ sqrt.(metric)
    ax = Axis(fig[i,1])
    heatmap!(ax, v1, colormap = :balance)
    # ax2 = Axis(fig[i,2])
    # lines!(ax2, log10.(S))
    ax3 = Axis(fig[i,2])
    lines!(ax3, V[:,i+shift])
end
display(fig)