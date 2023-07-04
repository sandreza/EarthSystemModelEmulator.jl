using EarthSystemModelEmulator
using NCDatasets
using GLMakie
using Statistics
using Printf

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

file_path = data_directory * files[end-1]
ds = Dataset(file_path)
t1_new = ds["ts"][:,:,:]
times_new = ds["time"][:]
close(ds)
##
GLMakie.activate!(inline=false)
fig = Figure(resolution = (3000,1000))

metric = reshape(cos.(lat / 180 * π), (1, length(lat)))
metric /= mean(metric)
colormap = :thermometer

sliderrange = minimum([size(t1_old)[3],minimum(size(t1_new)[3])])
slider = Slider(fig[2, 1:3], range = 1:sliderrange, value = 1)
obs = slider.value
ax = Axis(fig[1, 1]; title = @lift(string(times_old[$obs])[1:10]))
heatmap!(ax, lon, lat, @lift(t1_old[:, :, $obs]), colorrange = t1_r, colormap = colormap)
display(fig)

ax2 = Axis(fig[1, 2]; title =  @lift(string(times_new[$obs])[1:10]))
heatmap!(ax2, lon, lat, @lift(t1_new[:, :, $obs]), colorrange = t1_r, colormap = colormap)
display(fig)

DT = @lift(t1_new[:, :, $obs]-t1_old[:, :, $obs])
DT_bar = @lift(mean(metric .* $DT))
DT_std = @lift(sqrt.( mean(metric .* ($DT .^2)) .- mean(metric .* $DT)^2 ))
ΔT = @lift(maximum(abs.($DT /2)))
titlestring = @lift("Difference: ΔT global mean = " * @sprintf("%.1f", $DT_bar) * ", ΔT global std = " * @sprintf("%.1f", $DT_std) * ", ΔT local colorbar =" *  @sprintf("%.1f", $ΔT) )
ax3 = Axis(fig[1, 3]; title = titlestring)

colorrange = @lift( (-$ΔT, $ΔT) )
heatmap!(ax3, lon, lat, DT, colorrange = colorrange, colormap = :balance)
display(fig)
