rm *.webm
ffmpeg -i log10_number_density_%d.png -codec:v libvpx "snr_smbh_log10_number_density+`git rev-parse HEAD`.webm"
ffmpeg -i log10_temperature_%d.png -codec:v libvpx "snr_smbh_log10_temperature_`git rev-parse HEAD`.webm"
ffmpeg -i x_velocity_%d.png -codec:v libvpx "snr_smbh_x_velocity_`git rev-parse HEAD`.webm"
ffmpeg -i y_velocity_%d.png -codec:v libvpx "snr_smbh_y_velocity_`git rev-parse HEAD`.webm"
