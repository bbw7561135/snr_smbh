rm *.mp4
ffmpeg -i log10_number_density_%d.png snr_smbh_log10_number_density.mp4
ffmpeg -i log10_temperature_%d.png snr_smbh_log10_temperature.mp4
ffmpeg -i x_velocity_%d.png snr_smbh_x_velocity.mp4
ffmpeg -i y_velocity_%d.png snr_smbh_y_velocity.mp4
