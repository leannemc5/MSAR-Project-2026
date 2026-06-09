wavelengths = np.linspace(100e-9, 3000e-9, 5000)
temperatures = [3000, 5000, 7000, 10000]

# Wien's displacement constant (m·K)
b = 2.897771955e-3

for T in temperatures:
    # wein
    predicted_peak = b / T

    # Optional: Find the closest matching wavelength value in your grid array
    closest_grid_index = np.abs(wavelengths - predicted_peak).argmin()
    closest_grid_val = wavelengths[closest_grid_index]
    
    print(f"Temperature: {T:5d} K | Exact Peak: { predicted_peak*1e9:6.1f} nm | Grid Match: {closest_grid_val*1e9:6.1f} nm")
