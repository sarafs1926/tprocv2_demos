import numpy as np
import h5py
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import math
import os
import datetime

def calculate_qubit_temperature(frequency_mhz, ground_state_population, excited_state_population):
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck's constant in J·s
    frequency_hz = frequency_mhz * 1e6
    # T = (h * frequency_hz) / (k_B * np.log(ground_state_population / excited_state_population))
    # Check for invalid populations
    if excited_state_population <= 0 or ground_state_population <= 0:  # if one of them is zero can't calculate the temp
        print("Warning: Invalid population values encountered (<= 0). Skipping this dataset.")
        return None

    ratio = ground_state_population / excited_state_population
    if ratio <= 1:  # denominator would become zero at Pg=Pe
        print(f"Warning: Non-physical ratio (P_g/P_e = {ratio:.3f} <= 1) encountered. Skipping this dataset.")
        return None

    # If valid, calculate the temperature
    T = (h * frequency_hz) / (k_B * np.log(ratio))
    return T


def fit_double_gaussian_with_full_coverage(iq_data):
    gmm = GaussianMixture(n_components=2)
    gmm.fit(iq_data.reshape(-1, 1))

    means = gmm.means_.flatten()
    covariances = np.sqrt(gmm.covariances_).flatten()
    weights = gmm.weights_

    ground_gaussian = np.argmin(means)
    excited_gaussian = 1 - ground_gaussian

    # Generate x values to approximate the crossing point
    x_vals = np.linspace(means[ground_gaussian] - 3 * covariances[ground_gaussian],
                         means[excited_gaussian] + 3 * covariances[excited_gaussian], 1000)

    # Calculate Gaussian fits for each x value
    ground_gaussian_fit = weights[ground_gaussian] * (1 / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
        -0.5 * ((x_vals - means[ground_gaussian]) / covariances[ground_gaussian]) ** 2)
    excited_gaussian_fit = weights[excited_gaussian] * (
            1 / (np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
        -0.5 * ((x_vals - means[excited_gaussian]) / covariances[excited_gaussian]) ** 2)

    # Find the x value where the two Gaussian functions are closest
    crossing_point = x_vals[np.argmin(np.abs(ground_gaussian_fit - excited_gaussian_fit))]

    labels = gmm.predict(iq_data.reshape(-1, 1))

    ground_data = iq_data[(labels == ground_gaussian) & (iq_data < crossing_point)]
    excited_data = iq_data[(labels == excited_gaussian) & (iq_data > crossing_point)]

    ground_state_population = len(ground_data) / len(iq_data)
    excited_state_population_overlap = len(excited_data) / len(iq_data)

    return ground_state_population, excited_state_population_overlap, gmm, means, covariances, weights, crossing_point, ground_gaussian, excited_gaussian, ground_data, excited_data, iq_data


# Load data
base_dir = '/data/QICK_data/6transmon_run4/Non_RR_data'

# Create a "plots" folder if it doesnt already exist
plots_folder = os.path.join(base_dir, "plots")
os.makedirs(plots_folder, exist_ok=True)


# Lists to hold temperatures from Qubit_3 and Qubit_4
qubit_3_temperatures = []
qubit_4_temperatures = []

# Traverse date folders (e.g. '2024-11-13', '2024-11-14', etc.)
for date_folder in os.listdir(base_dir):
    date_folder_path = os.path.join(base_dir, date_folder)

    # Skip if it's not actually a directory
    if not os.path.isdir(date_folder_path):
        continue

    # Now specifically look for the "Qubit_Temps" subfolder in each date folder
    qubit_temps_folder = os.path.join(date_folder_path, "Qubit_Temps")
    if not os.path.isdir(qubit_temps_folder):
        print(f"Folder 'Qubit_Temps' not found under: {date_folder_path}. Skipping.")
        continue

    # Loop over .h5 files inside Qubit_Temps
    for filename in os.listdir(qubit_temps_folder):
        if not filename.endswith('.h5'):
            continue

        file_path = os.path.join(qubit_temps_folder, filename)

        # Check if it's Qubit 3 or Qubit 4

        #Qubit 3
        if filename.startswith('Qubit_3'):
            with h5py.File(file_path, 'r') as f:
                ig_new = f['ig_new'][:]
                qubit_frequency = f['qubit_frequency'][()]

            ground_state_population, excited_state_population_overlap, gmm, means, covariances, weights, crossing_point, ground_gaussian, excited_gaussian, ground_data, excited_data, iq_data= fit_double_gaussian_with_full_coverage(ig_new)

            # Calculate temperature
            temperature_k = calculate_qubit_temperature(qubit_frequency,
                                                        ground_state_population,
                                                        excited_state_population_overlap)
            # Check if the temperature_k is valid (not None) and below 500 mK (0.5 K)
            if temperature_k is not None and temperature_k < 0.5:
                temperature_mk = temperature_k * 1e3
                qubit_3_temperatures.append(temperature_mk)
                print(f"[{filename}]  Qubit 3: {temperature_mk:.2f} mK (Saved)")

                #-----------------------------------------PLOTS FOR Q3 FITS-----------------------------------------
                q_key = 2
                # Plotting double gaussian distributions and fitting
                xlims = [np.min(ig_new), np.max(ig_new)]
                plt.figure(figsize=(10, 6))

                # Plot histogram for `ig_new`
                steps = 3000
                numbins = round(math.sqrt(steps))
                n, bins, _ = plt.hist(ig_new, bins=numbins, range=xlims, density=False, alpha=0.5,
                                      label='Histogram of $I_g$',
                                      color='gray')
                # print(numbins)
                # Use the midpoints of bins to create boolean masks
                bin_centers = (bins[:-1] + bins[1:]) / 2
                ground_region = (bin_centers < crossing_point)
                excited_region = (bin_centers >= crossing_point)

                # Calculate scaling factors for each region
                scaling_factor_ground = max(n[ground_region]) / max(
                    (weights[ground_gaussian] / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
                        -0.5 * ((bin_centers[ground_region] - means[ground_gaussian]) / covariances[
                            ground_gaussian]) ** 2))

                scaling_factor_excited = max(n[excited_region]) / max(
                    (weights[excited_gaussian] / (np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
                        -0.5 * ((bin_centers[excited_region] - means[excited_gaussian]) / covariances[
                            excited_gaussian]) ** 2))

                # Generate x values for plotting Gaussian components
                x = np.linspace(xlims[0], xlims[1], 1000)
                ground_gaussian_fit = scaling_factor_ground * (
                        weights[ground_gaussian] / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
                    -0.5 * ((x - means[ground_gaussian]) / covariances[ground_gaussian]) ** 2)
                excited_gaussian_fit = scaling_factor_excited * (
                        weights[excited_gaussian] / (np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
                    -0.5 * ((x - means[excited_gaussian]) / covariances[excited_gaussian]) ** 2)

                plt.plot(x, ground_gaussian_fit, label='Ground Gaussian Fit', color='blue', linewidth=2)
                plt.plot(x, excited_gaussian_fit, label='Excited (leakage) Gaussian Fit', color='red', linewidth=2)
                plt.axvline(crossing_point, color='black', linestyle='--', linewidth=1,
                            label=f'Crossing Point ({crossing_point:.2f})')

                # Add shading for ground and excited state regions
                x_vals = np.linspace(np.min(ig_new), np.max(ig_new), 1000)

                # Add shading for ground_data points
                plt.hist(
                    ground_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                    alpha=0.5, color="blue", label="Ground Data Region", zorder=2
                )

                # Add shading for excited_data points
                plt.hist(
                    excited_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                    alpha=0.5, color="red", label="Excited Data Region", zorder=3
                )

                # plt.hist(
                #     iq_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                #     alpha=0.2, color="green", label="All IQ Data Region", zorder=1
                # )

                plt.title(
                    f"Fidelity Histogram and Double Gaussian Fit ; Qubit {q_key + 1}; Temperature = {temperature_mk:2f} mK")
                plt.xlabel('$I_g$', fontsize=14)
                #plt.ylabel('Probability Density', fontsize=14) or is it counts?
                plt.legend()
                #plt.show()

                # Save the plot to the Temperatures folder
                plot_filename = os.path.join(plots_folder,
                                             f"Qubit{q_key + 1}_Fidelityhist_gaussianfit_{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}.png")
                plt.savefig(plot_filename)
                print(f"Plot saved to {plot_filename}")
                plt.close()

                # Cleanup
                del ig_new, ground_gaussian_fit, excited_gaussian_fit, n, bins, x, gmm, \
                    bin_centers, ground_region, excited_region, scaling_factor_ground, \
                    scaling_factor_excited, weights, means, covariances, ground_gaussian, \
                    excited_gaussian, crossing_point

            else:
                print(f"[{filename}]  Qubit 3: Temperature >= 500 mK or invalid. Skipping.")

        #Qubit 4
        elif filename.startswith('Qubit_4'):
            with h5py.File(file_path, 'r') as f:
                ig_new = f['ig_new'][:]
                qubit_frequency = f['qubit_frequency'][()]

            ground_state_population, excited_state_population_overlap, gmm, means, covariances, weights, crossing_point, ground_gaussian, excited_gaussian, ground_data, excited_data, iq_data = fit_double_gaussian_with_full_coverage(ig_new)

            # Calculate temperature
            temperature_k = calculate_qubit_temperature(qubit_frequency,
                                                        ground_state_population,
                                                        excited_state_population_overlap)
            # Check if the temperature_k is valid (not None) and below 500 mK (0.5 K)
            if temperature_k is not None and temperature_k < 0.5:
                temperature_mk = temperature_k * 1e3
                qubit_4_temperatures.append(temperature_mk)
                print(f"[{filename}]  Qubit 4: {temperature_mk:.2f} mK (Saved)")

                # -----------------------------------------PLOTS FOR Q4 FITS-----------------------------------------
                q_key = 3
                # Plotting double gaussian distributions and fitting
                xlims = [np.min(ig_new), np.max(ig_new)]
                plt.figure(figsize=(10, 6))

                # Plot histogram for `ig_new`
                steps = 3000
                numbins = round(math.sqrt(steps))
                n, bins, _ = plt.hist(ig_new, bins=numbins, range=xlims, density=False, alpha=0.5,
                                      label='Histogram of $I_g$',
                                      color='gray')
                # print(numbins)
                # Use the midpoints of bins to create boolean masks
                bin_centers = (bins[:-1] + bins[1:]) / 2
                ground_region = (bin_centers < crossing_point)
                excited_region = (bin_centers >= crossing_point)

                # Calculate scaling factors for each region
                scaling_factor_ground = max(n[ground_region]) / max(
                    (weights[ground_gaussian] / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
                        -0.5 * ((bin_centers[ground_region] - means[ground_gaussian]) / covariances[
                            ground_gaussian]) ** 2))

                scaling_factor_excited = max(n[excited_region]) / max(
                    (weights[excited_gaussian] / (np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
                        -0.5 * ((bin_centers[excited_region] - means[excited_gaussian]) / covariances[
                            excited_gaussian]) ** 2))

                # Generate x values for plotting Gaussian components
                x = np.linspace(xlims[0], xlims[1], 1000)
                ground_gaussian_fit = scaling_factor_ground * (
                        weights[ground_gaussian] / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
                    -0.5 * ((x - means[ground_gaussian]) / covariances[ground_gaussian]) ** 2)
                excited_gaussian_fit = scaling_factor_excited * (
                        weights[excited_gaussian] / (np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
                    -0.5 * ((x - means[excited_gaussian]) / covariances[excited_gaussian]) ** 2)

                plt.plot(x, ground_gaussian_fit, label='Ground Gaussian Fit', color='blue', linewidth=2)
                plt.plot(x, excited_gaussian_fit, label='Excited (leakage) Gaussian Fit', color='red', linewidth=2)
                plt.axvline(crossing_point, color='black', linestyle='--', linewidth=1,
                            label=f'Crossing Point ({crossing_point:.2f})')

                # Add shading for ground and excited state regions
                x_vals = np.linspace(np.min(ig_new), np.max(ig_new), 1000)

                # Add shading for ground_data points
                plt.hist(
                    ground_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                    alpha=0.5, color="blue", label="Ground Data Region", zorder=2
                )

                # Add shading for excited_data points
                plt.hist(
                    excited_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                    alpha=0.5, color="red", label="Excited Data Region", zorder=3
                )

                # plt.hist(
                #     iq_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                #     alpha=0.2, color="green", label="All IQ Data Region", zorder=1
                # )

                plt.title(
                    f"Fidelity Histogram and Double Gaussian Fit ; Qubit {q_key + 1}; Temperature = {temperature_mk:2f} mK")
                plt.xlabel('$I_g$', fontsize=14)
                plt.ylabel('Probability Density', fontsize=14)
                plt.legend()
                #plt.show()

                plot_filename = os.path.join(plots_folder,
                                             f"Qubit{q_key + 1}_Fidelityhist_gaussianfit_{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}.png")
                plt.savefig(plot_filename)
                print(f"Plot saved to {plot_filename}")
                plt.close()

                # Cleanup
                del ig_new, ground_gaussian_fit, excited_gaussian_fit, n, bins, x, gmm, \
                    bin_centers, ground_region, excited_region, scaling_factor_ground, \
                    scaling_factor_excited, weights, means, covariances, ground_gaussian, \
                    excited_gaussian, crossing_point

            else:
                print(f"[{filename}]  Qubit 4: Temperature >= 500 mK or invalid. Skipping.")

#scatter plots
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.scatter(range(len(qubit_3_temperatures)), qubit_3_temperatures, alpha=0.7, color='blue')
plt.title("Qubit 3 Temperatures")
plt.xlabel("Measurement Index")
plt.ylabel("Temperature (mK)")

plt.subplot(1, 2, 2)
plt.scatter(range(len(qubit_4_temperatures)), qubit_4_temperatures, alpha=0.7, color='red')
plt.title("Qubit 4 Temperatures")
plt.xlabel("Measurement Index")
plt.ylabel("Temperature (mK)")

plt.tight_layout()

# Save scatter plot to the "plots" folder
scatter_plot_path = os.path.join(plots_folder, "qubit_temperatures_scatter.png")
plt.savefig(scatter_plot_path, dpi=150)  # Increase dpi if you want higher resolution

#plt.show()
plt.close()

#Histograms
plt.figure(figsize=(10, 5))

# Subplot 1: Qubit 3 histogram
plt.subplot(1, 2, 1)
plt.hist(qubit_3_temperatures, bins=30, alpha=0.7, color='blue')
plt.title("Qubit 3 Temperature Distribution")
plt.xlabel("Temperature (mK)")
plt.ylabel("Count")

# Subplot 2: Qubit 4 histogram
plt.subplot(1, 2, 2)
plt.hist(qubit_4_temperatures, bins=30, alpha=0.7, color='red')
plt.title("Qubit 4 Temperature Distribution")
plt.xlabel("Temperature (mK)")
plt.ylabel("Count")

plt.tight_layout()

# Save scatter plot to the "plots" folder
scatter_plot_path = os.path.join(plots_folder, "qubit_temperatures_hist.png")
plt.savefig(scatter_plot_path, dpi=150)  # Increase dpi if you want higher resolution

#plt.show()
plt.close()


