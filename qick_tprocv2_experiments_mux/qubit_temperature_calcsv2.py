import numpy as np
import h5py
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import math

#this script is outdated,check out analysis_014 script instead

def calculate_qubit_temperature(frequency_mhz, ground_state_population, excited_state_population):
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck's constant in J·s
    frequency_hz = frequency_mhz * 1e6
    T = (h * frequency_hz) / (k_B * np.log(ground_state_population / excited_state_population))
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

    return ground_state_population, excited_state_population_overlap, gmm, means, covariances, weights, crossing_point, ground_gaussian, excited_gaussian


# Load data
file_path = '/data/QICK_data/6transmon_run4/2024-11-13/Qubit_Temps/Qubit_4_temperatureSSdata_20241113_132158.h5'
with h5py.File(file_path, 'r') as f:
    ig_new = f['ig_new'][:]
    qubit_frequency = f['qubit_frequency'][()]
    saved_temperature_mk = f['temperature_millikelvin'][()]
    saved_fidelity = f['fidelity'][()]

ground_state_population, excited_state_population_overlap, gmm, means, covariances, weights, crossing_point, ground_gaussian, excited_gaussian = fit_double_gaussian_with_full_coverage(
    ig_new)
temperature_k = calculate_qubit_temperature(qubit_frequency, ground_state_population, excited_state_population_overlap)
temperature_mk = temperature_k * 1e3

print(f"Ground state population: {ground_state_population}")
print(f"Excited state (leakage) population: {excited_state_population_overlap}")
print(f"Recalculated Temperature: {temperature_mk:.2f} mK")
print(f"Saved Temperature in file: {saved_temperature_mk:.2f} mK")

# Plotting
xlims = [np.min(ig_new), np.max(ig_new)]
plt.figure(figsize=(10, 6))

# Plot histogram for `ig_new`
steps = 3000
numbins = round(math.sqrt(steps))
n, bins, _ = plt.hist(ig_new, bins=numbins, range=xlims, density=False, alpha=0.5, label='Histogram of $I_g$',
                      color='gray')
#print(numbins)
# Use the midpoints of bins to create boolean masks
bin_centers = (bins[:-1] + bins[1:]) / 2
ground_region = (bin_centers < crossing_point)
excited_region = (bin_centers >= crossing_point)

# Calculate scaling factors for each region
scaling_factor_ground = max(n[ground_region]) / max(
    (weights[ground_gaussian] / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
        -0.5 * ((bin_centers[ground_region] - means[ground_gaussian]) / covariances[ground_gaussian]) ** 2))

scaling_factor_excited = max(n[excited_region]) / max(
    (weights[excited_gaussian] / (np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
        -0.5 * ((bin_centers[excited_region] - means[excited_gaussian]) / covariances[excited_gaussian]) ** 2))

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
plt.axvline(crossing_point, color='black', linestyle='--', linewidth=1, label=f'Crossing Point ({crossing_point:.2f})')

plt.title(f"Fidelity Histogram and Double Gaussian Fit ; Fidelity = {saved_fidelity * 100:.2f}%")
plt.xlabel('$I_g$', fontsize=14)
#plt.ylabel('Probability Density', fontsize=14)
plt.legend()
plt.show()

# Cleanup
del ig_new, ground_gaussian_fit, excited_gaussian_fit, n, bins, x, gmm, \
    bin_centers, ground_region, excited_region, scaling_factor_ground, \
    scaling_factor_excited, weights, means, covariances, ground_gaussian, \
    excited_gaussian, crossing_point