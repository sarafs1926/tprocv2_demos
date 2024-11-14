import numpy as np
import h5py
from sklearn.mixture import GaussianMixture

def calculate_qubit_temperature(frequency_mhz, ground_state_population, excited_state_population):
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck's constant in J·s
    frequency_hz = frequency_mhz * 1e6
    T = (h * frequency_hz) / (k_B * np.log(ground_state_population / excited_state_population))
    return T

def fit_double_gaussian_with_overlap_exclusion(iq_data):
    # Initialize and fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2)
    gmm.fit(iq_data.reshape(-1, 1))

    # Get the means and covariances to find the "true" ground and excited (leakage) states
    means = gmm.means_.flatten()
    covariances = np.sqrt(gmm.covariances_).flatten()

    # Identify the ground state Gaussian as the one with the lower mean
    ground_gaussian = np.argmin(means)
    excited_gaussian = 1 - ground_gaussian  # Other Gaussian is considered the excited state overlap. Since there are only two Gaussians, the excited Gaussian will have the other index (either 0 or 1)

    # Calculate the threshold (midpoint) between the two Gaussian means
    midpoint = (means[ground_gaussian] + means[excited_gaussian]) / 2

    # Label each data point according to which Gaussian it belongs to
    labels = gmm.predict(iq_data.reshape(-1, 1))

    # Exclude data points within the overlapping region defined by the midpoint
    ground_data = iq_data[(labels == ground_gaussian) & (iq_data < midpoint)]
    excited_data = iq_data[(labels == excited_gaussian) & (iq_data > midpoint)]

    # Calculate populations as the fraction of data points in each Gaussian, excluding the shared region
    ground_state_population = len(ground_data) / len(iq_data)
    excited_state_population_overlap = len(excited_data) / len(iq_data)

    #print(f"Ground state population (true): {ground_state_population}")
    #print(f"Excited state population (overlap): {excited_state_population_overlap}")

    return ground_state_population, excited_state_population_overlap


#Select file you want to look at
file_path = '/data/QICK_data/6transmon_run4a/2024-11-13/Qubit_Temps/Qubit_4_temperatureSSdata_20241113_132158.h5'

# Open the HDF5 file and extract necessary data
with h5py.File(file_path, 'r') as f:
    ig_new = f['ig_new'][:]  # Load rotated ground-state I data
    qubit_frequency = f['qubit_frequency'][()]  # Load qubit frequency
    saved_temperature_mk = f['temperature_millikelvin'][()]  # Load saved temperature to compare with recalculated one

# Recalculate ground and excited state populations
ground_state_population, excited_state_population_overlap = fit_double_gaussian_with_overlap_exclusion(ig_new)

# Recalculate and print the temperature
temperature_k = calculate_qubit_temperature(qubit_frequency, ground_state_population, excited_state_population_overlap)
temperature_mk = temperature_k * 1e3  # Convert to millikelvin

print(f"Ground state population: {ground_state_population}")
print(f"Excited state (leakage) population: {excited_state_population_overlap}")
print(f"Recalculated Temperature: {temperature_mk:.2f} mK")
print(f"Saved Temperature in file: {saved_temperature_mk:.2f} mK")
