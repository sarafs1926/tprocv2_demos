import matplotlib.pyplot as plt
import numpy as np
import os

class PlotMetricDependencies:
    def __init__(self,run_name, number_of_qubits, final_figure_quality):
        self.run_name = run_name
        self.number_of_qubits = number_of_qubits
        self.final_figure_quality = final_figure_quality

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        if not os.path.exists(folder):
            os.makedirs(folder)

    def plot(self, date_times_1, metric_1, date_times_2, metric_2, metric_1_label,metric_2_label):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/metric_interdependencies/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('T2 Values vs Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime
        for i, ax in enumerate(axes):
            ax.set_title(titles[i], fontsize = font)

            datetime_objects_1 = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in date_times_1[i]]
            datetime_objects_2 = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in
                                  date_times_2[i]]

            combined_1 = list(zip(datetime_objects_1, metric_1[i]))
            combined_2 = list(zip(datetime_objects_2, metric_2[i]))

            combined_1.sort(reverse=True, key=lambda item: item[0])
            combined_2.sort(reverse=True, key=lambda item: item[0])

            sorted_date_times_1, sorted_metric_1 = zip(*combined_1)
            sorted_date_times_2, sorted_metric_2 = zip(*combined_2)

            sorted_date_times_1 = np.array(sorted_date_times_1)
            sorted_date_times_2 = np.array(sorted_date_times_2)

            sorted_metric_1 = np.array(sorted_metric_1)
            sorted_metric_2 = np.array(sorted_metric_2)

            sorted_date_times_1 = np.array(sorted_date_times_1)
            sorted_date_times_2 = np.array(sorted_date_times_2)
            sorted_metric_1 = np.array(sorted_metric_1)
            sorted_metric_2 = np.array(sorted_metric_2)

            # Step 1: Pick the shorter array as reference
            if len(sorted_date_times_1) <= len(sorted_date_times_2):
                ref_times = sorted_date_times_1
                ref_metrics = sorted_metric_1
                other_times = sorted_date_times_2
                other_mets = sorted_metric_2
            else:
                ref_times = sorted_date_times_2
                ref_metrics = sorted_metric_2
                other_times = sorted_date_times_1
                other_mets = sorted_metric_1

            # Step 2 & 3: For each timestamp in the reference array, find the closest
            # timestamp in the other array, and collect matched values
            matched_ref_metrics = []
            matched_other_metrics = []

            for t_ref, m_ref in zip(ref_times, ref_metrics):
                # Use np.argmin() on the absolute difference to find closest time in `other_times`
                idx_closest = np.argmin(np.abs(other_times - t_ref))

                # Append the reference metric and the matched metric from the other array
                matched_ref_metrics.append(m_ref)
                matched_other_metrics.append(other_mets[idx_closest])

            # Convert to NumPy arrays (optional but useful for subsequent processing)
            matched_ref_metrics = np.array(matched_ref_metrics)
            matched_other_metrics = np.array(matched_other_metrics)

            # Step 4: Now that both metrics have matching lengths and times, you can plot
            ax.scatter(matched_ref_metrics, matched_other_metrics, color='blue')

            ax.set_xlabel(metric_1_label, fontsize=font-2)
            ax.set_ylabel(metric_2_label, fontsize=font-2)


        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(analysis_folder + f'{metric_1_label}_vs_{metric_2_label}_correlation.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()
