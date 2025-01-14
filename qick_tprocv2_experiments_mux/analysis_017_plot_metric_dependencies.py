import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime

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
            ax.set_title(titles[i], fontsize=font)

            datetime_objects_1 = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S")
                                  for date_string in date_times_1[i]]
            datetime_objects_2 = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S")
                                  for date_string in date_times_2[i]]

            combined_1 = list(zip(datetime_objects_1, metric_1[i]))  # frequency
            combined_2 = list(zip(datetime_objects_2, metric_2[i]))  # T1

            combined_1.sort(reverse=True, key=lambda item: item[0])
            combined_2.sort(reverse=True, key=lambda item: item[0])

            sorted_date_times_1, sorted_metric_1 = zip(*combined_1)
            sorted_date_times_2, sorted_metric_2 = zip(*combined_2)

            # Convert to NumPy arrays
            sorted_date_times_1 = np.array(sorted_date_times_1)
            sorted_date_times_2 = np.array(sorted_date_times_2)
            sorted_metric_1 = np.array(sorted_metric_1)  # frequency
            sorted_metric_2 = np.array(sorted_metric_2)  # T1

            # ------------------------------------------------------
            # Remove the "pick the shorter array" logic. We assume:
            #   metric_1 == qubit frequency
            #   metric_2 == T1
            # ------------------------------------------------------
            ref_times = sorted_date_times_1
            ref_metrics = sorted_metric_1  # frequency
            other_times = sorted_date_times_2
            other_mets = sorted_metric_2  # T1

            # Now match times so that for each timestamp in ref_times,
            # we find the closest timestamp in other_times.
            matched_ref_metrics = []
            matched_other_metrics = []

            for t_ref, m_ref in zip(ref_times, ref_metrics):
                # Index of closest timestamp in other_times
                idx_closest = np.argmin(np.abs(other_times - t_ref))
                matched_ref_metrics.append(m_ref)
                matched_other_metrics.append(other_mets[idx_closest])

            matched_ref_metrics = np.array(matched_ref_metrics)
            matched_other_metrics = np.array(matched_other_metrics)

            # ------------------------------------------------------
            # Always plot qubit frequency on the x-axis (metric_1)
            # and T1 on the y-axis (metric_2).
            # ------------------------------------------------------
            ax.scatter(matched_ref_metrics, matched_other_metrics, color='blue')
            ax.set_xlabel("Q Freq", fontsize=font - 2)
            ax.set_ylabel("T1 (us)", fontsize=font - 2)


        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(analysis_folder + f'{metric_1_label}_vs_{metric_2_label}_correlation.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

    # -------------------------------- Single scatter plot for two metrics ---------------------------------------------
    def plot_single_pair(self, date_times_1, metric_1, date_times_2, metric_2, metric_1_label="Qubit 1 T1", metric_2_label="Qubit 3 T1"):
        """
        Creates a SINGLE scatter plot of metric_1 vs metric_2,
        matching data points by the same 'closest timestamp' logic.
        """
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)

        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/correlations_singleplots/"
        self.create_folder_if_not_exists(analysis_folder)

        #Converts timestamps from strings to datetime objects
        datetime_objects_1 = [datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")
                              for dt_str in date_times_1]
        datetime_objects_2 = [datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")
                              for dt_str in date_times_2]

        #Zip up times+metrics and sort descending by time
        combined_1 = list(zip(datetime_objects_1, metric_1))
        combined_2 = list(zip(datetime_objects_2, metric_2))

        combined_1.sort(reverse=True, key=lambda item: item[0])
        combined_2.sort(reverse=True, key=lambda item: item[0])

        #Unzip
        sorted_date_times_1, sorted_metric_1 = zip(*combined_1) if combined_1 else ([], [])
        sorted_date_times_2, sorted_metric_2 = zip(*combined_2) if combined_2 else ([], [])

        #Convert to np.array for easy math
        sorted_date_times_1 = np.array(sorted_date_times_1)
        sorted_metric_1 = np.array(sorted_metric_1)
        sorted_date_times_2 = np.array(sorted_date_times_2)
        sorted_metric_2 = np.array(sorted_metric_2)

        #Pick the shorter array as reference
        if len(sorted_date_times_1) <= len(sorted_date_times_2):
            ref_times = sorted_date_times_1
            ref_metrics = sorted_metric_1
            other_times = sorted_date_times_2
            other_metrics = sorted_metric_2
            x_label = metric_1_label
            y_label = metric_2_label
        else:
            ref_times = sorted_date_times_2
            ref_metrics = sorted_metric_2
            other_times = sorted_date_times_1
            other_metrics = sorted_metric_1
            x_label = metric_2_label
            y_label = metric_1_label

        #Match each ref timestamp with the closest timestamp in other_times
        matched_ref = []
        matched_other = []

        for t_ref, m_ref in zip(ref_times, ref_metrics):
            idx_closest = np.argmin(np.abs(other_times - t_ref))
            matched_ref.append(m_ref)
            matched_other.append(other_metrics[idx_closest])

        matched_ref = np.array(matched_ref)
        matched_other = np.array(matched_other)

        # Creates a SINGLE scatter plot
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(matched_ref, matched_other, color='blue')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(f"{metric_1_label} vs {metric_2_label}")

        plt.tight_layout()
        #plt.show()
        plt.savefig(analysis_folder + f'{metric_1_label}_vs_{metric_2_label}_correlation.png', transparent=False,
                    dpi=self.final_figure_quality)


    #----------------------------------Plots Qubit 1 Temperature and T1 vs. Time --------------------------------------------
    def plot_q1_temp_and_t1(
            self,
            q1_temp_times,
            q1_temps,
            q1_t1_times,
            q1_t1_vals,
            temp_label="mK",
            t1_label="T1 (µs)",
            magcan_dates=None,
            magcan_temps=None,
            magcan_label="mK",
            mcp2_dates=None,
            mcp2_temps=None,
            mcp2_label="mK",
            Q1_freqs=None,
            Q1_dates_spec=None,
            qspec_label="Q1 Frequency (MHz)",
            date_times_pi_amps_Q1 = None,
            pi_amps_Q1 = None,
            pi_amps_label = "Pi Amp (a.u.)"):
        """
        Plots Qubit 1's temperature vs. time and T1 vs. time on the same figure,
        using two y-axes (a left y-axis for T1, and a right y-axis for temperature).
        Also plots thermometry data, Pi Amp, and qubit frequency data during this time frame if provided.
        """

        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)

        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/correlations_singleplots/"
        self.create_folder_if_not_exists(analysis_folder)

        # If timestamps are strings, convert them to datetime objects
        def ensure_datetime(ts_list):
            if not ts_list:
                return []
            if isinstance(ts_list[0], str):
                return [datetime.strptime(t_str, "%Y-%m-%d %H:%M:%S") for t_str in ts_list]
            return ts_list

        q1_temp_times_dt = ensure_datetime(q1_temp_times)
        q1_t1_times_dt = ensure_datetime(q1_t1_times)

        if mcp2_dates is not None:
            mcp2_dates_dt = ensure_datetime(mcp2_dates)
        else:
            mcp2_dates_dt = []

        if magcan_dates is not None:
            magcan_dates_dt = ensure_datetime(magcan_dates)
        else:
            magcan_dates_dt = []

        if Q1_dates_spec is not None:
            Q1_dates_spec_dt = ensure_datetime(Q1_dates_spec)
        else:
            Q1_dates_spec_dt = []

        if date_times_pi_amps_Q1 is not None:
            date_times_pi_amps_Q1_dt = ensure_datetime(date_times_pi_amps_Q1)
        else:
            date_times_pi_amps_Q1_dt = []

        # figure
        fig, ax1 = plt.subplots(figsize=(14, 6))
        ax2 = ax1.twinx()  # second y-axis (right side) for temperature

        #Plots T1 data on ax1 (left y-axis)
        ax1.scatter(q1_t1_times_dt, q1_t1_vals, color='blue', alpha=0.8, label=t1_label, edgecolor='black')
        ax1.set_xlabel("Time")
        ax1.set_ylabel(t1_label, color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')
        ax1.tick_params(axis='x', rotation=45)
        ax1.grid(True, alpha=0.3)

        ############################# This sets a limit on the x axis if you only want to look at specific dates
        year = 2024
        month = 12
        day = 20
        start_time = datetime(year, month, day, 8, 0, 0)
        end_time = datetime(year, month, day, 18, 0, 0)
        ax1.set_xlim(start_time, end_time)  # Set the x-axis limits
        #############################

        #Plots temperature data on ax2 (right y-axis)
        ax2.scatter(q1_temp_times_dt, q1_temps, color='red', alpha=0.8, label=temp_label, edgecolor='black')
        ax2.set_ylabel("Temperature (mK)", color='black')
        ax2.tick_params(axis='y', labelcolor='black')

        # plot thermometry data if available
        if mcp2_dates_dt and mcp2_temps:
            ax2.scatter(
                mcp2_dates_dt, mcp2_temps,
                color='orange', alpha=0.8, label=mcp2_label
            )

        if magcan_dates_dt and magcan_temps:
            ax2.scatter(
                magcan_dates_dt, magcan_temps,
                color='green', alpha=0.8, label=magcan_label
            )

        # Plots frequency data on ax3
        if Q1_dates_spec_dt and Q1_freqs:
            ax3 = ax1.twinx()  # third y-axis (right side) for frequency
            ax3.spines["right"].set_position(("outward", 60))  # Offset the third axis
            ax3.scatter(Q1_dates_spec_dt, Q1_freqs, color='orchid', alpha=0.8, label=qspec_label, edgecolor='black')
            ax3.set_ylabel("Frequency (MHz)", color='purple')
            ax3.tick_params(axis='y', labelcolor='purple')

        # Plots frequency data on ax4
        if date_times_pi_amps_Q1_dt and pi_amps_Q1:
            ax4 = ax1.twinx()  # fourth y-axis (right side) for Pi Amp (a.u.)
            ax4.spines["right"].set_position(("outward", 120))  # Offset the third axis
            ax4.scatter(date_times_pi_amps_Q1_dt, pi_amps_Q1, color='gold', alpha=0.8, label=pi_amps_label, edgecolor='black')
            ax4.set_ylabel("Pi Amp (a.u.)", color='goldenrod')
            ax4.tick_params(axis='y', labelcolor='goldenrod')

        plt.title("Qubit 1: T1 , Freq, Pi Amp, and Temperature vs. Time")

        # Collect legend info from ax1 and ax2
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        handles3, labels3 = ax3.get_legend_handles_labels()
        handles4, labels4= ax4.get_legend_handles_labels()

        # Make one combined legend on ax1 (or plt)
        ax1.legend(handles1 + handles2 + handles3 + handles4, labels1 + labels2 + labels3 + labels4, loc='upper left', fontsize=10, markerscale=0.8, handletextpad=0.6, borderpad=0.7, labelspacing=0.4)

        fig.tight_layout()

        unique_str = datetime.now().strftime('%Y%m%d%H%M%S')
        plot_filename = os.path.join(analysis_folder, f"Q1_Temp_and_T1_vs_Time_{unique_str}.png")
        plt.savefig(plot_filename, transparent=False, dpi=self.final_figure_quality)
        #plt.show()
        #plt.close()