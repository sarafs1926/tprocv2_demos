from datetime import datetime
import matplotlib.pyplot as plt

x = ["2023-10-02 14:00:00", "2023-05-06 16:45:00", "2022-01-12 09:30:00", "2023-05-06 17:45:00"]
y = [10, 1, 5, 8]

#
# # Convert strings to datetime objects.
# datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]
#
# # Combine datetime objects and y values into a list of tuples and sort by datetime.
# combined = list(zip(datetime_objects, y))
# combined.sort(reverse=True, key=lambda x: x[0])
#
# # Unpack them back into separate lists, in order from latest to most recent.
# sorted_x, sorted_y = zip(*combined)
# ax.scatter(sorted_x, sorted_y, color=colors[i])
#
# sorted_x = np.asarray(sorted(x))
#
# num_points = 5
# indices = np.linspace(0, len(sorted_x) - 1, num_points, dtype=int)
#
# # Set new x-ticks using the datetime objects at the selected indices
# ax.set_xticks(sorted_x[indices])
# ax.set_xticklabels([dt for dt in sorted_x[indices]], rotation=45)
#
# ax.scatter(x, y, color=colors[i])
# if show_legends:
#     ax.legend(edgecolor='black')
# ax.set_xlabel('Time (Days)', fontsize=font)
# ax.set_ylabel('Resonator Center (MHz)', fontsize=font)
# ax.tick_params(axis='both', which='major', labelsize=font)
