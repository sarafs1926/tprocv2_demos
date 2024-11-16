import numpy as np

def convert_to_float_list(data):
    if isinstance(data, np.ndarray) and data.dtype == '|S8':
        return [float(x.decode('utf-8')) for x in data] #Corrected to handle the np.ndarray directly.
    elif isinstance(data, (list, tuple)):
        return [convert_to_float_list(item) for item in data]
    elif isinstance(data, dict):
        return {key: convert_to_float_list(value) for key, value in data.items()}
    else:
        return data

# Example input
data = [[np.bytes_(b'-3.5'), np.bytes_(b'-3.38'), np.bytes_(b'-3.26')]]  # Add your complete list here
result = convert_to_float_list(data)

print(result)