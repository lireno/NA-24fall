import os
import matplotlib.pyplot as plt

def plot_multiple_splines(input_dir, output_path):
    # Initialize the plot
    plt.figure(figsize=(10, 8))
    
    # Iterate over all .txt files in the input directory
    for file_name in sorted(os.listdir(input_dir)):  # Sorting ensures consistent order
        if file_name.endswith(".txt"):
            file_path = os.path.join(input_dir, file_name)
            x, y = [], []
            
            # Read data from file
            with open(file_path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) == 2:
                        x.append(float(parts[0]))
                        y.append(float(parts[1]))
            
            # Plot the data with a label for the file
            plt.plot(x, y, label=file_name.replace(".txt", ""))
    
    # Add title, labels, and legend
    plt.title("Spline Curves")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend(title="Curves")
    plt.grid(True)
    
    # Save the combined plot
    plt.savefig(output_path)
    plt.close()
    print(f"Saved combined plot to {output_path}")

if __name__ == "__main__":
    input_directory = "./plotting/data/"
    output_path = "./plotting/figures/spline_curves_combined.png"
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Plot and save the combined splines
    plot_multiple_splines(input_directory, output_path)
