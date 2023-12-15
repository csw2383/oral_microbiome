import matplotlib.pyplot as plt
import os


def extract(file_path):
    genus_data = {}
    sum = 0.0

    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        parts = line.strip().split('\t')

        if len(parts) > 1:
            taxa_info = parts[0]

            taxa_elements = taxa_info.split('|')
            genuses = [elem for elem in taxa_elements if elem.startswith('g__')]

            if genuses:
                genus_name = genuses[0].split('__')[1]
                if 'incert' in genus_name:
                    genus_name = 'Unknown'  

                val = float(parts[1])

                if genus_name in genus_data:
                    genus_data[genus_name] += val
                else:
                    genus_data[genus_name] = val

                sum += val

    genus_data = {genus: val / sum for genus, val in genus_data.items()}

    return genus_data

def pie(genus_data,output_file):
    sorted_genus_data = sorted(genus_data.items(), key=lambda x: x[1]) 
    labels = [item[0] for item in sorted_genus_data]
    sizes = [item[1] for item in sorted_genus_data]

    plt.figure(figsize=(10, 10))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
    plt.axis('equal')
    plt.savefig(output_file)


def main():
    directory = '/Users/weiwei/cu/cs4775/motus_output'  # Replace this with your directory path

    for file in os.listdir(directory):
        if file.endswith('.motus.ratio.genus'):  
            file_path = os.path.join(directory, file)
            extracted_data = extract(file_path)
            output_file = os.path.splitext(file)[0] + '.png'  
            pie(extracted_data, output_file)

if __name__ == "__main__":
    main()