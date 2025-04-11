import csv
import json

def parse_block(block):
    """
    Given a block of rows (each row is a list of strings from the CSV)
    representing one enzyme, parse the enzyme name, substrate, product,
    kinetic law equation, and the parameters into a dictionary.
    """
    enzyme_data = {}

    # 1. The enzyme name is the first cell of the first row.
    enzyme_data["enzyme"] = block[0][0].strip()

    # Initialize fields in case they are missing.
    enzyme_data["substrate"] = ""
    enzyme_data["product"] = ""
    enzyme_data["equation"] = ""
    enzyme_data["parameters"] = {}

    # iterate through each row in the block
    for i, row in enumerate(block):
        # Use lower() to allow for variations in casing.
        first_cell = row[0].strip().lower() if row[0] else ""

        if first_cell == "substrate" and len(row) > 1:
            enzyme_data["substrate"] = row[1].strip()
        elif first_cell == "product" and len(row) > 1:
            enzyme_data["product"] = row[1].strip()
        elif row[0].strip() == "Kinetic Law":
            # We assume the next two rows exist:
            #   next row is the header row (with 'type','formula',...)
            #   following row holds the kinetic law details.
            if i + 2 < len(block):
                kinetic_row = block[i + 2]
                if len(kinetic_row) > 1:
                    enzyme_data["equation"] = kinetic_row[1].strip()
        elif row[0].strip() == "Parameter":
            # The next row is the parameter header.
            # All subsequent non-empty rows are the parameter rows.
            param_index = i + 2  # parameters start two rows after "Parameter"
            while param_index < len(block):
                param_row = block[param_index]
                # if the row is completely empty, we can stop processing parameters.
                if not any(cell.strip() for cell in param_row):
                    break
                # Extract parameter info.
                # We expect: name, type, species, start val., unit
                if len(param_row) >= 4:
                    name = param_row[0].strip()
                    param_type = param_row[1].strip()
                    species = param_row[2].strip()
                    start_val = param_row[3].strip()
                    # Only add parameters with a non-empty name.
                    if name:
                        enzyme_data["parameters"][name] = {
                            "type": param_type,
                            "species": species,
                            "start_val": start_val
                        }
                param_index += 1
            break

    return enzyme_data

def parse_csv(filename):
    """
    Reads the CSV file and splits it into enzyme blocks. Two consecutive
    completely empty rows are used as the delimiter between enzymes.
    Returns a list of enzyme dictionaries.
    """
    enzymes = []
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        block = []
        empty_line_count = 0

        for row in reader:
            if not any(cell.strip() for cell in row):
                empty_line_count += 1
                if empty_line_count >= 2 and block:
                    enzyme = parse_block(block)
                    enzymes.append(enzyme)
                    block = []
                    empty_line_count = 0
            else:
                empty_line_count = 0
                block.append(row)

        if block:
            enzyme = parse_block(block)
            enzymes.append(enzyme)

    return enzymes

if __name__ == '__main__':
    # Replace 'input.csv' with your CSV file path.
    filename = "glycolysis.csv"
    enzyme_list = parse_csv(filename)
    
    # Print the output in JSON format for clarity.
    with open("glycolysis.json", "w") as f:
        json.dump(enzyme_list, f, indent=4)
    
    #print(json.dumps(enzyme_list, indent=4))
