# T_robinson_valid_drawing

script to compute a valid drawing for a Robinson matrix

## Requirements
- pandas
- matplotlib
- pulp

## Usage

<code>python caterpillar_metric.py  [OPTIONS] [FILE]<\code>

### Example
  - Use the third predefined exmaple and plot it 
    python caterpillar_metric.py -e 3 -v
  - 

### Parameters 
-e, --example_number=INT
                   use the predefined example INT

-f, --input_file=FILE
                   use the FILE (.csv) as input

-g, --generate_random=INT
                   generate a random matrix with size INT

-v
                   verbose: print informationn and plots

-d
                   use a distance matrix instead of a dissimilarity

-p
                   forse solution to be a path instead of centipede
