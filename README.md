# null-model
Repository to hold work on a new null model for HMMER and Infernal

File Descriptions: [last updated January 24 2021]

count_occurences.c -- calculates pattern counts and prints the first few to the terminal. does not output any files.
caculate_probabilities.c -- uses the same code as count_occurences.c to count pattern occurences and uses the counts to calculate conditional probabilities for each residue. does not produce any output files. prints first few probabilities as well as numerator and denominator for the fractions.
calculate_probabilities_with_output.c -- does the same as calculate_probabilities.c and outputs the pattern_counts into a file to compare with python (outputs pattern counts, not probabilities)
