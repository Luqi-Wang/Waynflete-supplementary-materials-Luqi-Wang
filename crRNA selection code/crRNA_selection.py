# Python program coded by Luqi Wang, Magdalen College School, for selecting potential crRNA target sequences to target the 
#SARS-CoV-2 RdRp gene. Skeleton code could be easily adapted into selecting potential crRNA target sequences for other genes

SPACER_SEQUENCE_LENGTH = 22 # Declares a constant which stores the length of the spacer sequence 
# A function that calculates the annealing temperature of the crRNA and the target seuqence and returns the value
def calculate_annealing_temperature(GC_percentage, SPACER_SEQUENCE_LENGTH): # A 
    Annealing_temperature = 81.5 + 0.4 * GC_percentage - (675/SPACER_SEQUENCE_LENGTH) - 8 # This formula is cited in the paper
    return Annealing_temperature

crRNA_array = [] # An array used to store all the crRNA spacer sequences that have been selected
GC_Crietria = False
Temp_Criteria = False      # The three criterias that have to be True (satisfied) for the sequence to be selected
PolyU_Criteria = False
GC_count = 0  
GC_percentage = 0
Annealing_temperature = 0
start_index = 0
end_index = 22  
Total_crRNA = 0 # A count of all the crRNA sequences, including the ones that are not selected

# Opening the SARS-CoV-2 RdRp gene sequence file derived from SnapGene. 
genome_file =  open("SARS-CoV-2-RdRp.txt", 'r+', encoding = 'utf-8') 
genome = genome_file.read().replace('\n', '') # Reads the file and turns the file into a very long string
# This is a loop which goes through the entire RdRp gene sequence to gain all the possible 22 nucleotide spacer sequences.
while end_index < len(genome) + 1: 
    GC_Crietria = False
    Temp_Criteria = False      # All the values are set to False or 0 at the start of the loop
    PolyU_Criteria = False
    GC_count = 0
    GC_percentage = 0
    Annealing_temperature = 0
    crRNA = list(genome[start_index:end_index])

    for nucleotide in crRNA:   # This loops through all the nucleotides of the crRNA spacer sequence counting the number of Gs and Cs
        if nucleotide == 'G' or nucleotide == 'C':
            GC_count = GC_count + 1
    GC_percentage = GC_count / SPACER_SEQUENCE_LENGTH * 100 # Calculates the GC percentage of the entire 22 nucleotide spacer sequence
    
    if GC_percentage < 60 and GC_percentage > 40: # Checks if the GC percentage lies within 40 and 60
        GC_Crietria = True # If it does lie in between 40 and 60 I set the GC_Criteria to True.

    counter = 0
    while counter < len(crRNA): # I loop through all the nucleotides again
        if counter == len(crRNA) - 2:
            break                    # This time I check for poly-U sequences (more than 3 consecutive Uracil nucelotides together)
        if crRNA[counter] == 'U' and crRNA[counter+1] == 'U' and crRNA[counter+2] == 'U': 
            PolyU_Criteria = False
            break
        else:
            PolyU_Criteria = True # If the sequence does not contain any, I set the PolyU_Criteria to True
        counter = counter + 1
                        # I call the function to calculate the annealing temperature here.
    Annealing_temperature = calculate_annealing_temperature(GC_percentage, SPACER_SEQUENCE_LENGTH) 
    if Annealing_temperature > 55:  # I check if the annealing temperature is over 55 degrees Celsius, if yes I set Temp_crietria to True
        Temp_Criteria = True 

    if Temp_Criteria == True and PolyU_Criteria == True and GC_Crietria == True: # I check if the sequence satisfied all crieterias
        crRNA_array.append(crRNA) # If yes, I add the sequence to the crRNA_array

    start_index = start_index + 1   # I increase the start index and the end index to get the next 22 nucleotide sequence from the 
    end_index = end_index + 1       # RdRp gene.
    Total_crRNA = Total_crRNA + 1 # Counts the total number of crRNA including not selected ones 
genome_file.close()
# This prints out all of the selected crRNAs in the crRNA array

crRNA_array_file = open("crRNA-array-for-SARS-CoV-2-RdRp.txt", 'w', encoding = 'utf-8')
DIVIDER = '>'
for crRNA in crRNA_array: 
    nucleotide_sequence = ''
    for nucleotide in crRNA:
        nucleotide_sequence = nucleotide_sequence + nucleotide
    print(nucleotide_sequence)
    crRNA_array_file.write(DIVIDER)        # Formats the sequences in the crRNA array into a file in the format FASTA for 
    crRNA_array_file.write('\n')           # Multi query BLAST.
    crRNA_array_file.write(nucleotide_sequence)    
    crRNA_array_file.write('\n')
crRNA_array_file.close()
print('Target gene: SARS-CoV-2 RdRp')
print('Number of nucleotides in SARS-CoV-2 RdRp gene:')
print(len(genome))
print('Total number of crRNA spacer sequences possible for targeting SARS-CoV-2 RdRp: ')
print(Total_crRNA)
print('Total number of selected crRNAs: ')
print(len(crRNA_array))


