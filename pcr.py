# import required libraries
import random
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('tkagg')

# param: a DNA strand
# Return: complement from given strand
def getComp(strand):
    comp = strand.upper()  # Ensure all caps to prevent errors
    comp = comp.replace("A", "X")
    comp = comp.replace("T", "A")
    comp = comp.replace("X", "T")
    comp = comp.replace("C", "X")
    comp = comp.replace("G", "C")
    comp = comp.replace("X", "G")

    return comp


# param: null
# return: a tuple of two tuples where each tuple contain primers sequence, GC, start, and end position.
def getPrimers():
    # Primer pair (Sequence, GC, start, end)
    forwardPrimer = ("TGTACTCATTCGTTTCGGAA", 0.5, 7, 28)
    reversePrimer = ("AAGGACTAGAAGACCAGATT", 0.5, 134, 115)
    primers = (forwardPrimer, reversePrimer)

    return primers


# param: a list of tuples of 2 strs, representing double stranded dna segments
# return: a list of single strand dna segments
def denaturation(dna_segments):

    # list to store single strands DNA
    single_strands_DNA = []

    # iterate through all DNA items in DNA segments
    for item in dna_segments:

        # append the single strands DNA to the list
        single_strands_DNA.append(item[0])
        single_strands_DNA.append(item[1])

    # return the list of single strands DNA
    return single_strands_DNA


# param: a list of single strand DNAs, primer pair, and fall of rate
# return: a list of tuples of 2 strands of DNA
def annealing_elongation(singleStrandDNAs, primers, fall_of_rate=50):

    # For the annealing and elongation, Primers are needed
    # Pass the primers to the single strands of DNA
    # Bind them and extend until the fall of rate or the end of DNA segment

    # Get the primer sequence from primers tuple

    # Get forward primer
    forward_primer = primers[0][0]

    # Get reverse primer
    reverse_primer = primers[1][0]

    # get complement of forward and reverse primers
    comp_forward_primer = getComp(forward_primer)
    comp_reverse_primer = getComp(reverse_primer)
    comp_reverse_primer = comp_reverse_primer[::-1]

    # set the length of primers
    primer_length = 20

    # list to store DNA result after annealing and elongation part
    resultDNA = []

    # iterate through every single strand to anneal and elongate the DNA strand
    for item in singleStrandDNAs:

        # set the rate
        # d + r
        # d = 200, r = random number between (-50, 50)
        rate = 200 + random.randint(-fall_of_rate, fall_of_rate)

        # copy the item iterating to the firstDNAStrand variable
        # first strand of the list is the first signle strand to be worked on
        firstDNAStrand = item

        # if the firstDNAStrand is empty, continue
        if firstDNAStrand == "":
            continue

        secondDNAStrand = ""

        # if the complementary reverse primer is found in firstDNAStrand
        # then the work is to be done is with the second
        if firstDNAStrand.find(comp_reverse_primer) != -1:

            # Complement the reverse primer to find the end index
            checkPrimer = comp_reverse_primer

            # Get 5 -> 3 strand
            secondDNAStrand = firstDNAStrand

            # get the end index of the DNA strand by checking with with the complement of reverse primer
            endIndex = secondDNAStrand.index(checkPrimer)

            # Get the strand we need the complement of
            # Get the complement
            # Use endIndex to do so.
            secondDNAStrand = secondDNAStrand[endIndex:endIndex + primer_length + rate]

            # Get complement of the secondDNAStrand
            secondDNAStrand = getComp(secondDNAStrand)

            # Reverse the secondDNAStrand
            secondDNAStrand = secondDNAStrand[::-1]

        # if the first condition is not true
        # which is if complementary of forward primer is found in firstDNAStrand
        elif firstDNAStrand.find(comp_forward_primer) != -1:

            # Complement the forward primer to find the end index
            checkPrimer = comp_forward_primer

            secondDNAStrand = firstDNAStrand

            # get the start index of the DNA strand by checking with with the complement of reverse primer
            startIndex = secondDNAStrand.index(checkPrimer)

            # Get the strand we need the complement of
            # Get the complement
            # Use startIndex to do so.
            secondDNAStrand = secondDNAStrand[startIndex: startIndex + primer_length + rate]

            # Get complement of the secondDNAStrand
            secondDNAStrand = getComp(secondDNAStrand)

            # Reverse the secondDNAStrand
            secondDNAStrand = secondDNAStrand[::-1]

        # Append the result DNA with the result obtained from annealing and elongation
        resultDNA.append((firstDNAStrand, secondDNAStrand))

    # return the final resultDNA obtained after annealing and elongation.
    return resultDNA


# param: gene to be copied (a tuple of 2 strs), fall of rate of DNA polymerase (int), and num_cycles to run PCR (int)
# return: a list of double stranded dna segments
def PCR(dna_segment_to_be_copied, fall_of_rate, num_cycles):

    # get the primers from getPrimers() function
    primers = getPrimers()

    # initialize the cycles to 0 before the start of the PCR process
    cycles = 0

    # Store the DNA segments to be copied in PCRproducts list
    PCRproducts = [dna_segment_to_be_copied]

    # Begin the PCR cycles
    # Run the cycles for the desired number of cycles number passed as argument
    while cycles < num_cycles:

        # Denature the list of double-stranded DNAs at first
        # And store the denatured single strands of DNAs in singleStrandDNAs list
        singleStrandDNAs = denaturation(PCRproducts)

        # Begin the anealing and elongation process
        # store the obtained result in PCRproducts
        PCRproducts = annealing_elongation(singleStrandDNAs, primers, fall_of_rate)

        # Print the number of cycles completed after completion of each cyle
        print(" Cycle " + str(cycles + 1) + " completed \n")

        cycles += 1

    # Return the PCR products obtained.
    return PCRproducts

# param: Replicated list of DNA strands
# return: void
def stats(replicatedDNAs):

    # Get the statistics of the PCR results
    # Calculate:
    # Total Strands found, Average GC, Max length, Min length, Average Length
    # Distribution of the lengths of DNA fragments in graph using matplotlib

    # list of lengths of DNA segments
    lengthsOfSegments = []

    # list of GC contents
    GCcontents = []

    # iterate through each pair of replicated DNAs and each strand in a pair
    for DNApair in replicatedDNAs:
        for DNAstrand in DNApair:

            # if the strand is not empty, calculate the statistics
            if DNAstrand != '':
                lengthsOfSegments.append(len(DNAstrand))

                # Find GC contents of the strands
                # Count number of C
                num_of_c = DNAstrand.count('C')

                # Count number of G
                num_of_g = DNAstrand.count('G')

                # Get the total of G and C
                GCcontent = num_of_c + num_of_g

                # Append it to the GCcontents list
                GCcontents.append(GCcontent)

    # Get the number of total strands obtained
    totalStrands = len(lengthsOfSegments)

    # Get the maximum length of the strand
    maximum = max(lengthsOfSegments)

    # Get the minimum length of the strand
    minimum = min(lengthsOfSegments)

    # Get the average length of the strands
    avg = sum(lengthsOfSegments) / len(lengthsOfSegments)

    # Get average GC content
    avg_gc = (sum(GCcontents) / len(GCcontents)) / avg

    # plot graph using matplotlib
    # x-axis: Lengths of strands
    # y-axis: Frequency of strands occurance
    plt.hist(lengthsOfSegments)
    plt.xlabel('Lengths of strands')
    plt.ylabel('Frequency of strands')
    plt.title('Distribution of the Lengths of DNA Fragments')

    # Display the statistics
    print(f'Number of DNA fragments found:{totalStrands}')
    print(f'Average GC Content of the segments in percentage:{avg_gc * 100}%', )
    print(f'Maximum length of the DNA fragments:{maximum}')
    print(f'Min length of the DNA fragments:{minimum}')
    print(f'Average length of the DNA fragments:{avg}')
    plt.show()

    # return nothing on this function
    pass
