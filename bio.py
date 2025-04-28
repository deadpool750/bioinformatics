import tkinter as tk
import numpy as np
from Bio import SeqIO

def maxvalue_needleman(matrix, x_pos, y_pos, match_score, mismatch_penalty, gap_penalty):
    """
    Calculate and set maximum score for a given cell in the matrix
    based on match, mismatch, and gap values provided.
    """

    #first we are checking the diagonal we have a condition based on the matching or mismatching parts of the sequence
    if matrix[x_pos][0] == matrix[0][y_pos]: #match
        diagonal = matrix[x_pos - 1][y_pos - 1] + match_score
    else:
        diagonal = matrix[x_pos - 1][y_pos - 1] + mismatch_penalty #mismatch

    #we calculate the left value
    left = matrix[x_pos - 1][y_pos] + gap_penalty

    #we calculate the top value
    top = matrix[x_pos][y_pos - 1] + gap_penalty

    #the final value set in the matrix is the maximum value out of the three
    matrix[x_pos][y_pos] = max(left, top, diagonal)

def match_or_mismatch(matrix, x_pos, y_pos, match_score, mismatch_penalty):
    """
    Returns match or mismatch value based on the nucleotides provided in the 0th row and column.
    """

    #if the values are the same return match score else return mismatch penalty
    return match_score if matrix[x_pos][0] == matrix[0][y_pos] else mismatch_penalty

def traceback_needleman(matrix, match_score, mismatch_penalty, gap_penalty):
    """
    Traces back the path through the matrix and returns the path values
    and their coordinates.
    """

    #initializing lists for the coordinates and the values
    path_list = []
    path_coordinates = []

    #getting the length and width of the matrix for iteratign through it
    max_length = matrix.shape[0]
    max_width = matrix.shape[1]

    #subtracting 1 from it since we are operating on matrices
    position_x = max_length - 1
    position_y = max_width - 1

    #while loop that stops after we reach coordinates x < 1 or y < 1
    while position_x > 1 or position_y > 1:

        #adding the results to the list
        path_list.append(matrix[position_x, position_y])
        path_coordinates.append((position_x, position_y))

        #first we check the diagonal
        if matrix[position_x][position_y] == matrix[position_x - 1][position_y - 1] + match_or_mismatch(
                matrix, position_x, position_y, match_score, mismatch_penalty):
            position_x -= 1
            position_y -= 1

        #if not the diagonal then top value
        elif matrix[position_x][position_y] == matrix[position_x - 1][position_y] + gap_penalty:
            position_x -= 1

        #if not diagonal and top then left value
        else:
            position_y -= 1

    #adding the [1,1] and [0,0] coordinates and values manually because they are the same everytime
    path_list.append(matrix[1, 1])
    path_coordinates.append((1, 1))

    path_list.append(matrix[0, 0])
    path_coordinates.append((0, 0))

    #coordinates given in the
    return path_list, path_coordinates

def generate_array(strand1, strand2, match_value, mismatch_value, gap_value):
    """
    Generates the Needleman-Wunsch alignment matrix based on the input strands and scoring values.
    """

    #creating appropriate length and width for the matrix taking into account the space for gap penalties and 0 values
    row_length = len(strand1) + 2
    column_length = len(strand2) + 2

    #validating the strand nucleotides
    bases = {'A', 'C', 'G', 'T', 'U'}
    if not set(strand1).issubset(bases) or not set(strand2).issubset(bases):
        print("Invalid bases detected!")
        return None

    #creating the array originally all with zeros
    matrix = np.zeros((row_length, column_length), dtype='object')

    #placing the two strands into the matrix
    matrix[0, 2:2 + len(strand2)] = list(strand2)
    matrix[2:2 + len(strand1), 0] = list(strand1)

    #placing the zero manually as it will be always in the same place
    matrix[1, 1] = 0

    #placing the gap values in every cell of the 0th row
    n = gap_value
    for i in range(2, row_length):
        matrix[i, 1] = n
        n += gap_value

    #placing the gap values in every cell of the 0th column
    n = gap_value
    for i in range(2, column_length):
        matrix[1, i] = n
        n += gap_value

    #filling in the values using previously created maxvalue_needleman function
    for i in range(2, row_length):
        for j in range(2, column_length):
            maxvalue_needleman(matrix, i, j, match_value, mismatch_value, gap_value)

    return matrix

def get_from_fasta_file(file_name):
    """
    Reads two sequences from a fasta file with the name provided and returns them as a tuple.
    (I am using Biopython library here)
    """

    #initializing the list for the strands
    sequences = []

    #extracting 2 strands from the file
    for record in SeqIO.parse(file_name, "fasta"):
        sequences.append(str(record.seq))

    #checking if there are exactly 2 values in the file
    if len(sequences) != 2:
        raise ValueError("Two values in the fasta file expected!")

    return sequences[0], sequences[1]

def analyze_alignment(matrix, match_score, mismatch_penalty, gap_penalty, strand1, strand2):
    #initializing lists to store the aligned sequences
    aligned1 = []
    aligned2 = []

    #setting starting positions at the bottom-right of the matrix
    i = len(strand1) + 1
    j = len(strand2) + 1

    #tracing back through the matrix until reaching the top-left corner
    while i > 1 or j > 1:
        #checking for a diagonal move (match or mismatch)
        if i > 1 and j > 1:
            match = match_score if strand1[i - 2] == strand2[j - 2] else mismatch_penalty
            if matrix[i][j] == matrix[i - 1][j - 1] + match:
                aligned1.append(strand1[i - 2])
                aligned2.append(strand2[j - 2])
                i -= 1
                j -= 1
                continue

        #checking for a move from the top (gap in strand2)
        if i > 1 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            aligned1.append(strand1[i - 2])
            aligned2.append("-")
            i -= 1
            continue

        #checking for a move from the left (gap in strand1)
        if j > 1 and matrix[i][j] == matrix[i][j - 1] + gap_penalty:
            aligned1.append("-")
            aligned2.append(strand2[j - 2])
            j -= 1
            continue

    #reversing the sequences to get the final alignment
    aligned1.reverse()
    aligned2.reverse()

    #returning the aligned sequences
    return aligned1, aligned2


def save_alignment_to_file(filename, aligned1, aligned2, match, mismatch, gap, alignment_length, match_count, gap_count, identity_percentage):
    """
    saves all relevant alignment results to a text file
    """

    #structure how it should be saved into a .txt file
    with open(filename, "w") as f:
        f.write("Needleman-Wunsch Alignment Result\n")
        f.write("===============================\n\n")
        f.write("Aligned Sequences:\n")
        f.write("Strand 1: " + ''.join(aligned1) + "\n")
        f.write("Strand 2: " + ''.join(aligned2) + "\n\n")
        f.write(f"Scoring:\n  Match = {match}, Mismatch = {mismatch}, Gap = {gap}\n\n")
        f.write(f"Alignment Length: {alignment_length}\n")
        f.write(f"Number of Matches: {match_count}\n")
        f.write(f"Number of Gaps: {gap_count}\n")
        f.write(f"Identity Percentage: {identity_percentage}%\n")


class AlignmentGUI:
    """
    GUI class for visualizing the algorithm
    """
    def __init__(self, root):
        #initialize the window and default values
        self.root = root
        self.root.title("Welsch-Needleman Algorithm")

        self.matrix = None
        self.path = None
        self.path_coordinates = None

        #creating all GUI components
        self.create_input_fields()
        self.create_matrix_display()
        self.create_path_display()
        self.create_result_display()

    def create_input_fields(self):
        """
        Creates fields and buttons for the user.
        """
        #creating input for first strand
        tk.Label(self.root, text="Strand 1:").grid(row=0, column=0)
        self.strand1_entry = tk.Entry(self.root, width=10)
        self.strand1_entry.grid(row=0, column=1)

        #creating input for second strand
        tk.Label(self.root, text="Strand 2:").grid(row=0, column=2)
        self.strand2_entry = tk.Entry(self.root, width=10)
        self.strand2_entry.grid(row=0, column=3)

        #match value input
        tk.Label(self.root, text="Match:").grid(row=1, column=0)
        self.match_entry = tk.Entry(self.root, width=5)
        self.match_entry.insert(0, "1")
        self.match_entry.grid(row=1, column=1)

        #mismatch value input
        tk.Label(self.root, text="Mismatch:").grid(row=1, column=2)
        self.mismatch_entry = tk.Entry(self.root, width=5)
        self.mismatch_entry.insert(0, "-1")
        self.mismatch_entry.grid(row=1, column=3)

        #gap value input
        tk.Label(self.root, text="Gap:").grid(row=1, column=4)
        self.gap_entry = tk.Entry(self.root, width=5)
        self.gap_entry.insert(0, "-1")
        self.gap_entry.grid(row=1, column=5)

        #run button to start the algorithm
        self.run_button = tk.Button(self.root, text="Run", command=self.update_matrix)
        self.run_button.grid(row=1, column=6, padx=5, pady=5)

        #fasta file name input and button to run the program
        tk.Label(self.root, text="FASTA file:").grid(row=2, column=0)
        self.fasta_entry = tk.Entry(self.root, width=20)
        self.fasta_entry.grid(row=2, column=1, columnspan=2)

        #get from fasta button
        self.fasta_button = tk.Button(self.root, text="Get from FASTA", command=self.load_from_fasta)
        self.fasta_button.grid(row=2, column=3, padx=5, pady=5)

    def create_matrix_display(self):
        """
        Creates a frame for showing the matrix.
        """
        #frame for the matrix table
        self.matrix_frame = tk.Frame(self.root)
        self.matrix_frame.grid(row=3, column=0, columnspan=7)

    def create_path_display(self):
        """
        Creates a frame for showing the path
        """
        #frame for showing coordinates of the optimal path
        self.path_frame = tk.Frame(self.root)
        self.path_frame.grid(row=4, column=0, columnspan=7)

    def display_matrix(self):
        """
        displays the matrix and highlights the path
        """

        #removing any previously displayed matrix from the frame
        for widget in self.matrix_frame.winfo_children():
            widget.destroy()

        #exiting if no matrix is available to display
        if self.matrix is None:
            return

        #iterating through each cell in the matrix to display it
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i])):

                #highlighting the cell if it is part of the optimal path, excluding the top-left cell
                if (i, j) in self.path_coordinates and (i, j) != (0, 0):
                    label = tk.Label(self.matrix_frame, text=str(self.matrix[i][j]), borderwidth=1, relief="solid",
                                     width=5, height=2, bg="green")
                #displaying a normal cell if it is not part of the path
                else:
                    label = tk.Label(self.matrix_frame, text=str(self.matrix[i][j]), borderwidth=1, relief="solid",
                                     width=5, height=2)

                #placing the label widget in the appropriate position in the grid
                label.grid(row=i, column=j, padx=2, pady=2)

    def display_path(self):
        """
        Displays the coordinates of the path under the drawn matrix.
        """
        #clearing previous path display
        for widget in self.path_frame.winfo_children():
            widget.destroy()

        if self.path is None or self.path_coordinates is None:
            return

        #building a string of coordinates
        path_str = '->'.join([f"({p[0]},{p[1]})" for p in self.path_coordinates])
        path_label = tk.Label(self.path_frame, text=f"Path (row, column): {path_str}")
        path_label.grid(row=0, column=0, padx=5, pady=5)

    def update_matrix(self):
        """
        updates the matrix, path, and alignment details after clicking run button
        """
        try:
            #retrieve input values from the user interface
            strand1 = self.strand1_entry.get().upper()
            strand2 = self.strand2_entry.get().upper()
            match = int(self.match_entry.get())
            mismatch = int(self.mismatch_entry.get())
            gap = int(self.gap_entry.get())

            #generate the scoring matrix based on the input sequences and scoring values
            self.matrix = generate_array(strand1, strand2, match, mismatch, gap)

            #trace back through the matrix to find the optimal alignment path
            self.path, self.path_coordinates = traceback_needleman(self.matrix, match, mismatch, gap)

            #display the updated matrix and path in the GUI
            self.display_matrix()
            self.display_path()

            #perform the alignment based on the calculated path and scoring
            aligned1, aligned2 = analyze_alignment(self.matrix, match, mismatch, gap, strand1, strand2)

            #calculate alignment metrics
            align_len = len(aligned1)
            match_count = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != '-' and b != '-')
            gap_count = aligned1.count('-') + aligned2.count('-')
            identity = round((match_count / align_len) * 100, 2)

            #update the GUI labels with the alignment results
            self.align_len_label.config(text=f"Alignment Length: {align_len}")
            self.match_count_label.config(text=f"Number of Matches: {match_count}")
            self.gap_count_label.config(text=f"Number of Gaps: {gap_count}")
            self.identity_label.config(text=f"Identity Percentage: {identity}%")

        except ValueError:
            #handle invalid input values
            print("Please enter valid values for match, mismatch, and gap.")

    def load_from_fasta(self):
        """
        Function for getting the sequences from a fasta file and updating the fields.
        """
        try:
            #getting fasta filename
            file_name = self.fasta_entry.get()

            #reading sequences from file
            strand1, strand2 = get_from_fasta_file(file_name)

            #updating input fields with loaded sequences
            self.strand1_entry.delete(0, tk.END)
            self.strand1_entry.insert(0, strand1)

            self.strand2_entry.delete(0, tk.END)
            self.strand2_entry.insert(0, strand2)

            #running update if scoring values already provided
            if all([self.match_entry.get(), self.mismatch_entry.get(), self.gap_entry.get()]):
                self.update_matrix()

        except Exception as e:
            print(f"Error reading FASTA file: {e}")

    def save_result(self):
        if self.matrix is None or self.path_coordinates is None:
            return

        strand1 = self.strand1_entry.get().upper()
        strand2 = self.strand2_entry.get().upper()
        match = int(self.match_entry.get())
        mismatch = int(self.mismatch_entry.get())
        gap = int(self.gap_entry.get())

        aligned1, aligned2 = analyze_alignment(self.matrix, match, mismatch, gap, strand1, strand2)

        if not aligned1 or not aligned2:
            return

        align_len = len(aligned1)
        match_count = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != '-' and b != '-')
        gap_count = aligned1.count('-') + aligned2.count('-')
        identity = round((match_count / align_len) * 100, 2)

        # Update GUI labels
        self.align_len_label.config(text=f"Alignment Length: {align_len}")
        self.match_count_label.config(text=f"Number of Matches: {match_count}")
        self.gap_count_label.config(text=f"Number of Gaps: {gap_count}")
        self.identity_label.config(text=f"Identity Percentage: {identity}%")

        # Save to file
        save_alignment_to_file("alignment_result.txt", aligned1, aligned2, match, mismatch, gap, align_len, match_count,
                               gap_count, identity)

    def create_result_display(self):
        """
        creates labels to display alignment results and adds the save button below them
        """

        #create a frame for the result labels and the save button
        self.result_frame = tk.Frame(self.root)
        self.result_frame.grid(row=6, column=0, columnspan=7, pady=5)

        #create empty labels for alignment results, updated after clicking run
        self.align_len_label = tk.Label(self.result_frame, text="")
        self.align_len_label.pack()

        self.match_count_label = tk.Label(self.result_frame, text="")
        self.match_count_label.pack()

        self.gap_count_label = tk.Label(self.result_frame, text="")
        self.gap_count_label.pack()

        self.identity_label = tk.Label(self.result_frame, text="")
        self.identity_label.pack()

        #create the save button
        self.save_button = tk.Button(self.result_frame, text="Save Result", command=self.save_result)
        self.save_button.pack(pady=5)


#launching the GUI
root = tk.Tk()
app = AlignmentGUI(root)
root.mainloop()