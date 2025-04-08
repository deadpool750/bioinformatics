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
        #clearing previous matrix display
        for widget in self.matrix_frame.winfo_children():
            widget.destroy()

        if self.matrix is None:
            return

        #drawing each cell of the matrix
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i])):
                #highlighting cells that are in the path in green
                if (i, j) in self.path_coordinates:
                    label = tk.Label(self.matrix_frame, text=str(self.matrix[i][j]), borderwidth=1, relief="solid",
                                     width=5, height=2, bg="green")
                else:
                    label = tk.Label(self.matrix_frame, text=str(self.matrix[i][j]), borderwidth=1, relief="solid",
                                     width=5, height=2)
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
        path_label = tk.Label(self.path_frame, text=f"Path: {path_str}")
        path_label.grid(row=0, column=0, padx=5, pady=5)

    def update_matrix(self):
        """
        Updates the matrix and path after clicking run button.
        """
        try:
            #getting all user inputs
            strand1 = self.strand1_entry.get().upper()
            strand2 = self.strand2_entry.get().upper()
            match = int(self.match_entry.get())
            mismatch = int(self.mismatch_entry.get())
            gap = int(self.gap_entry.get())

            print(f"Running algorithm with: Strand1={strand1}, Strand2={strand2}, Match={match}, Mismatch={mismatch}, Gap={gap}")

            #generating matrix and path
            self.matrix = generate_array(strand1, strand2, match, mismatch, gap)
            self.path, self.path_coordinates = traceback_needleman(self.matrix, match, mismatch, gap)

            #displaying results
            self.display_matrix()
            self.display_path()

        except ValueError:
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

#launching the GUI
root = tk.Tk()
app = AlignmentGUI(root)
root.mainloop()