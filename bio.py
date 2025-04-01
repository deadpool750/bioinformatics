import tkinter as tk
import numpy as np


def maxvalue_needleman(matrix, x_pos, y_pos, match_score, mismatch_penalty, gap_penalty):
    if matrix[x_pos][0] == matrix[0][y_pos]:
        diagonal = matrix[x_pos - 1][y_pos - 1] + match_score
    else:
        diagonal = matrix[x_pos - 1][y_pos - 1] + mismatch_penalty

    left = matrix[x_pos - 1][y_pos] + gap_penalty
    top = matrix[x_pos][y_pos - 1] + gap_penalty

    matrix[x_pos][y_pos] = max(left, top, diagonal)

def match_or_mismatch(matrix, x_pos, y_pos, match_score, mismatch_penalty):
    value = 0
    if matrix[x_pos][0] == matrix[0][y_pos]:
        value = match_score
    else:
        value = mismatch_penalty

    return value


def traceback_needleman(matrix, match_score, mismatch_penalty, gap_penalty):
    path_list = []

    max_length = matrix.shape[0]
    max_width = matrix.shape[1]

    position_x = max_length
    position_y = max_width
    if position_x == 1 and position_y == 1 and matrix[position_x][position_y] == 0:
        return path_list
    else:
        # moving diagonal
        if matrix[position_x][position_y] == matrix[position_x - 1][position_y - 1] + match_or_mismatch(matrix,
                                                                                                        position_x,
                                                                                                        position_y,
                                                                                                        match_score,
                                                                                                        mismatch_penalty):
            path_list.append([position_x - 1, position_y - 1])
            position_x -= 1
            position_y -= 1
        # moving up
        elif matrix[position_x][position_y] == matrix[position_x - 1][position_y] + gap_penalty:
            path_list.append([position_x - 1, position_y])
            position_x -=  1
        # moving left
        else:
            path_list.append([position_x, position_y - 1])
            position_y -=  1

def generate_array(strand1, strand2, match_value, mismatch_value, gap_value):
    row_length = len(strand1) + 2
    column_length = len(strand2) + 2

    bases = {'A', 'C', 'G', 'T'}
    if not set(strand1).issubset(bases) or not set(strand2).issubset(bases):
        print("Invalid bases detected!")
        return None

    matrix = np.zeros((row_length, column_length), dtype='object')

    # Fill headers
    matrix[0, 2:2 + len(strand2)] = list(strand2)
    matrix[2:2 + len(strand1), 0] = list(strand1)

    matrix[1, 1] = 0

    # Fill gap penalties
    n = gap_value
    for i in range(2, row_length):
        matrix[i, 1] = n
        n += gap_value

    n = gap_value
    for i in range(2, column_length):
        matrix[1, i] = n
        n += gap_value

    # Compute matrix values
    for i in range(2, row_length):
        for j in range(2, column_length):
            maxvalue_needleman(matrix, i, j, match_value, mismatch_value, gap_value)

    return matrix


class AlignmentGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Welsch-Needleman Alignment")

        self.matrix = None

        # Create UI Elements
        self.create_input_fields()
        self.create_matrix_display()

    def create_input_fields(self):
        """Creates entry fields for strands and scoring values."""
        tk.Label(self.root, text="Strand 1:").grid(row=0, column=0)
        self.strand1_entry = tk.Entry(self.root, width=10)
        self.strand1_entry.grid(row=0, column=1)

        tk.Label(self.root, text="Strand 2:").grid(row=0, column=2)
        self.strand2_entry = tk.Entry(self.root, width=10)
        self.strand2_entry.grid(row=0, column=3)

        tk.Label(self.root, text="Match:").grid(row=1, column=0)
        self.match_entry = tk.Entry(self.root, width=5)
        self.match_entry.grid(row=1, column=1)

        tk.Label(self.root, text="Mismatch:").grid(row=1, column=2)
        self.mismatch_entry = tk.Entry(self.root, width=5)
        self.mismatch_entry.grid(row=1, column=3)

        tk.Label(self.root, text="Gap:").grid(row=1, column=4)
        self.gap_entry = tk.Entry(self.root, width=5)
        self.gap_entry.grid(row=1, column=5)

        self.run_button = tk.Button(self.root, text="Run", command=self.update_matrix)
        self.run_button.grid(row=1, column=6, padx=5, pady=5)

    def create_matrix_display(self):
        """Placeholder frame for matrix display."""
        self.matrix_frame = tk.Frame(self.root)
        self.matrix_frame.grid(row=2, column=0, columnspan=7)

    def display_matrix(self):
        """Displays the computed matrix."""
        # Clear previous matrix
        for widget in self.matrix_frame.winfo_children():
            widget.destroy()

        if self.matrix is None:
            return

        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[i])):
                label = tk.Label(self.matrix_frame, text=str(self.matrix[i][j]), borderwidth=1, relief="solid", width=5, height=2)
                label.grid(row=i, column=j, padx=2, pady=2)

    def update_matrix(self):
        try:
            strand1 = self.strand1_entry.get().upper()
            strand2 = self.strand2_entry.get().upper()

            match = int(self.match_entry.get())
            mismatch = int(self.mismatch_entry.get())
            gap = int(self.gap_entry.get())

            print(f"Running algorithm with: Strand1={strand1}, Strand2={strand2}, Match={match}, Mismatch={mismatch}, Gap={gap}")

            self.matrix = generate_array(strand1, strand2, match, mismatch, gap)
            self.display_matrix()

        except ValueError:
            print("Please enter valid numbers for match, mismatch, and gap.")


root = tk.Tk()
app = AlignmentGUI(root)
root.mainloop()
